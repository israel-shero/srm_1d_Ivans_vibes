#!/usr/bin/env python3
"""Run manifested Chunc limiting-case checks without changing defaults."""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from cases.loader import PROJECT_ROOT
from scripts.run_numerical_verification import _extract_metrics
from srm_1d.burn_rate import burn_rate_cell, saint_robert_from_tabs
from srm_1d.openmotor_adapter import convert_propellant, load_ric, run_from_ric
from srm_1d.run_artifacts import artifact_dir, verify_run_health


CONFIG_PATH = PROJECT_ROOT / "cases" / "baseline_configs.json"
STEADY_INTERVAL_S = (0.6, 1.8)
STEADY_PRESSURE_CV_LIMIT = 0.05
STEADY_PRESSURE_SLOPE_LIMIT_PERCENT_PER_S = 5.0


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_value(*args: str) -> str:
    completed = subprocess.run(
        ["git", *args], cwd=PROJECT_ROOT, check=False,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else ""


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _steady_interval_metrics(
    result: dict,
    start_s: float = STEADY_INTERVAL_S[0],
    end_s: float = STEADY_INTERVAL_S[1],
) -> dict:
    time = np.asarray(result["time"], dtype=float)
    pressure = np.asarray(result["P_head"], dtype=float)
    mask = (time >= start_s) & (time <= end_s)
    if np.count_nonzero(mask) < 2:
        return {
            "status": "not_evaluated",
            "reason": "fewer than two history points in declared interval",
            "start_s": start_s,
            "end_s": end_s,
        }

    interval_time = time[mask]
    interval_pressure = pressure[mask]
    mean_pressure = float(np.mean(interval_pressure))
    centered_time = interval_time - np.mean(interval_time)
    denominator = float(np.dot(centered_time, centered_time))
    slope = 0.0
    if denominator > 0.0:
        slope = float(
            np.dot(centered_time, interval_pressure - mean_pressure)
            / denominator
        )
    cv = float(np.std(interval_pressure) / mean_pressure)
    slope_percent_per_s = float(slope / mean_pressure * 100.0)
    screen_passed = (
        cv <= STEADY_PRESSURE_CV_LIMIT
        and abs(slope_percent_per_s)
        <= STEADY_PRESSURE_SLOPE_LIMIT_PERCENT_PER_S
    )
    return {
        "status": "evaluated",
        "start_s": start_s,
        "end_s": end_s,
        "samples": int(np.count_nonzero(mask)),
        "mean_pressure_mpa": mean_pressure / 1.0e6,
        "pressure_cv": cv,
        "pressure_slope_mpa_per_s": slope / 1.0e6,
        "pressure_slope_percent_per_s": slope_percent_per_s,
        "pressure_range_percent": float(
            (np.max(interval_pressure) - np.min(interval_pressure))
            / mean_pressure * 100.0
        ),
        "screen": {
            "pressure_cv_limit": STEADY_PRESSURE_CV_LIMIT,
            "absolute_pressure_slope_limit_percent_per_s": (
                STEADY_PRESSURE_SLOPE_LIMIT_PERCENT_PER_S
            ),
            "passed": bool(screen_passed),
            "qualification": (
                "Numerical stationarity screen, not experimental steady-state "
                "validation."
            ),
        },
    }


def _relative_change(reference: float, candidate: float) -> float | None:
    if reference == 0.0:
        return None
    return (candidate - reference) / abs(reference) * 100.0


def _closure_limiting_cases(propellant) -> dict:
    tab = propellant.representative_tab()
    tab_min, tab_max, tab_a, tab_n = propellant.tab_arrays()
    pressure = min(max(5.0e6, tab.min_pressure + 1.0), tab.max_pressure - 1.0)
    r0 = saint_robert_from_tabs(
        pressure, tab_min, tab_max, tab_a, tab_n, len(tab_a),
    )
    common = (
        pressure, 0.04, 0.5, 35.0e-6,
        propellant.Cp_gas * propellant.mu_gas / propellant.k_gas,
        propellant.k_gas, propellant.Cp_gas,
        tab.T_flame, propellant.T_surface,
        propellant.rho_propellant, propellant.Cps, propellant.T_initial,
        tab_min, tab_max, tab_a, tab_n, len(tab_a), 0.45,
    )
    zero_total, zero_erosive = burn_rate_cell(
        common[0], 0.0, *common[1:]
    )
    flow_total, flow_erosive = burn_rate_cell(
        common[0], 5.0e5, *common[1:]
    )
    normal_from_flow = flow_total - flow_erosive
    return {
        "pressure_pa": pressure,
        "normal_rate_m_s": float(r0),
        "zero_crossflow": {
            "total_rate_m_s": float(zero_total),
            "erosive_rate_m_s": float(zero_erosive),
            "normal_rate_preserved": bool(
                np.isclose(zero_total, r0, rtol=1.0e-12, atol=1.0e-15)
                and zero_erosive == 0.0
            ),
        },
        "finite_crossflow": {
            "reynolds_number": 5.0e5,
            "total_rate_m_s": float(flow_total),
            "erosive_rate_m_s": float(flow_erosive),
            "recovered_normal_rate_m_s": float(normal_from_flow),
            "normal_rate_preserved": bool(
                np.isclose(normal_from_flow, r0, rtol=1.0e-12, atol=1.0e-15)
            ),
        },
    }


def run_limiting_cases(
    output_root: Path | None = None,
    resume_path: Path | None = None,
) -> Path:
    baseline = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))["cases"][
        "chunc"
    ]
    motor_path = PROJECT_ROOT / baseline["motor_path"]
    options = dict(baseline["run_options"])
    options.update({
        "target_propellant_cells": 100,
        "cfl_target": 0.3,
        "print_interval": 0.0,
        "snapshot_interval": 0.1,
    })
    if resume_path is not None:
        result_path = resume_path.resolve()
        output = result_path.parent
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        provenance = payload.get("provenance", {})
        if provenance.get("motor_sha256") != _sha256(motor_path):
            raise ValueError("resume motor checksum does not match current input")
        if provenance.get("baseline_config_sha256") != _sha256(CONFIG_PATH):
            raise ValueError(
                "resume baseline-config checksum does not match current input"
            )
        payload["status"] = "running"
        payload["resumed_utc"] = datetime.now(timezone.utc).isoformat()
    else:
        output = artifact_dir(
            "verification_limiting_cases",
            root=output_root if output_root is not None else PROJECT_ROOT,
        )
        result_path = output / "limiting_cases.json"
        payload = {
            "schema_version": 1,
            "status": "running",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "study": "chunc_limiting_cases",
            "configuration": options,
            "steady_interval_s": list(STEADY_INTERVAL_S),
            "provenance": {
                "python": sys.version,
                "platform": platform.platform(),
                "git_commit": _git_value("rev-parse", "HEAD"),
                "git_status": _git_value("status", "--short"),
                "motor_path": str(motor_path.relative_to(PROJECT_ROOT)),
                "motor_sha256": _sha256(motor_path),
                "baseline_config_sha256": _sha256(CONFIG_PATH),
            },
            "claims": {
                "research_note_results_reproduced": False,
                "limiting_cases_validate_full_physics": False,
            },
            "runs": {},
        }
    _write_json(result_path, payload)

    variants = (
        ("baseline", {}),
        ("erosive_disabled", {"diagnostic_disable_erosive": True}),
    )
    propellant = None
    for name, overrides in variants:
        if name in payload["runs"]:
            continue
        result, performance, _, _, run_propellant = run_from_ric(
            str(motor_path), verbose=False, **options, **overrides,
        )
        verify_run_health(
            result, motor_name=f"chunc-{name}", min_t_burn_s=1.9,
            raise_on_fail=True,
        )
        if propellant is None:
            propellant = run_propellant
        payload["runs"][name] = {
            "overrides": overrides,
            "metrics": _extract_metrics(result, performance),
            "steady_interval": _steady_interval_metrics(result),
        }
        _write_json(result_path, payload)
        del result, performance
        gc.collect()

    if propellant is None:
        propellant = convert_propellant(load_ric(str(motor_path))["propellant"])
    payload["closure_limiting_cases"] = _closure_limiting_cases(propellant)
    reference = payload["runs"]["baseline"]["metrics"]
    disabled = payload["runs"]["erosive_disabled"]["metrics"]
    payload["erosive_disabled_comparison"] = {
        key: _relative_change(reference[key], disabled[key])
        for key in (
            "peak_pressure_mpa", "peak_time_s", "total_impulse_ns",
            "performance_burn_time_s",
        )
    }
    payload["status"] = "complete"
    payload["completed_utc"] = datetime.now(timezone.utc).isoformat()
    _write_json(result_path, payload)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    output = run_limiting_cases(args.output_root, args.resume)
    print(f"Limiting-case artifacts: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
