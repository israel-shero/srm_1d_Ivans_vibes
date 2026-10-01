#!/usr/bin/env python3
"""Reconcile historical and current Chunc startup-fill configurations."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from cases.loader import PROJECT_ROOT
from scripts.run_numerical_verification import _extract_metrics
from srm_1d.openmotor_adapter import load_pyrogen, run_from_ric
from srm_1d.run_artifacts import artifact_dir, verify_run_health


CONFIG_PATH = PROJECT_ROOT / "cases" / "baseline_configs.json"
HISTORICAL_PROBE_PATH = (
    PROJECT_ROOT / "docs" / "v0_7_4" / "probes" / "chunc_mach_convergence.py"
)
GRID_LEVELS = (50, 100, 200)
STARTUP_DURATION_S = 0.05
LEGACY_BPNV_MARKER = "bpnv_pre_bfc2f3f"
HISTORICAL_RECORDED = {
    50: {"max_fill_mach": 3.36, "peak_pressure_mpa": 11.77},
    100: {"max_fill_mach": 4.81, "peak_pressure_mpa": 12.65},
    200: {"max_fill_mach": 12.04, "peak_pressure_mpa": 13.22},
}


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


def _variants(current_options: dict) -> dict[str, dict]:
    """Return the predeclared factorial comparison configurations."""
    current = dict(current_options)
    current.update({
        "t_max": STARTUP_DURATION_S,
        "snapshot_interval": 0.01,
        "print_interval": 0.0,
    })
    historical_knobs = {
        "roughness": 32.0e-6,
        "kappa": 0.44,
        "T_ignition": 756.0,
        "P_cutoff": 0.01e6,
    }
    return {
        "current_mtv": dict(current),
        "current_bpnv": dict(current, pyrogen="bpnv"),
        "current_legacy_bpnv": dict(current, pyrogen=LEGACY_BPNV_MARKER),
        "historical_mtv": dict(current, pyrogen="mtv", **historical_knobs),
        "historical_bpnv": dict(current, pyrogen="bpnv", **historical_knobs),
        "historical_legacy_bpnv": dict(
            current, pyrogen=LEGACY_BPNV_MARKER, **historical_knobs
        ),
    }


def _resolve_options(recorded_options: dict) -> dict:
    """Resolve the serializable legacy marker to the historical gas state."""
    options = dict(recorded_options)
    if options.get("pyrogen") == LEGACY_BPNV_MARKER:
        options["pyrogen"] = replace(
            load_pyrogen("bpnv"),
            name="BPNV legacy pre-bfc2f3f",
            T_flame=2800.0,
            M=0.030,
            gamma=1.25,
            gas_mass_fraction=1.0,
        )
    return options


def _relative_change(reference: float, candidate: float) -> float | None:
    if reference == 0.0:
        return None
    return (candidate - reference) / abs(reference) * 100.0


def run_comparison(output_root: Path | None = None) -> Path:
    baseline = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))["cases"][
        "chunc"
    ]
    motor_path = PROJECT_ROOT / baseline["motor_path"]
    variants = _variants(baseline["run_options"])
    output = artifact_dir(
        "verification_fill_configuration",
        root=output_root if output_root is not None else PROJECT_ROOT,
    )
    result_path = output / "fill_configuration.json"
    payload = {
        "schema_version": 1,
        "status": "running",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "study": "chunc_historical_bpnv_vs_current_mtv_startup_fill",
        "scope": {
            "duration_s": STARTUP_DURATION_S,
            "fill_window_s": 0.003,
            "target_propellant_cells": list(GRID_LEVELS),
            "cfl_target": 0.3,
            "design": (
                "2x2 pyrogen-by-nonpyrogen-configuration comparison; "
                "historical values are transcribed from the named probe"
            ),
        },
        "provenance": {
            "python": sys.version,
            "platform": platform.platform(),
            "git_commit": _git_value("rev-parse", "HEAD"),
            "git_status": _git_value("status", "--short"),
            "motor_path": str(motor_path.relative_to(PROJECT_ROOT)),
            "motor_sha256": _sha256(motor_path),
            "baseline_config_sha256": _sha256(CONFIG_PATH),
            "historical_probe_path": str(
                HISTORICAL_PROBE_PATH.relative_to(PROJECT_ROOT)
            ),
            "historical_probe_sha256": _sha256(HISTORICAL_PROBE_PATH),
        },
        "claims": {
            "historical_results_reproduced": False,
            "research_note_results_reproduced": False,
            "physics_validated": False,
        },
        "historical_recorded_results": {
            "source": (
                "docs/v0_7_4/IGNITION_SPIKE_REOPENED.md section 6, added at "
                "commit 16bc527 before BPNV gas-property commit bfc2f3f"
            ),
            "qualification": (
                "Repository-reported values for comparison; exact reproduction "
                "is not claimed."
            ),
            "grid": HISTORICAL_RECORDED,
        },
        "variants": {},
    }
    result_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    for name, fixed_options in variants.items():
        rows = []
        for cells in GRID_LEVELS:
            options = dict(
                _resolve_options(fixed_options),
                target_propellant_cells=cells,
                cfl_target=0.3,
            )
            result, performance, *_ = run_from_ric(
                str(motor_path), verbose=False, **options,
            )
            verify_run_health(
                result,
                motor_name=f"chunc-{name}-cells={cells}",
                min_t_burn_s=0.049,
                raise_on_fail=True,
            )
            rows.append({
                "target_propellant_cells": cells,
                **_extract_metrics(result, performance),
            })
        payload["variants"][name] = {
            "fixed_options": fixed_options,
            "grid": rows,
        }
        result_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )

    comparisons = {}
    for cells in GRID_LEVELS:
        by_name = {
            name: next(
                row for row in item["grid"]
                if row["target_propellant_cells"] == cells
            )
            for name, item in payload["variants"].items()
        }
        reference = by_name["current_mtv"]
        comparisons[str(cells)] = {
            name: {
                metric: _relative_change(reference[metric], row[metric])
                for metric in ("max_fill_mach", "peak_pressure_mpa", "peak_time_s")
            }
            for name, row in by_name.items()
            if name != "current_mtv"
        }
    payload["comparisons_vs_current_mtv_percent"] = comparisons
    legacy_rows = payload["variants"]["historical_legacy_bpnv"]["grid"]
    payload["historical_recorded_differences"] = {
        str(row["target_propellant_cells"]): {
            metric: row[metric] - HISTORICAL_RECORDED[
                row["target_propellant_cells"]
            ][metric]
            for metric in ("max_fill_mach", "peak_pressure_mpa")
        }
        for row in legacy_rows
    }
    payload["status"] = "complete"
    payload["completed_utc"] = datetime.now(timezone.utc).isoformat()
    result_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return output


def main() -> int:
    output = run_comparison()
    print(f"Fill-configuration artifacts: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
