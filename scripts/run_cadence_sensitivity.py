#!/usr/bin/env python3
"""Run a bounded, resumable Chunc startup cadence sensitivity matrix."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

PROJECT_ROOT = SCRIPT_ROOT
CONFIG_PATH = PROJECT_ROOT / "cases" / "baseline_configs.json"
TARGET_CELLS = 800
CFL_TARGET = 0.075
T_MAX_S = 0.05
CADENCE_MATRIX = ((161, 161), (80, 161), (161, 80), (80, 80))
REFERENCE_CADENCE = (161, 161)
SENSITIVITY_SCREEN_PERCENT = 2.0
COMPARISON_METRICS = (
    "peak_pressure_mpa",
    "peak_time_s",
    "max_fill_mach",
    "minimum_pressure_pa",
    "integrated_abs_clipping_energy_j",
    "ignition_rate_refresh_delay_max_s",
)
SOURCE_PATHS = (
    Path("scripts/run_cadence_sensitivity.py"),
    Path("scripts/run_numerical_verification.py"),
    Path("srm_1d/simulation.py"),
    Path("srm_1d/solver.py"),
)
DEPENDENCY_MANIFESTS = (
    Path("pyproject.toml"),
    Path("requirements/workspace-lock.json"),
    Path("requirements/py310-macos-x86_64.lock.txt"),
)
RUNTIME_DISTRIBUTIONS = (
    "numpy", "scipy", "numba", "PyYAML", "matplotlib", "scikit-fmm",
    "scikit-image",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_value(*args: str, cwd: Path = PROJECT_ROOT) -> str:
    completed = subprocess.run(
        ["git", *args], cwd=cwd, check=False,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    return completed.stdout.strip() if completed.returncode == 0 else ""


def _runtime_versions() -> dict[str, str | None]:
    versions = {}
    for distribution in RUNTIME_DISTRIBUTIONS:
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = None
    return versions


def _build_run_matrix() -> list[dict]:
    return [
        {
            "target_propellant_cells": TARGET_CELLS,
            "burn_update_interval": burn_interval,
            "geometry_update_interval": geometry_interval,
        }
        for burn_interval, geometry_interval in CADENCE_MATRIX
    ]


def _relative_change_percent(value, reference):
    if value is None or reference is None or abs(float(reference)) == 0.0:
        return None
    return abs(float(value) - float(reference)) / abs(float(reference)) * 100.0


def _build_comparisons(rows: list[dict]) -> list[dict]:
    lookup = {
        (int(row["burn_update_interval"]), int(row["geometry_update_interval"])): row
        for row in rows
    }
    reference = lookup.get(REFERENCE_CADENCE, {})
    reference_ok = reference.get("status") == "complete"
    comparisons = []
    for burn_interval, geometry_interval in CADENCE_MATRIX:
        if (burn_interval, geometry_interval) == REFERENCE_CADENCE:
            continue
        row = lookup.get((burn_interval, geometry_interval), {})
        run_ok = row.get("status") == "complete"
        changes = {}
        for metric in COMPARISON_METRICS:
            value = row.get("metrics", {}).get(metric) if run_ok else None
            reference_value = (
                reference.get("metrics", {}).get(metric) if reference_ok else None
            )
            relative_change = _relative_change_percent(value, reference_value)
            changes[metric] = {
                "value": value,
                "reference_value": reference_value,
                "absolute_change": (
                    abs(float(value) - float(reference_value))
                    if value is not None and reference_value is not None else None
                ),
                "relative_change_percent": relative_change,
                "exceeds_screen": (
                    relative_change > SENSITIVITY_SCREEN_PERCENT
                    if relative_change is not None else None
                ),
            }
        comparisons.append({
            "burn_update_interval": burn_interval,
            "geometry_update_interval": geometry_interval,
            "reference_burn_update_interval": REFERENCE_CADENCE[0],
            "reference_geometry_update_interval": REFERENCE_CADENCE[1],
            "status": "compared" if run_ok and reference_ok else "unavailable",
            "metrics": changes,
        })
    return comparisons


def _summarize_run_statuses(rows: list[dict]) -> dict:
    counts = {}
    for row in rows:
        status = row["status"]
        counts[status] = counts.get(status, 0) + 1
    all_complete = len(rows) == len(CADENCE_MATRIX) and all(
        row["status"] == "complete" for row in rows
    )
    return {
        "all_runs_complete": all_complete,
        "counts": counts,
        "study_status": "complete" if all_complete else "complete_with_failures",
    }


def _expected_configuration(baseline: dict) -> dict:
    options = dict(baseline["run_options"])
    options.update({
        "target_propellant_cells": TARGET_CELLS,
        "t_max": T_MAX_S,
        "cfl_target": CFL_TARGET,
        "snapshot_interval": 0.01,
        "print_interval": 0.0,
    })
    return options


def _validate_resume(payload: dict, motor_path: Path, configuration: dict) -> None:
    provenance = payload.get("provenance", {})
    if provenance.get("motor_sha256") != _sha256(motor_path):
        raise ValueError("resume motor checksum does not match current input")
    if provenance.get("baseline_config_sha256") != _sha256(CONFIG_PATH):
        raise ValueError("resume baseline-config checksum does not match")
    if payload.get("configuration") != configuration:
        raise ValueError("resume configuration does not match current run options")
    if payload.get("scope", {}).get("cadence_matrix") != [
        list(pair) for pair in CADENCE_MATRIX
    ]:
        raise ValueError("resume cadence matrix does not match")
    expected_motor = str(motor_path.relative_to(PROJECT_ROOT))
    if provenance.get("motor_path") != expected_motor:
        raise ValueError("resume motor path does not match current input")


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _new_payload(motor_path: Path, configuration: dict) -> dict:
    openmotor_env = os.environ.get("SRM1D_OPENMOTOR_PATH")
    openmotor_path = Path(openmotor_env).resolve() if openmotor_env else None
    return {
        "schema_version": 1,
        "status": "running",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "study": "chunc_startup_cadence_sensitivity",
        "scope": {
            "target_propellant_cells": TARGET_CELLS,
            "cfl_target": CFL_TARGET,
            "t_max_s": T_MAX_S,
            "cadence_matrix": [list(pair) for pair in CADENCE_MATRIX],
            "reference_cadence": list(REFERENCE_CADENCE),
            "metrics": list(COMPARISON_METRICS),
            "sensitivity_screen_percent": SENSITIVITY_SCREEN_PERCENT,
            "qualification": (
                "Bounded startup cadence sensitivity screen; cadence effects are "
                "not evidence of numerical convergence or physical validation."
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
            "source_sha256": {
                str(path): _sha256(PROJECT_ROOT / path) for path in SOURCE_PATHS
            },
            "dependency_manifest_sha256": {
                str(path): _sha256(PROJECT_ROOT / path)
                for path in DEPENDENCY_MANIFESTS
            },
            "runtime_package_versions": _runtime_versions(),
            "openmotor": {
                "configured_path": str(openmotor_path) if openmotor_path else None,
                "git_commit": (
                    _git_value("rev-parse", "HEAD", cwd=openmotor_path)
                    if openmotor_path else None
                ),
                "git_status": (
                    _git_value("status", "--short", cwd=openmotor_path)
                    if openmotor_path else None
                ),
            },
        },
        "configuration": configuration,
        "claims": {
            "convergence_established": False,
            "physics_validated": False,
            "research_note_results_reproduced": False,
        },
        "runs": [
            {**row, "status": "pending", "error": None, "metrics": None}
            for row in _build_run_matrix()
        ],
        "comparisons": [],
    }


def run_sensitivity(
    output_root: Path | None = None,
    resume_path: Path | None = None,
) -> Path:
    from scripts.run_numerical_verification import _extract_metrics
    from srm_1d.openmotor_adapter import run_from_ric
    from srm_1d.run_artifacts import artifact_dir, verify_run_health

    baseline = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))["cases"]["chunc"]
    motor_path = PROJECT_ROOT / baseline["motor_path"]
    configuration = _expected_configuration(baseline)
    if resume_path is not None:
        result_path = resume_path.resolve()
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        _validate_resume(payload, motor_path, configuration)
        payload["status"] = "running"
        payload["resumed_utc"] = datetime.now(timezone.utc).isoformat()
    else:
        output = artifact_dir(
            "verification_chunc_cadence_sensitivity",
            root=output_root if output_root is not None else PROJECT_ROOT,
        )
        result_path = output / "cadence_sensitivity.json"
        payload = _new_payload(motor_path, configuration)
    _write_json(result_path, payload)

    for row in payload["runs"]:
        if row["status"] not in ("pending", "running"):
            continue
        row["status"] = "running"
        row["error"] = None
        _write_json(result_path, payload)
        run_config = {
            key: row[key] for key in (
                "target_propellant_cells", "burn_update_interval",
                "geometry_update_interval",
            )
        }
        try:
            result, performance, *_ = run_from_ric(
                str(motor_path), verbose=False,
                **dict(configuration, **run_config),
            )
            row["health_check_passed"] = verify_run_health(
                result,
                motor_name=(
                    "chunc-cadence-sensitivity "
                    f"burn={row['burn_update_interval']} "
                    f"geometry={row['geometry_update_interval']}"
                ),
                min_t_burn_s=T_MAX_S - 0.001,
                raise_on_fail=False,
            )
            row["metrics"] = _extract_metrics(result, performance)
            termination_code = int(result["summary"]["termination_code"])
            row["expected_termination_reached"] = termination_code == 0
            row["pressure_floor_activated"] = (
                row["metrics"].get("limit_pressure_floor_total_activations", 0) > 0
            )
            row["status"] = (
                "complete"
                if row["health_check_passed"] and row["expected_termination_reached"]
                else "unhealthy"
            )
            if not row["health_check_passed"]:
                row["error"] = "Run failed the existing health check."
            elif not row["expected_termination_reached"]:
                row["error"] = (
                    "Run did not reach the expected t_max termination "
                    f"(termination_code={termination_code})."
                )
        except Exception as exc:  # Preserve the rest of the bounded matrix.
            row["status"] = "failed"
            row["error"] = f"{type(exc).__name__}: {exc}"
        payload["comparisons"] = _build_comparisons(payload["runs"])
        _write_json(result_path, payload)

    payload["comparisons"] = _build_comparisons(payload["runs"])
    outcome = _summarize_run_statuses(payload["runs"])
    payload["status"] = outcome["study_status"]
    payload["run_outcomes"] = {
        "all_runs_complete": outcome["all_runs_complete"],
        "counts": outcome["counts"],
    }
    payload["completed_utc"] = datetime.now(timezone.utc).isoformat()
    _write_json(result_path, payload)
    return result_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    result_path = run_sensitivity(args.output_root, args.resume)
    print(f"Cadence-sensitivity artifact: {result_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
