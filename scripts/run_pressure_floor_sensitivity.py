#!/usr/bin/env python3
"""Run bounded Chunc pressure-floor sensitivity cases."""

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

from cases.loader import PROJECT_ROOT
from scripts.run_numerical_verification import _extract_metrics
from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.run_artifacts import artifact_dir, verify_run_health


CONFIG_PATH = PROJECT_ROOT / "cases" / "baseline_configs.json"
TARGET_CELLS = (800, 1600)
PRESSURE_FLOORS_PA = (100.0, 1000.0, 10000.0)
REFERENCE_PRESSURE_FLOOR_PA = 1000.0
SENSITIVITY_METRICS = (
    "peak_pressure_mpa",
    "peak_time_s",
    "max_fill_mach",
)
SENSITIVITY_SCREEN_PERCENT = 2.0
SOURCE_PATHS = (
    Path("scripts/run_pressure_floor_sensitivity.py"),
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
    "numpy",
    "scipy",
    "numba",
    "PyYAML",
    "matplotlib",
    "scikit-fmm",
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


def _build_run_matrix(
    target_cells: tuple[int, ...] = TARGET_CELLS,
    pressure_floors_pa: tuple[float, ...] = PRESSURE_FLOORS_PA,
) -> list[dict]:
    return [
        {
            "target_propellant_cells": int(cells),
            "pressure_floor_pa": float(floor_pa),
        }
        for cells in target_cells
        for floor_pa in pressure_floors_pa
    ]


def _relative_change_percent(value: float | None, reference: float | None):
    if value is None or reference is None or abs(float(reference)) == 0.0:
        return None
    return abs(float(value) - float(reference)) / abs(float(reference)) * 100.0


def _build_comparisons(rows: list[dict]) -> list[dict]:
    """Compare each floor to the same-grid 1000 Pa reference."""
    lookup = {
        (int(row["target_propellant_cells"]), float(row["pressure_floor_pa"])): row
        for row in rows
    }
    comparisons = []
    for config in _build_run_matrix():
        cells = config["target_propellant_cells"]
        floor_pa = config["pressure_floor_pa"]
        if floor_pa == REFERENCE_PRESSURE_FLOOR_PA:
            continue
        row = lookup.get((cells, floor_pa), {})
        reference = lookup.get((cells, REFERENCE_PRESSURE_FLOOR_PA), {})
        reference_ok = reference.get("status") == "complete"
        run_ok = row.get("status") == "complete"
        changes = {}
        for metric in SENSITIVITY_METRICS:
            value = row.get("metrics", {}).get(metric) if run_ok else None
            reference_value = (
                reference.get("metrics", {}).get(metric) if reference_ok else None
            )
            relative = _relative_change_percent(value, reference_value)
            changes[metric] = {
                "value": value,
                "reference_value": reference_value,
                "relative_change_percent": relative,
                "exceeds_screen": (
                    relative > SENSITIVITY_SCREEN_PERCENT
                    if relative is not None else None
                ),
            }
        comparisons.append({
            "target_propellant_cells": cells,
            "pressure_floor_pa": floor_pa,
            "reference_pressure_floor_pa": REFERENCE_PRESSURE_FLOOR_PA,
            "status": "compared" if run_ok and reference_ok else "unavailable",
            "metrics": changes,
        })
    return comparisons


def _summarize_run_statuses(rows: list[dict]) -> dict:
    counts = {}
    for row in rows:
        status = row["status"]
        counts[status] = counts.get(status, 0) + 1
    all_complete = len(rows) == len(_build_run_matrix()) and all(
        row["status"] == "complete" for row in rows
    )
    return {
        "all_runs_complete": all_complete,
        "counts": counts,
        "study_status": "complete" if all_complete else "complete_with_failures",
    }


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run_sensitivity(output_root: Path | None = None) -> Path:
    baseline = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))["cases"]["chunc"]
    motor_path = PROJECT_ROOT / baseline["motor_path"]
    openmotor_env = os.environ.get("SRM1D_OPENMOTOR_PATH")
    openmotor_path = Path(openmotor_env).resolve() if openmotor_env else None
    options = dict(baseline["run_options"])
    options.update({
        "t_max": 0.05,
        "cfl_target": 0.3,
        "snapshot_interval": 0.01,
        "print_interval": 0.0,
    })
    output = artifact_dir(
        "verification_pressure_floor_sensitivity",
        root=output_root if output_root is not None else PROJECT_ROOT,
    )
    result_path = output / "pressure_floor_sensitivity.json"
    payload = {
        "schema_version": 1,
        "status": "running",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "study": "chunc_pressure_floor_sensitivity",
        "scope": {
            "target_propellant_cells": list(TARGET_CELLS),
            "pressure_floor_pa": list(PRESSURE_FLOORS_PA),
            "reference_pressure_floor_pa": REFERENCE_PRESSURE_FLOOR_PA,
            "metrics": list(SENSITIVITY_METRICS),
            "sensitivity_screen_percent": SENSITIVITY_SCREEN_PERCENT,
            "qualification": (
                "Bounded numerical sensitivity screen only; it does not establish "
                "convergence or physical validation."
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
                "configured_path": (
                    str(openmotor_path) if openmotor_path is not None else None
                ),
                "git_commit": (
                    _git_value("rev-parse", "HEAD", cwd=openmotor_path)
                    if openmotor_path is not None else None
                ),
                "git_status": (
                    _git_value("status", "--short", cwd=openmotor_path)
                    if openmotor_path is not None else None
                ),
            },
        },
        "configuration": options,
        "claims": {
            "convergence_established": False,
            "physical_validation_established": False,
            "research_note_results_reproduced": False,
        },
        "runs": [],
        "comparisons": [],
    }
    _write_json(result_path, payload)

    for run_config in _build_run_matrix():
        row = {
            **run_config,
            "status": "running",
            "error": None,
            "metrics": None,
        }
        payload["runs"].append(row)
        _write_json(result_path, payload)
        try:
            result, performance, *_ = run_from_ric(
                str(motor_path), verbose=False,
                **dict(options, **run_config),
            )
            row["metrics"] = _extract_metrics(result, performance)
            row["health_check_passed"] = verify_run_health(
                result,
                motor_name=(
                    "chunc-pressure-floor-sensitivity "
                    f"cells={run_config['target_propellant_cells']} "
                    f"floor={run_config['pressure_floor_pa']} Pa"
                ),
                min_t_burn_s=0.049,
                raise_on_fail=False,
            )
            termination_code = int(result["summary"]["termination_code"])
            row["expected_termination_reached"] = termination_code == 0
            row["status"] = (
                "complete"
                if (row["health_check_passed"]
                    and row["expected_termination_reached"])
                else "unhealthy"
            )
            if not row["health_check_passed"]:
                row["error"] = "Run failed the existing health check."
            elif not row["expected_termination_reached"]:
                row["error"] = (
                    "Run did not reach the expected t_max termination "
                    f"(termination_code={termination_code})."
                )
        except Exception as exc:  # Continue the bounded matrix after a case failure.
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
    args = parser.parse_args()
    result_path = run_sensitivity(args.output_root)
    print(f"Pressure-floor sensitivity artifact: {result_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
