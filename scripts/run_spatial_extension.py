#!/usr/bin/env python3
"""Extend the current Chunc startup spatial refinement to 800 cells."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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
TARGET_CELLS = (100, 200, 400, 800, 1600)
METRICS = (
    "peak_pressure_mpa",
    "peak_time_s",
    "max_fill_mach",
    "minimum_pressure_pa",
)
SUCCESSIVE_CHANGE_LIMIT_PERCENT = 2.0


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


def _successive_changes(rows: list[dict], metric: str) -> list[dict]:
    changes = []
    for coarse, fine in zip(rows, rows[1:]):
        coarse_value = float(coarse[metric])
        fine_value = float(fine[metric])
        relative = (
            abs(fine_value - coarse_value) / abs(fine_value) * 100.0
            if fine_value != 0.0 else None
        )
        changes.append({
            "from_target_cells": coarse["target_propellant_cells"],
            "to_target_cells": fine["target_propellant_cells"],
            "relative_change_percent": relative,
        })
    return changes


def _approximate_observed_order(values: list[float]) -> dict:
    """Return factor-two observed orders only for monotone triplets."""
    triplets = []
    for index in range(len(values) - 2):
        coarse, medium, fine = (float(value) for value in values[index:index + 3])
        first_difference = medium - coarse
        second_difference = fine - medium
        order = None
        extrapolated = None
        if first_difference * second_difference > 0.0 and second_difference != 0.0:
            ratio = abs(first_difference / second_difference)
            if ratio > 0.0:
                candidate_order = math.log(ratio, 2.0)
                denominator = 2.0 ** candidate_order - 1.0
                if candidate_order > 0.0 and denominator != 0.0:
                    order = candidate_order
                    extrapolated = fine + second_difference / denominator
        triplets.append({
            "indices": [index, index + 1, index + 2],
            "monotone": first_difference * second_difference > 0.0,
            "approximate_order": order,
            "approximate_extrapolated_value": extrapolated,
            "qualification": (
                "Uses target-cell factor 2; actual cell counts include fixed "
                "gap cells and are recorded separately."
            ),
        })
    return {"triplets": triplets}


def run_extension(
    output_root: Path | None = None,
    resume_path: Path | None = None,
) -> Path:
    baseline = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))["cases"][
        "chunc"
    ]
    motor_path = PROJECT_ROOT / baseline["motor_path"]
    options = dict(baseline["run_options"])
    options.update({
        "t_max": 0.05,
        "cfl_target": 0.3,
        "snapshot_interval": 0.01,
        "print_interval": 0.0,
    })
    if resume_path is not None:
        result_path = resume_path.resolve()
        output = result_path.parent
        payload = json.loads(result_path.read_text(encoding="utf-8"))
        provenance = payload.get("provenance", {})
        if provenance.get("motor_sha256") != _sha256(motor_path):
            raise ValueError("resume motor checksum does not match current input")
        if provenance.get("baseline_config_sha256") != _sha256(CONFIG_PATH):
            raise ValueError("resume baseline-config checksum does not match")
        payload["status"] = "running"
        payload["resumed_utc"] = datetime.now(timezone.utc).isoformat()
        payload["scope"]["target_propellant_cells"] = list(TARGET_CELLS)
    else:
        output = artifact_dir(
            "verification_chunc_spatial_extension",
            root=output_root if output_root is not None else PROJECT_ROOT,
        )
        result_path = output / "spatial_extension.json"
        payload = {
            "schema_version": 1,
            "status": "running",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "study": "chunc_current_mtv_startup_spatial_extension",
            "scope": {
                "target_propellant_cells": list(TARGET_CELLS),
                "metrics": list(METRICS),
                "successive_change_screen_percent": (
                    SUCCESSIVE_CHANGE_LIMIT_PERCENT
                ),
                "qualification": (
                    "Startup numerical screening only; not physics validation "
                    "or a full-burn convergence claim."
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
            },
            "configuration": options,
            "claims": {
                "spatially_converged_reference_established": False,
                "research_note_results_reproduced": False,
                "physics_validated": False,
            },
            "grid": [],
        }
    result_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    completed_cells = {
        int(row["target_propellant_cells"]) for row in payload["grid"]
    }
    for cells in TARGET_CELLS:
        if cells in completed_cells:
            continue
        result, performance, *_ = run_from_ric(
            str(motor_path), verbose=False,
            **dict(options, target_propellant_cells=cells),
        )
        verify_run_health(
            result, motor_name=f"chunc-spatial-extension cells={cells}",
            min_t_burn_s=0.049, raise_on_fail=True,
        )
        payload["grid"].append({
            "target_propellant_cells": cells,
            **_extract_metrics(result, performance),
        })
        result_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )

    payload["analysis"] = {}
    for metric in METRICS:
        changes = _successive_changes(payload["grid"], metric)
        order = _approximate_observed_order(
            [row[metric] for row in payload["grid"]]
        )
        payload["analysis"][metric] = {
            "screen_applicable": True,
            "successive_changes": changes,
            "finest_change_below_screen": (
                changes[-1]["relative_change_percent"] is not None
                and changes[-1]["relative_change_percent"]
                <= SUCCESSIVE_CHANGE_LIMIT_PERCENT
            ),
            "all_successive_changes_below_screen": all(
                item["relative_change_percent"] is not None
                and item["relative_change_percent"]
                <= SUCCESSIVE_CHANGE_LIMIT_PERCENT
                for item in changes
            ),
            **order,
        }
    floor_rows = [
        row for row in payload["grid"]
        if row.get("limit_pressure_floor_total_activations", 0) > 0
    ]
    payload["analysis"]["minimum_pressure_floor_censoring"] = {
        "censored": bool(floor_rows),
        "first_target_cells": (
            floor_rows[0]["target_propellant_cells"] if floor_rows else None
        ),
        "qualification": (
            "Minimum pressure at and beyond this grid is a limiter threshold, "
            "not an unconstrained solution value."
        ),
    }
    if floor_rows:
        pressure_analysis = payload["analysis"]["minimum_pressure_pa"]
        pressure_analysis["screen_applicable"] = False
        pressure_analysis["finest_change_below_screen"] = None
        pressure_analysis["all_successive_changes_below_screen"] = None
    payload["status"] = "complete"
    payload["completed_utc"] = datetime.now(timezone.utc).isoformat()
    result_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    output = run_extension(args.output_root, args.resume)
    print(f"Spatial-extension artifacts: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
