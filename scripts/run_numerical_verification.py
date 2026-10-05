#!/usr/bin/env python3
"""Run portable Chunc startup or full-burn grid/CFL verification."""

from __future__ import annotations

import argparse
import csv
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
from srm_1d.openmotor_adapter import run_from_ric
from srm_1d.run_artifacts import artifact_dir, verify_run_health
from srm_1d.tools.ignition_diagnostics import sensible_enthalpy_audit


CONFIG_PATH = PROJECT_ROOT / "cases" / "baseline_configs.json"
GRID_LEVELS = (50, 100, 200)
STARTUP_CFL_LEVELS = (0.3, 0.15, 0.075)
FULL_CFL_LEVELS = (0.3, 0.2, 0.1)
REFERENCE_CELLS = 100
REFERENCE_CFL = 0.3
STARTUP_DURATION_S = 0.05
FILL_WINDOW_S = 0.003
PRESSURE_SCREENING_TOLERANCE_PERCENT = 2.0


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


def _write_limit_csv(path: Path, result: dict) -> None:
    limits = result["numerical_limits"]
    names = ("temperature_floor", "temperature_ceiling", "pressure_floor")
    fieldnames = ["cell_index", "x_m"]
    for name in names:
        fieldnames.extend((
            f"{name}_activations",
            f"{name}_duration_s",
            f"{name}_first_activation_time_s",
            f"{name}_last_activation_time_s",
            f"{name}_absolute_correction_energy_j",
        ))
    fieldnames.append("pressure_floor_maximum_raw_deficit_pa")
    pressure_floor_events = limits["pressure_floor"].get(
        "maximum_deficit_event", {}
    )
    fieldnames.extend(
        f"pressure_floor_maximum_deficit_event_{field}"
        for field in pressure_floor_events
    )
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for index, x_m in enumerate(limits["x_m"]):
            row = {"cell_index": index, "x_m": float(x_m)}
            for name in names:
                diagnostic = limits[name]
                row[f"{name}_activations"] = int(
                    diagnostic["activation_count_by_cell"][index]
                )
                row[f"{name}_duration_s"] = float(
                    diagnostic["duration_s_by_cell"][index]
                )
                row[f"{name}_first_activation_time_s"] = float(
                    diagnostic["first_activation_time_s_by_cell"][index]
                )
                row[f"{name}_last_activation_time_s"] = float(
                    diagnostic["last_activation_time_s_by_cell"][index]
                )
                energy = diagnostic.get("absolute_correction_energy_j_by_cell")
                row[f"{name}_absolute_correction_energy_j"] = (
                    float(energy[index]) if energy is not None else 0.0
                )
            deficit = limits["pressure_floor"].get(
                "maximum_raw_deficit_pa_by_cell"
            )
            row["pressure_floor_maximum_raw_deficit_pa"] = (
                float(deficit[index]) if deficit is not None else 0.0
            )
            for field, values in pressure_floor_events.items():
                row[f"pressure_floor_maximum_deficit_event_{field}"] = (
                    np.asarray(values)[index].item()
                )
            writer.writerow(row)


def _write_mach_limit_csv(path: Path, result: dict) -> None:
    diagnostic = result["numerical_limits"]["port_mach_cap"]
    fieldnames = [
        "interior_face_index", "x_m", "activations", "duration_s",
        "first_activation_time_s", "last_activation_time_s",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for index, x_m in enumerate(diagnostic["x_m"]):
            writer.writerow({
                "interior_face_index": index + 1,
                "x_m": float(x_m),
                "activations": int(diagnostic["activation_count_by_face"][index]),
                "duration_s": float(diagnostic["duration_s_by_face"][index]),
                "first_activation_time_s": float(
                    diagnostic["first_activation_time_s_by_face"][index]
                ),
                "last_activation_time_s": float(
                    diagnostic["last_activation_time_s_by_face"][index]
                ),
            })


def _ignition_refresh_delays(result: dict, interval: int) -> np.ndarray:
    """Return each observed ignition's delay to the next rate refresh."""
    time = np.asarray(result["time"], dtype=float)
    ignition_times = np.asarray(result["ignition_time_by_cell"], dtype=float)
    finite = ignition_times[
        np.isfinite(ignition_times) & (ignition_times < 1.0e9)
    ]
    delays = []
    for ignition_time in finite:
        ignition_step = int(np.searchsorted(time, ignition_time, side="left"))
        next_refresh_step = (ignition_step // interval + 1) * interval
        if next_refresh_step < len(time):
            delays.append(max(0.0, float(time[next_refresh_step] - ignition_time)))
    return np.asarray(delays, dtype=float)


def _extract_metrics(result: dict, performance: dict | None = None) -> dict:
    time = np.asarray(result["time"], dtype=float)
    pressure = np.asarray(result["P_head"], dtype=float)
    mach = np.asarray(result["max_mach"], dtype=float)
    fill = time <= FILL_WINDOW_S
    peak_index = int(np.argmax(pressure))
    clipping = np.asarray(result["clipping_correction_power"], dtype=float)
    residual = np.asarray(result["energy_residual"], dtype=float)
    thermal = np.asarray(result["thermal_source_power"], dtype=float)
    residual_integral = float(np.trapz(np.abs(residual), time))
    thermal_integral = float(np.trapz(np.abs(thermal), time))
    clipping_integral = float(np.trapz(np.abs(clipping), time))
    clipping_steps = int(np.count_nonzero(clipping))
    summary = result["summary"]
    energy_audit = sensible_enthalpy_audit(result)
    source_mismatch = energy_audit["source_ledger_mismatch_power"]
    nozzle_mismatch = energy_audit["nozzle_boundary_mismatch_power"]
    storage_jump = energy_audit["between_step_storage_jump"]
    metrics = {
        "actual_cells": int(summary["cells"]),
        "steps": int(summary["steps"]),
        "termination": summary["termination"],
        "peak_pressure_mpa": float(pressure[peak_index] / 1e6),
        "peak_time_s": float(time[peak_index]),
        "max_fill_mach": float(np.max(mach[fill])) if np.any(fill) else None,
        "minimum_pressure_pa": float(np.min(result["min_pressure"])),
        "maximum_temperature_k": float(np.max(result["max_gas_temperature"])),
        "clipping_active_steps": clipping_steps,
        "clipping_active_fraction": clipping_steps / len(clipping),
        "max_abs_clipping_power_w": float(np.max(np.abs(clipping))),
        "integrated_abs_clipping_energy_j": clipping_integral,
        "integrated_abs_clipping_over_thermal_source": (
            clipping_integral / thermal_integral if thermal_integral > 0.0 else None
        ),
        "max_abs_energy_residual_w": float(np.max(np.abs(residual))),
        "integrated_abs_energy_residual_j": residual_integral,
        "integrated_abs_residual_over_thermal_source": (
            residual_integral / thermal_integral if thermal_integral > 0.0 else None
        ),
        "energy_balance_classification": summary.get("energy_convention", {}).get(
            "classification", "sensible-enthalpy scalar transport balance"
        ),
        "max_abs_source_ledger_mismatch_w": float(np.max(np.abs(source_mismatch))),
        "integrated_abs_source_ledger_mismatch_j": float(
            np.trapz(np.abs(source_mismatch), time)
        ),
        "max_abs_nozzle_boundary_mismatch_w": float(np.max(np.abs(nozzle_mismatch))),
        "integrated_abs_nozzle_boundary_mismatch_j": float(
            np.trapz(np.abs(nozzle_mismatch), time)
        ),
        "max_abs_between_step_storage_jump_j": float(np.max(np.abs(storage_jump))),
        "sum_abs_between_step_storage_jump_j": float(np.sum(np.abs(storage_jump))),
        "mass_balance_error_fraction": float(summary["mass_balance_error"]),
        "dt_min_s": float(summary["dt_min"]),
        "dt_median_s": float(summary["dt_median"]),
        "wall_time_s": float(summary["wall_time"]),
        "first_burnout_time_s": summary.get("t_first_burnout"),
        "final_throat_diameter_m": summary.get("D_throat_final"),
    }
    burn_update_interval = summary.get("burn_update_interval")
    if burn_update_interval is not None:
        burn_update_interval = int(burn_update_interval)
        delays = _ignition_refresh_delays(result, burn_update_interval)
        metrics["burn_update_interval_steps"] = burn_update_interval
        metrics["ignition_rate_refresh_delay_count"] = int(len(delays))
        metrics["ignition_rate_refresh_delay_max_s"] = (
            float(np.max(delays)) if len(delays) else None
        )
        metrics["ignition_rate_refresh_delay_median_s"] = (
            float(np.median(delays)) if len(delays) else None
        )
    geometry_update_interval = summary.get("geometry_update_interval")
    if geometry_update_interval is not None:
        metrics["geometry_update_interval_steps"] = int(
            geometry_update_interval
        )
    if performance is not None:
        metrics.update({
            "total_impulse_ns": performance.get("total_impulse"),
            "average_thrust_n": performance.get("average_thrust"),
            "peak_thrust_n": performance.get("peak_thrust"),
            "performance_burn_time_s": performance.get("burn_time"),
            "average_isp_s": performance.get("average_Isp"),
        })
    limits = result.get("numerical_limits")
    if limits is not None:
        x_m = np.asarray(limits["x_m"], dtype=float)
        for name in ("temperature_floor", "temperature_ceiling", "pressure_floor"):
            diagnostic = limits[name]
            counts = np.asarray(diagnostic["activation_count_by_cell"], dtype=np.int64)
            durations = np.asarray(diagnostic["duration_s_by_cell"], dtype=float)
            max_index = int(np.argmax(durations))
            prefix = f"limit_{name}"
            metrics[f"{prefix}_active_cells"] = int(np.count_nonzero(counts))
            metrics[f"{prefix}_total_activations"] = int(np.sum(counts))
            metrics[f"{prefix}_total_cell_duration_s"] = float(np.sum(durations))
            metrics[f"{prefix}_max_cell_duration_s"] = float(durations[max_index])
            metrics[f"{prefix}_max_duration_x_m"] = float(x_m[max_index])
            first_times = np.asarray(
                diagnostic["first_activation_time_s_by_cell"], dtype=float
            )
            last_times = np.asarray(
                diagnostic["last_activation_time_s_by_cell"], dtype=float
            )
            metrics[f"{prefix}_first_activation_time_s"] = (
                float(np.nanmin(first_times)) if np.any(np.isfinite(first_times)) else None
            )
            metrics[f"{prefix}_last_activation_time_s"] = (
                float(np.nanmax(last_times)) if np.any(np.isfinite(last_times)) else None
            )
            if name == "pressure_floor":
                metrics[f"{prefix}_threshold_pa"] = float(
                    diagnostic["threshold_pa"]
                )
            energy = diagnostic.get("absolute_correction_energy_j_by_cell")
            if energy is not None:
                energy = np.asarray(energy, dtype=float)
                energy_index = int(np.argmax(energy))
                metrics[f"{prefix}_total_abs_energy_j"] = float(np.sum(energy))
                metrics[f"{prefix}_max_energy_x_m"] = float(x_m[energy_index])
            deficit = diagnostic.get("maximum_raw_deficit_pa_by_cell")
            if deficit is not None:
                deficit = np.asarray(deficit, dtype=float)
                deficit_index = int(np.argmax(deficit))
                metrics[f"{prefix}_maximum_raw_deficit_pa"] = float(
                    deficit[deficit_index]
                )
                metrics[f"{prefix}_maximum_raw_deficit_x_m"] = (
                    float(x_m[deficit_index])
                    if deficit[deficit_index] > 0.0 else None
                )
                event = diagnostic.get("maximum_deficit_event", {})
                for field, values in event.items():
                    value = np.asarray(values)[deficit_index].item()
                    metrics[
                        f"{prefix}_maximum_deficit_event_{field}"
                    ] = value if deficit[deficit_index] > 0.0 else None
        diagnostic = limits["port_mach_cap"]
        counts = np.asarray(diagnostic["activation_count_by_face"], dtype=np.int64)
        durations = np.asarray(diagnostic["duration_s_by_face"], dtype=float)
        first_times = np.asarray(
            diagnostic["first_activation_time_s_by_face"], dtype=float
        )
        last_times = np.asarray(
            diagnostic["last_activation_time_s_by_face"], dtype=float
        )
        metrics["limit_port_mach_cap_threshold_mach"] = float(
            diagnostic["threshold_mach"]
        )
        metrics["limit_port_mach_cap_active_faces"] = int(np.count_nonzero(counts))
        metrics["limit_port_mach_cap_total_activations"] = int(np.sum(counts))
        metrics["limit_port_mach_cap_total_face_duration_s"] = float(
            np.sum(durations)
        )
        metrics["limit_port_mach_cap_first_activation_time_s"] = (
            float(np.nanmin(first_times)) if np.any(np.isfinite(first_times)) else None
        )
        metrics["limit_port_mach_cap_last_activation_time_s"] = (
            float(np.nanmax(last_times)) if np.any(np.isfinite(last_times)) else None
        )
    return metrics


def _adjacent_changes(rows: list[dict], key: str) -> list[dict]:
    changes = []
    for coarse, fine in zip(rows, rows[1:]):
        denominator = abs(float(fine[key]))
        relative = (
            abs(float(fine[key]) - float(coarse[key])) / denominator * 100.0
            if denominator > 0.0 else None
        )
        changes.append({
            "from": coarse["label"],
            "to": fine["label"],
            "metric": key,
            "relative_change_percent": relative,
        })
    return changes


def run_single_point(
    profile: str,
    cells: int,
    cfl: float,
    history_capacity: int | None = None,
    output_root: Path | None = None,
    port_mach_cap: float | None = None,
    burn_update_interval: int | None = None,
    geometry_update_interval: int | None = None,
) -> Path:
    if profile not in {"startup", "full"}:
        raise ValueError(f"Unknown verification profile: {profile!r}")
    baseline = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))["cases"]["chunc"]
    motor_path = PROJECT_ROOT / baseline["motor_path"]
    options = dict(baseline["run_options"])
    if profile == "startup":
        options.update({
            "t_max": STARTUP_DURATION_S,
            "snapshot_interval": 0.01,
            "print_interval": 0.0,
        })
        minimum_run_time = 0.049
    else:
        options["print_interval"] = 0.0
        minimum_run_time = float(options["t_max"]) * 0.95
    options.update({
        "target_propellant_cells": cells,
        "cfl_target": cfl,
    })
    if history_capacity is not None:
        options["history_capacity"] = history_capacity
    if port_mach_cap is not None:
        options["port_mach_cap"] = port_mach_cap
    if burn_update_interval is not None:
        options["burn_update_interval"] = burn_update_interval
    if geometry_update_interval is not None:
        options["geometry_update_interval"] = geometry_update_interval

    output = artifact_dir(
        f"verification_chunc_{profile}_point",
        root=output_root if output_root is not None else PROJECT_ROOT,
    )
    point_path = output / "point.json"
    payload = {
        "schema_version": 1,
        "status": "running",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "study": f"chunc_{profile}_single_point",
        "configuration": options,
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
            "full_burn_numerically_converged": False,
            "research_note_results_reproduced": False,
        },
    }
    _write_json(point_path, payload)
    try:
        result, performance, *_ = run_from_ric(
            str(motor_path), verbose=False, **options,
        )
        verify_run_health(
            result,
            motor_name=f"chunc-{profile} cells={cells} cfl={cfl}",
            min_t_burn_s=minimum_run_time,
            raise_on_fail=True,
        )
        payload["status"] = "complete"
        payload["completed_utc"] = datetime.now(timezone.utc).isoformat()
        payload["metrics"] = _extract_metrics(result, performance)
        limit_path = output / "numerical_limits_by_cell.csv"
        _write_limit_csv(limit_path, result)
        mach_limit_path = output / "mach_limit_by_face.csv"
        _write_mach_limit_csv(mach_limit_path, result)
        payload["outputs"] = {
            limit_path.name: _sha256(limit_path),
            mach_limit_path.name: _sha256(mach_limit_path),
        }
        _write_json(point_path, payload)
        return output
    except Exception as exc:
        payload["status"] = "failed"
        payload["completed_utc"] = datetime.now(timezone.utc).isoformat()
        payload["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
        }
        _write_json(point_path, payload)
        raise


def run_study(profile: str = "startup", output_root: Path | None = None) -> Path:
    if profile not in {"startup", "full"}:
        raise ValueError(f"Unknown verification profile: {profile!r}")
    baseline = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))["cases"]["chunc"]
    motor_path = PROJECT_ROOT / baseline["motor_path"]
    options = dict(baseline["run_options"])
    if profile == "startup":
        cfl_levels = STARTUP_CFL_LEVELS
        options.update({
            "t_max": STARTUP_DURATION_S,
            "snapshot_interval": 0.01,
            "print_interval": 0.0,
        })
        minimum_run_time = 0.049
        qualification = (
            "Startup screening only; this is not a full-burn convergence claim."
        )
    else:
        cfl_levels = FULL_CFL_LEVELS
        options["print_interval"] = 0.0
        minimum_run_time = float(options["t_max"]) * 0.95
        qualification = (
            "Full configured duration; convergence applies only to recorded metrics."
        )

    output = artifact_dir(
        f"verification_chunc_{profile}",
        root=output_root if output_root is not None else PROJECT_ROOT,
    )
    progress_path = output / "progress.json"
    progress = {
        "schema_version": 1,
        "status": "running",
        "profile": profile,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "completed_points": [],
        "failed_points": [],
    }
    _write_json(progress_path, progress)

    points = {(cells, REFERENCE_CFL) for cells in GRID_LEVELS}
    points.update((REFERENCE_CELLS, cfl) for cfl in cfl_levels)
    results = {}
    for cells, cfl in sorted(points, key=lambda item: (item[0], -item[1])):
        point_options = dict(options, target_propellant_cells=cells, cfl_target=cfl)
        try:
            result, performance, *_ = run_from_ric(
                str(motor_path), verbose=False, **point_options,
            )
            verify_run_health(
                result, motor_name=f"chunc-{profile} cells={cells} cfl={cfl}",
                min_t_burn_s=minimum_run_time, raise_on_fail=True,
            )
            point_result = {
                "target_propellant_cells": cells,
                "cfl_target": cfl,
                **_extract_metrics(result, performance),
            }
            results[(cells, cfl)] = point_result
            progress["completed_points"].append(point_result)
            _write_json(progress_path, progress)
        except Exception as exc:
            progress["status"] = "failed"
            progress["failed_points"].append({
                "target_propellant_cells": cells,
                "cfl_target": cfl,
                "error_type": type(exc).__name__,
                "error": str(exc),
            })
            _write_json(progress_path, progress)
            raise

    spatial = []
    for cells in GRID_LEVELS:
        row = dict(results[(cells, REFERENCE_CFL)])
        row["label"] = f"cells={cells}"
        spatial.append(row)
    temporal = []
    for cfl in cfl_levels:
        row = dict(results[(REFERENCE_CELLS, cfl)])
        row["label"] = f"cfl={cfl:g}"
        temporal.append(row)

    comparisons = {
        "spatial_peak_pressure": _adjacent_changes(spatial, "peak_pressure_mpa"),
        "spatial_peak_time": _adjacent_changes(spatial, "peak_time_s"),
        "spatial_fill_mach": _adjacent_changes(spatial, "max_fill_mach"),
        "spatial_minimum_pressure": _adjacent_changes(spatial, "minimum_pressure_pa"),
        "temporal_peak_pressure": _adjacent_changes(temporal, "peak_pressure_mpa"),
        "temporal_peak_time": _adjacent_changes(temporal, "peak_time_s"),
        "temporal_fill_mach": _adjacent_changes(temporal, "max_fill_mach"),
        "temporal_minimum_pressure": _adjacent_changes(temporal, "minimum_pressure_pa"),
    }
    if profile == "full":
        for prefix, rows in (("spatial", spatial), ("temporal", temporal)):
            for metric in (
                "total_impulse_ns",
                "performance_burn_time_s",
                "first_burnout_time_s",
                "final_throat_diameter_m",
            ):
                comparisons[f"{prefix}_{metric}"] = _adjacent_changes(rows, metric)
    pressure_changes = (
        comparisons["spatial_peak_pressure"] + comparisons["temporal_peak_pressure"]
    )
    pressure_screen_passed = all(
        item["relative_change_percent"] <= PRESSURE_SCREENING_TOLERANCE_PERCENT
        for item in pressure_changes
    )

    payload = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "study": f"chunc_{profile}_spatial_and_cfl_refinement",
        "scope": {
            "profile": profile,
            "duration_s": options["t_max"],
            "fill_window_s": FILL_WINDOW_S,
            "spatial_levels": list(GRID_LEVELS),
            "cfl_levels": list(cfl_levels),
            "pressure_screening_tolerance_percent": PRESSURE_SCREENING_TOLERANCE_PERCENT,
            "qualification": qualification,
        },
        "provenance": {
            "python": sys.version,
            "platform": platform.platform(),
            "git_commit": _git_value("rev-parse", "HEAD"),
            "git_status": _git_value("status", "--short"),
            "motor_path": str(motor_path.relative_to(PROJECT_ROOT)),
            "motor_sha256": _sha256(motor_path),
            "baseline_config_sha256": _sha256(CONFIG_PATH),
            "baseline_configuration_source": baseline["configuration_source"],
        },
        "fixed_run_options": options,
        "spatial": spatial,
        "temporal": temporal,
        "comparisons": comparisons,
        "screening": {
            "all_peak_pressure_successive_changes_below_tolerance": pressure_screen_passed,
            "finest_spatial_peak_pressure_change_below_tolerance": (
                comparisons["spatial_peak_pressure"][-1]["relative_change_percent"]
                <= PRESSURE_SCREENING_TOLERANCE_PERCENT
            ),
            "finest_temporal_peak_pressure_change_below_tolerance": (
                comparisons["temporal_peak_pressure"][-1]["relative_change_percent"]
                <= PRESSURE_SCREENING_TOLERANCE_PERCENT
            ),
            "tolerance_is_proposed_not_experimental_accuracy": True,
        },
        "claims": {
            "full_burn_numerically_converged": False,
            "research_note_results_reproduced": False,
        },
    }
    result_path = output / "verification.json"
    _write_json(result_path, payload)

    rows = [dict(study="spatial", **row) for row in spatial]
    rows += [dict(study="temporal", **row) for row in temporal]
    with (output / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    progress["status"] = "complete"
    progress["completed_utc"] = datetime.now(timezone.utc).isoformat()
    _write_json(progress_path, progress)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=["startup", "full"], default="startup")
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--single-cells", type=int)
    parser.add_argument("--single-cfl", type=float)
    parser.add_argument("--history-capacity", type=int)
    parser.add_argument("--single-port-mach-cap", type=float)
    parser.add_argument("--single-burn-update-interval", type=int)
    parser.add_argument("--single-geometry-update-interval", type=int)
    args = parser.parse_args()
    single_requested = args.single_cells is not None or args.single_cfl is not None
    if single_requested:
        if args.single_cells is None or args.single_cfl is None:
            parser.error("--single-cells and --single-cfl must be supplied together")
        output = run_single_point(
            args.profile,
            args.single_cells,
            args.single_cfl,
            args.history_capacity,
            args.output_root,
            args.single_port_mach_cap,
            args.single_burn_update_interval,
            args.single_geometry_update_interval,
        )
    else:
        if (args.history_capacity is not None
                or args.single_port_mach_cap is not None
                or args.single_burn_update_interval is not None
                or args.single_geometry_update_interval is not None):
            parser.error(
                "single-point overrides require --single-cells and --single-cfl"
            )
        output = run_study(args.profile, args.output_root)
    print(f"Verification artifacts: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
