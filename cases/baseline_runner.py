"""Run a registered baseline and save a complete provenance manifest."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import math
import os
import platform
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib

matplotlib.use("Agg")

from srm_1d.openmotor_adapter import (
    load_ric,
    ric_to_sim_args,
    run_from_ric,
    save_csv,
)
from srm_1d.plotting import plot_flow_snapshot, plot_pressure, plot_summary
from srm_1d.run_artifacts import artifact_dir, verify_run_health
from srm_1d.simulation import run_simulation

from .loader import PROJECT_ROOT, get_case, load_measurement


DEFAULT_CONFIGS = Path(__file__).with_name("baseline_configs.json")
OPENMOTOR_PATH_ENV = "SRM1D_OPENMOTOR_PATH"


def _json_safe(value: Any) -> Any:
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return repr(value)


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(
        json.dumps(_json_safe(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_state(repo: Path) -> dict:
    def git(*args: str) -> str | None:
        completed = subprocess.run(
            ["git", *args],
            cwd=repo,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return completed.stdout.strip() if completed.returncode == 0 else None

    status = git("status", "--short")
    return {
        "path": str(repo.resolve()),
        "commit": git("rev-parse", "HEAD"),
        "branch": git("branch", "--show-current"),
        "origin": git("remote", "get-url", "origin"),
        "dirty": bool(status),
        "status": status or "",
    }


def _call_defaults(function) -> dict:
    defaults = {}
    for name, parameter in inspect.signature(function).parameters.items():
        if parameter.default is not inspect.Parameter.empty:
            defaults[name] = _json_safe(parameter.default)
    return defaults


def _dependency_versions() -> dict[str, str]:
    packages = [
        "numpy",
        "scipy",
        "numba",
        "llvmlite",
        "PyYAML",
        "matplotlib",
        "scikit-fmm",
        "scikit-image",
        "pytest",
        "PyQt6",
    ]
    versions = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def _registered_files(case: dict) -> list[dict]:
    files = []
    for entry in case["files"]:
        item = dict(entry)
        path = PROJECT_ROOT / entry["path"]
        item["exists"] = path.is_file()
        if path.is_file():
            item["observed_sha256"] = _sha256(path)
        files.append(item)
    return files


def _comparison_metrics(result: dict, comparison: dict) -> dict:
    if not comparison.get("enabled", False):
        return {
            "status": "not_computed",
            "reason": comparison.get("reason", "comparison disabled"),
        }

    trace = load_measurement(comparison["measurement_case_id"])
    offset = float(comparison["time_offset_s"])
    sample_time = np.asarray(trace["time"], dtype=float) + offset
    sim_time = np.asarray(result["time"], dtype=float)
    sim_pressure_mpa = np.asarray(result["P_head"], dtype=float) / 1e6
    measured_pressure_mpa = np.asarray(trace["pressure"], dtype=float)
    mask = (sample_time >= sim_time[0]) & (sample_time <= sim_time[-1])
    if not np.any(mask):
        return {
            "status": "not_computed",
            "reason": "no aligned measurement samples overlap the simulation",
        }

    predicted = np.interp(sample_time[mask], sim_time, sim_pressure_mpa)
    residual = predicted - measured_pressure_mpa[mask]
    return {
        "status": "computed",
        "alignment_method": comparison["alignment_method"],
        "time_offset_s": offset,
        "samples_used": int(np.count_nonzero(mask)),
        "rmse_mpa": float(np.sqrt(np.mean(residual**2))),
        "mae_mpa": float(np.mean(np.abs(residual))),
        "bias_mpa": float(np.mean(residual)),
        "measured_peak_mpa": float(np.max(measured_pressure_mpa[mask])),
        "simulated_at_measurement_peak_mpa": float(
            predicted[np.argmax(measured_pressure_mpa[mask])]
        ),
        "qualification": "newly computed with a manually aligned digitized trace",
    }


def _metrics(result: dict, performance: dict, comparison: dict) -> dict:
    summary = result["summary"]
    return {
        "simulation": {
            "termination": summary.get("termination"),
            "termination_code": summary.get("termination_code"),
            "t_burn_s": summary.get("t_burn"),
            "t_first_burnout_s": summary.get("t_first_burnout"),
            "peak_pressure_pa": summary.get("P_peak"),
            "peak_pressure_time_s": summary.get("t_peak"),
            "mid_burn_pressure_pa": summary.get("P_mid"),
            "mass_balance_error_fraction": summary.get("mass_balance_error"),
            "propellant_mass_kg": summary.get("propellant_mass"),
            "mass_produced_kg": summary.get("mass_produced"),
            "mass_nozzle_kg": summary.get("mass_nozzle"),
            "steps": summary.get("steps"),
            "cells": summary.get("cells"),
            "dt_min_s": summary.get("dt_min"),
            "dt_median_s": summary.get("dt_median"),
            "dt_final_s": summary.get("dt_final"),
            "wall_time_s": summary.get("wall_time"),
        },
        "performance": {
            "designation": performance.get("designation"),
            "total_impulse_ns": performance.get("total_impulse"),
            "average_thrust_n": performance.get("average_thrust"),
            "peak_thrust_n": performance.get("peak_thrust"),
            "average_isp_s": performance.get("average_Isp"),
            "burn_time_s": performance.get("burn_time"),
            "c_star_m_s": performance.get("c_star"),
        },
        "measurement_comparison": _comparison_metrics(result, comparison),
    }


def _new_manifest(case_id: str, case: dict, config: dict, output: Path) -> dict:
    openmotor_value = os.environ.get(OPENMOTOR_PATH_ENV)
    openmotor_path = Path(openmotor_value).resolve() if openmotor_value else None
    run_defaults = _call_defaults(run_from_ric)
    simulation_defaults = _call_defaults(run_simulation)
    adapter_defaults = _call_defaults(ric_to_sim_args)
    motor = load_ric(str(PROJECT_ROOT / config["motor_path"]))
    resolved_simulation = dict(simulation_defaults)
    for name, value in config["run_options"].items():
        if name in resolved_simulation:
            resolved_simulation[name] = value
    resolved_simulation["P_ambient"] = motor.get("config", {}).get(
        "ambPressure", 101325.0
    )
    return {
        "schema_version": 1,
        "status": "prepared",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "case_id": case_id,
        "case_readiness_at_run": case["readiness"],
        "output_directory": str(output.resolve()),
        "command": [sys.executable, *sys.argv],
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "python_executable": sys.executable,
        },
        "repositories": {
            "srm_1d": _git_state(PROJECT_ROOT),
            "openMotor": (
                _git_state(openmotor_path)
                if openmotor_path and openmotor_path.is_dir()
                else {"path": openmotor_value, "available": False}
            ),
        },
        "dependencies": _dependency_versions(),
        "registered_files": _registered_files(case),
        "active_inputs": {
            "motor": config["motor_path"],
            "configuration_source": config["configuration_source"],
            "transport": config["transport"],
            "measurement_case_id": config["comparison"]["measurement_case_id"],
            "measurement_active": config["comparison"]["enabled"],
        },
        "configuration": {
            "requested": config,
            "run_from_ric_defaults": run_defaults,
            "ric_to_sim_args_defaults": adapter_defaults,
            "run_simulation_defaults": simulation_defaults,
            "resolved_run_options": {
                **run_defaults,
                **config["run_options"],
            },
            "resolved_simulation_options": resolved_simulation,
            "target_propellant_cells": adapter_defaults[
                "target_propellant_cells"
            ],
        },
        "claims": {
            "historical_research_results_reproduced": False,
            "note": "Only metrics in this run directory are established by this execution."
        }
    }


def run_baseline(
    case_id: str,
    *,
    configs_path: Path = DEFAULT_CONFIGS,
    output_root: Path | None = None,
    save_plots: bool = False,
    dry_run: bool = False,
) -> Path:
    configs = json.loads(configs_path.read_text(encoding="utf-8"))
    try:
        config = configs["cases"][case_id]
    except KeyError as exc:
        raise KeyError(f"No baseline configuration registered for {case_id!r}") from exc
    case = get_case(case_id)

    output = artifact_dir(
        f"baseline_{case_id}",
        root=output_root if output_root is not None else PROJECT_ROOT,
    )
    manifest_path = output / "run_manifest.json"
    manifest = _new_manifest(case_id, case, config, output)
    _write_json(manifest_path, manifest)
    if dry_run:
        manifest["status"] = "dry_run_complete"
        _write_json(manifest_path, manifest)
        return output

    options = dict(config["run_options"])
    options["verbose"] = True
    transport_path = config["transport"].get("transport_path")
    if transport_path is not None:
        options["transport_path"] = str(PROJECT_ROOT / transport_path)

    try:
        result, performance, nozzle, geometry, propellant = run_from_ric(
            str(PROJECT_ROOT / config["motor_path"]),
            **options,
        )
        healthy = verify_run_health(
            result,
            motor_name=f"baseline:{case_id}",
            raise_on_fail=True,
        )
        metrics = _metrics(result, performance, config["comparison"])
        metrics["run_health_passed"] = healthy
        _write_json(output / "metrics.json", metrics)
        manifest["resolved_runtime"] = {
            "geometry": {
                "motor_length_m": geometry.L_motor,
                "outer_diameter_m": geometry.D_outer,
                "segment_count": len(geometry.segments),
                "cell_count": geometry.N_cells,
                "cell_width_m": geometry.dx,
            },
            "nozzle": {
                "throat_diameter_m": nozzle.D_throat,
                "exit_diameter_m": nozzle.D_exit,
                "expansion_ratio": nozzle.expansion_ratio,
                "efficiency": nozzle.efficiency,
                "divergence_angle_deg": nozzle.div_angle,
                "convergence_angle_deg": nozzle.conv_angle,
                "throat_length_m": nozzle.throat_length,
                "erosion_coefficient_um_s_mpa": nozzle.erosion_coeff,
                "slag_coefficient_m_mpa_s": nozzle.slag_coeff,
            },
            "propellant": {
                "density_kg_m3": propellant.rho_propellant,
                "solid_conductivity_w_m_k": propellant.k_solid,
                "radiation_emissivity": propellant.radiation_emissivity,
                "flame_front_enabled": propellant.flame_front_enabled,
                "zn_enabled": propellant.zn_enabled,
            },
        }
        save_csv(
            str(output / "result.csv"),
            result,
            performance,
            geometry,
            propellant,
            dt_sample=0.001,
        )

        if save_plots:
            experimental = None
            offset = 0.0
            if config["comparison"].get("enabled", False):
                experimental = load_measurement(
                    config["comparison"]["measurement_case_id"]
                )
                offset = float(config["comparison"]["time_offset_s"])
            plot_pressure(
                result,
                experimental=experimental,
                time_offset=offset,
                title=f"Baseline: {case['label']}",
                save_path=str(output / "pressure.png"),
            )
            plot_flow_snapshot(
                result,
                t_target=config["plots"]["flow_snapshot_time_s"],
                save_path=str(output / "flow.png"),
            )
            plot_summary(
                result,
                performance=performance,
                experimental=experimental,
                time_offset=offset,
                title=f"Baseline: {case['label']}",
                save_path=str(output / "summary.png"),
            )

        manifest["status"] = "complete"
        manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["outputs"] = {
            path.name: _sha256(path)
            for path in sorted(output.iterdir())
            if path.is_file() and path.name != manifest_path.name
        }
        _write_json(manifest_path, manifest)
        return output
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        _write_json(manifest_path, manifest)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, choices=["hasegawa_a", "chunc"])
    parser.add_argument("--configs", type=Path, default=DEFAULT_CONFIGS)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--plots", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    output = run_baseline(
        args.case,
        configs_path=args.configs,
        output_root=args.output_root,
        save_plots=args.plots,
        dry_run=args.dry_run,
    )
    print(f"Baseline artifacts: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
