#!/usr/bin/env python3
"""Zerox startup pressure-shape comparison (protocol v0.1 compliant).

Runs the baseline Zerox motor (zerox.ric) with cross-motor nominal settings
(v0.7.5) and compares the simulated head-end pressure against the registered
native-trigger measurement in the startup window [0.0, 0.36] s.

Shape metrics only — amplitude scoring is blocked pending independent
manufacturer pressure-calibration documentation (see validation_case_registry.json
scoring_eligibility.absolute_pressure and ZEROX_TRANSDUCER_CALIBRATION_AUDIT.md).

No fitted model parameters are used here; the LHS-calibrated zerox_LHS.ric is
intentionally NOT used so this remains a protocol-compliant shape probe.

Outputs:
    docs/v0_8_0/data/zerox_startup_shape_comparison.json
"""

from __future__ import annotations

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

from cases.loader import load_measurement
from srm_1d.openmotor_adapter import run_from_ric


OUTPUT_PATH = SCRIPT_ROOT / "docs" / "v0_8_0" / "data" / "zerox_startup_shape_comparison.json"
MOTOR_PATH = SCRIPT_ROOT / "motors" / "zerox.ric"
STARTUP_WINDOW_S = (0.0, 0.36)

# v0.7.5 cross-motor nominal settings — identical to the Chunc matched-grid runs
# (SRM_Erosive_Prediction_Plan §5; docs/v0_7_5/RESULT.md)
NOMINAL_SETTINGS = dict(
    roughness=32e-6,
    kappa=0.44,
    pyrogen="bpnv",
    pyrogen_mass=None,
    T_ignition=756.0,
    P_cutoff=0.01e6,
    snapshot_interval=0.2,
    print_interval=0.0,
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_value(*args: str) -> str:
    result = subprocess.run(
        ["git", *args], cwd=SCRIPT_ROOT, check=False,
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else ""


def _interp(xs: list[float], ys: list[float], x: float) -> float:
    """Linear interpolation; clamps to boundary values outside range."""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            dx = xs[i + 1] - xs[i]
            if dx == 0.0:
                return ys[i]
            t = (x - xs[i]) / dx
            return ys[i] + t * (ys[i + 1] - ys[i])
    return ys[-1]


def run_shape_comparison() -> dict:
    """Run the shape comparison and return the full result payload."""
    # --- Load registered measurement ---
    meas = load_measurement("zerox")
    meas_time = meas["time"]
    meas_pres = meas["pressure"]  # MPa, gauge relative to pre-event baseline

    # Filter to startup window
    win_lo, win_hi = STARTUP_WINDOW_S
    meas_win_t = [t for t in meas_time if win_lo <= t <= win_hi]
    meas_win_p = [p for t, p in zip(meas_time, meas_pres) if win_lo <= t <= win_hi]
    n_meas = len(meas_win_t)

    if n_meas < 2:
        raise RuntimeError(
            f"Only {n_meas} measurement samples in startup window; cannot compute shape metrics."
        )

    # --- Run simulator (baseline motor, nominal settings, t_max just past window) ---
    print(f"Running baseline Zerox simulation (t_max=0.5 s)...")
    result, perf, nozzle, geo, prop = run_from_ric(
        str(MOTOR_PATH),
        t_max=0.5,
        **NOMINAL_SETTINGS,
    )

    # Extract head-end pressure time series from result
    # result.time and result.pressure are the standard attributes
    sim_time = list(result.time)
    sim_pres = list(result.pressure)  # Pa (absolute)

    # --- Convert sim pressure to gauge MPa ---
    # Measurement is gauge relative to pre-event atmospheric baseline.
    # The simulation starts at P_cutoff=0.01 MPa; subtract atmospheric (0.101325 MPa).
    # However: the measurement is gauged to the pre-event pressure mean (9.844 psi = 0.0679 MPa gauge).
    # We use sim_gauge = sim_abs_MPa - 0.101325 MPa (standard atmosphere).
    P_atm_mpa = 0.101325
    sim_pres_gauge_mpa = [(p / 1e6) - P_atm_mpa for p in sim_pres]

    # --- Interpolate simulator onto measurement time grid ---
    sim_at_meas = [_interp(sim_time, sim_pres_gauge_mpa, t) for t in meas_win_t]

    # --- Compute shape metrics (window [0.0, 0.36] s) ---
    # RMSE and MAE (in MPa) — these are shape metrics on the absolute pressure scale.
    # Amplitude scoring is blocked; these values are not validation scores.
    residuals = [s - m for s, m in zip(sim_at_meas, meas_win_p)]
    rmse = math.sqrt(sum(r**2 for r in residuals) / len(residuals))
    mae = sum(abs(r) for r in residuals) / len(residuals)

    # Normalized shape metrics: divide by measurement peak in window
    meas_peak = max(abs(p) for p in meas_win_p)
    rmse_normalized = rmse / meas_peak if meas_peak > 0 else None
    mae_normalized = mae / meas_peak if meas_peak > 0 else None

    # Peak pressure comparison
    sim_win_t = [t for t in sim_time if win_lo <= t <= win_hi]
    sim_win_p = [p for t, p in zip(sim_time, sim_pres_gauge_mpa) if win_lo <= t <= win_hi]
    sim_peak = max(sim_win_p) if sim_win_p else None
    sim_peak_t = sim_win_t[sim_win_p.index(sim_peak)] if sim_win_p else None
    meas_peak_t = meas_win_t[meas_win_p.index(max(meas_win_p))]

    # --- Build result payload ---
    payload = {
        "schema_version": 1,
        "study": "zerox_startup_pressure_shape_comparison",
        "recorded_date": datetime.now(timezone.utc).date().isoformat(),
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_version": "v0.1",
        "status": "complete",
        "scope": {
            "case": "zerox",
            "window_s": list(STARTUP_WINDOW_S),
            "window_basis": (
                "Registered startup window from scoring_eligibility.case_windows.startup "
                "in validation_case_registry.json; set without model predictions."
            ),
            "shape_only": True,
            "amplitude_scoring_blocked": True,
            "amplitude_blocked_reason": (
                "Independent manufacturer pressure-calibration documentation absent; "
                "see ZEROX_TRANSDUCER_CALIBRATION_AUDIT.md and scoring_eligibility.absolute_pressure."
            ),
        },
        "motor": {
            "file": "motors/zerox.ric",
            "sha256": _sha256_file(MOTOR_PATH),
            "note": (
                "Baseline (un-LHS-calibrated) motor. The LHS-calibrated zerox_LHS.ric "
                "is intentionally excluded to prevent using fitted parameters."
            ),
        },
        "simulation_settings": NOMINAL_SETTINGS,
        "provenance": {
            "git_commit": _git_value("rev-parse", "HEAD"),
            "git_status": _git_value("status", "--short"),
            "platform": platform.platform(),
            "python": sys.version,
        },
        "measurement": {
            "source": "load_measurement('zerox') → static_fire_data/zerox_high_resolution.csv",
            "qualification": meas.get("qualification", ""),
            "sample_count_in_window": n_meas,
            "time_range_s": [meas_win_t[0], meas_win_t[-1]],
            "peak_mpa": max(meas_win_p),
            "peak_time_s": meas_peak_t,
        },
        "simulation": {
            "peak_gauge_mpa": sim_peak,
            "peak_time_s": sim_peak_t,
            "total_steps": len(sim_time),
        },
        "shape_metrics": {
            "window_s": list(STARTUP_WINDOW_S),
            "n_measurement_samples": n_meas,
            "rmse_mpa": rmse,
            "mae_mpa": mae,
            "rmse_normalized_by_meas_peak": rmse_normalized,
            "mae_normalized_by_meas_peak": mae_normalized,
            "peak_pressure_sim_mpa": sim_peak,
            "peak_pressure_meas_mpa": max(meas_win_p),
            "peak_time_sim_s": sim_peak_t,
            "peak_time_meas_s": meas_peak_t,
            "peak_pressure_ratio_sim_over_meas": (
                sim_peak / max(meas_win_p) if max(meas_win_p) > 0 else None
            ),
        },
        "claims": {
            "shape_comparison_complete": True,
            "amplitude_scoring": False,
            "formal_validation": False,
            "fit_parameters_used": False,
            "measurement_aligned_to_simulation": False,
        },
        "limitations": [
            "Only 28 measurement samples in [0.0, 0.36] s startup window (coarse ~10–33 ms intervals).",
            "Measurement amplitude is gauge from pre-event baseline; manufacturer calibration absent.",
            "Simulation gauge uses standard atmosphere offset (0.101325 MPa); mismatch possible.",
            "Baseline motor (zerox.ric) may differ from LHS-calibrated geometry; amplitude not meaningful.",
            "70 ms timing uncertainty in trigger transition bracket; time alignment unverified.",
            "Shape metrics are informational only; not a validation or fitness score.",
        ],
    }

    return payload


def main() -> int:
    payload = run_shape_comparison()
    OUTPUT_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"\nShape comparison saved to: {OUTPUT_PATH}")
    print(f"RMSE: {payload['shape_metrics']['rmse_mpa']:.4f} MPa")
    print(f"MAE:  {payload['shape_metrics']['mae_mpa']:.4f} MPa")
    print(f"RMSE (normalized): {payload['shape_metrics']['rmse_normalized_by_meas_peak']:.4f}")
    print(f"Sim peak: {payload['shape_metrics']['peak_pressure_sim_mpa']:.4f} MPa @ "
          f"t={payload['shape_metrics']['peak_time_sim_s']:.4f} s")
    print(f"Meas peak: {payload['shape_metrics']['peak_pressure_meas_mpa']:.4f} MPa @ "
          f"t={payload['shape_metrics']['peak_time_meas_s']:.4f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
