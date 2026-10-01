"""
Ignition-spike diagnostics for v0.7.0 startup investigations.

The functions here are post-processing only. They consume a simulation
``result`` plus optional geometry/propellant objects and estimate which
source family is active when the head-end pressure spike occurs.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:  # pragma: no cover - exercised only without matplotlib
    HAS_MATPLOTLIB = False


def stack_snapshots(result: dict, key: str) -> np.ndarray:
    """Return a ``(n_snapshots, n_cells)`` array for a snapshot key."""
    snapshots = result.get("snapshots", [])
    if not snapshots:
        return np.empty((0, 0))
    if key not in snapshots[0]:
        return np.empty((len(snapshots), 0))
    return np.array([np.asarray(s[key]) for s in snapshots])


def _snapshot_times(result: dict) -> np.ndarray:
    return np.array([float(s["t"]) for s in result.get("snapshots", [])])


def _infer_dx(result: dict, geo=None) -> float:
    if geo is not None and getattr(geo, "dx", None) is not None:
        return float(geo.dx)
    snapshots = result.get("snapshots", [])
    if not snapshots:
        return 1.0
    x = np.asarray(snapshots[0].get("x", []), dtype=float)
    if len(x) > 1:
        return float(np.median(np.diff(x)))
    return 1.0


def _rho_propellant(propellant=None) -> float:
    if propellant is None:
        return 1.0
    return float(getattr(propellant, "rho_propellant", 1.0))


def _first_time_at_fraction(times: np.ndarray, frac: np.ndarray,
                            threshold: float) -> float:
    idx = np.flatnonzero(frac >= threshold)
    if idx.size == 0:
        return float("nan")
    return float(times[int(idx[0])])


def _pyrogen_active_end_from_history(result: dict, active_fraction=0.01) -> float:
    t = np.asarray(result.get("time", []), dtype=float)
    mdot = np.asarray(result.get("mdot_ig", []), dtype=float)
    if t.size == 0 or mdot.size == 0:
        return float("nan")
    mdot_peak = float(np.max(mdot))
    threshold = max(1.0e-9, float(active_fraction) * mdot_peak)
    active = np.flatnonzero(mdot > threshold)
    if active.size == 0:
        return float("nan")
    return float(t[int(active[-1])])


def pressure_landmarks(result: dict, startup_margin_s=0.10,
                       fallback_startup_window_s=0.25) -> dict[str, float]:
    """Global pressure peak plus separately labeled startup-window peak."""
    t = np.asarray(result.get("time", []), dtype=float)
    p = np.asarray(result.get("P_head", []), dtype=float)
    if t.size == 0 or p.size == 0:
        return {
            "global_peak_time_s": float("nan"),
            "global_peak_pressure_pa": float("nan"),
            "global_peak_pressure_mpa": float("nan"),
            "startup_window_end_s": float("nan"),
            "startup_window_peak_time_s": float("nan"),
            "startup_window_peak_pressure_pa": float("nan"),
            "startup_window_peak_pressure_mpa": float("nan"),
            "takeoff_time_s": float("nan"),
            "peak_time_s": float("nan"),
            "peak_pressure_pa": float("nan"),
            "peak_pressure_mpa": float("nan"),
            "post_peak_trough_time_s": float("nan"),
            "post_peak_trough_pressure_pa": float("nan"),
        }

    global_peak_idx = int(np.argmax(p))
    pyrogen_active_end = _pyrogen_active_end_from_history(result)
    if np.isfinite(pyrogen_active_end):
        startup_window_end = pyrogen_active_end + float(startup_margin_s)
    else:
        startup_window_end = float(fallback_startup_window_s)
    startup_window_end = min(startup_window_end, float(t[-1]))

    startup_candidates = np.flatnonzero(t <= startup_window_end)
    if startup_candidates.size == 0:
        startup_peak_idx = global_peak_idx
    else:
        local = int(np.argmax(p[startup_candidates]))
        startup_peak_idx = int(startup_candidates[local])

    ambient = float(result.get("P_ambient", p[0]))
    threshold = ambient + 0.05 * max(float(p[startup_peak_idx]) - ambient, 0.0)
    takeoff_candidates = np.flatnonzero(p >= threshold)
    takeoff_time = float(t[int(takeoff_candidates[0])]) if takeoff_candidates.size else float("nan")

    post = p[startup_peak_idx:]
    trough_rel = int(np.argmin(post)) if post.size else 0
    trough_idx = startup_peak_idx + trough_rel
    return {
        "global_peak_time_s": float(t[global_peak_idx]),
        "global_peak_pressure_pa": float(p[global_peak_idx]),
        "global_peak_pressure_mpa": float(p[global_peak_idx] / 1.0e6),
        "startup_window_end_s": float(startup_window_end),
        "startup_window_peak_time_s": float(t[startup_peak_idx]),
        "startup_window_peak_pressure_pa": float(p[startup_peak_idx]),
        "startup_window_peak_pressure_mpa": float(p[startup_peak_idx] / 1.0e6),
        # Backward-compatible aliases now refer to the startup diagnostic peak.
        "takeoff_time_s": takeoff_time,
        "peak_time_s": float(t[startup_peak_idx]),
        "peak_pressure_pa": float(p[startup_peak_idx]),
        "peak_pressure_mpa": float(p[startup_peak_idx] / 1.0e6),
        "post_peak_trough_time_s": float(t[trough_idx]),
        "post_peak_trough_pressure_pa": float(p[trough_idx]),
        "post_peak_trough_pressure_mpa": float(p[trough_idx] / 1.0e6),
    }


def pyrogen_landmarks(result: dict, peak_time_s: float | None = None) -> dict[str, float | bool]:
    """Pyrogen pressure/mass-flow timing from time-history arrays."""
    t = np.asarray(result.get("time", []), dtype=float)
    mdot = np.asarray(result.get("mdot_ig", []), dtype=float)
    pig = np.asarray(result.get("P_ig", []), dtype=float)
    if t.size == 0 or mdot.size == 0:
        return {
            "mdot_peak_time_s": float("nan"),
            "mdot_peak_kg_s": float("nan"),
            "pyrogen_active_end_s": float("nan"),
            "P_ig_peak_time_s": float("nan"),
            "P_ig_peak_pa": float("nan"),
            "pyrogen_active_at_peak": False,
            "mdot_at_pressure_peak_kg_s": float("nan"),
        }

    mdot_peak = float(np.max(mdot))
    mdot_threshold = max(1.0e-9, 0.01 * mdot_peak)
    active = mdot > mdot_threshold
    mdot_peak_idx = int(np.argmax(mdot))
    active_end = float(t[int(np.flatnonzero(active)[-1])]) if np.any(active) else float("nan")
    pig_peak_idx = int(np.argmax(pig)) if pig.size else 0

    mdot_at_peak = float("nan")
    active_at_peak = False
    if peak_time_s is not None and np.isfinite(peak_time_s):
        mdot_at_peak = float(np.interp(peak_time_s, t, mdot))
        active_at_peak = bool(mdot_at_peak > mdot_threshold)

    return {
        "mdot_peak_time_s": float(t[mdot_peak_idx]),
        "mdot_peak_kg_s": mdot_peak,
        "mdot_peak_g_s": mdot_peak * 1000.0,
        "pyrogen_active_end_s": active_end,
        "P_ig_peak_time_s": float(t[pig_peak_idx]) if pig.size else float("nan"),
        "P_ig_peak_pa": float(pig[pig_peak_idx]) if pig.size else float("nan"),
        "P_ig_peak_mpa": float(pig[pig_peak_idx] / 1.0e6) if pig.size else float("nan"),
        "pyrogen_active_at_peak": active_at_peak,
        "mdot_at_pressure_peak_kg_s": mdot_at_peak,
        "mdot_at_pressure_peak_g_s": mdot_at_peak * 1000.0 if np.isfinite(mdot_at_peak) else float("nan"),
    }


def ignition_spread_metrics(result: dict) -> dict[str, Any]:
    """Burning-fraction time history and inferred first-burning times."""
    snapshot_times = _snapshot_times(result)
    burning = stack_snapshots(result, "is_burning").astype(bool)
    is_grain = stack_snapshots(result, "is_grain").astype(bool)

    ignition_times = np.asarray(result.get("ignition_time_by_cell", []), dtype=float)
    n_cells = 0
    if ignition_times.size:
        n_cells = int(ignition_times.size)
    elif burning.size:
        n_cells = int(burning.shape[1])

    if n_cells == 0:
        return {
            "times_s": snapshot_times,
            "burning_fraction": np.array([]),
            "first_ignition_time_s": float("nan"),
            "t10_s": float("nan"),
            "t50_s": float("nan"),
            "t90_s": float("nan"),
            "t100_s": float("nan"),
            "spread_10_90_s": float("nan"),
            "instant_ignition_collapse": False,
            "first_burning_time_by_cell_s": np.array([]),
            "axial_ignition_order": np.array([], dtype=int),
            "spread_metric_source": "none",
            "exact_spread_metrics": False,
        }

    if is_grain.size:
        grain_mask = np.any(is_grain, axis=0)
    else:
        grain_mask = np.ones(n_cells, dtype=bool)
    if grain_mask.size != n_cells:
        fixed = np.zeros(n_cells, dtype=bool)
        n = min(n_cells, grain_mask.size)
        fixed[:n] = grain_mask[:n]
        if n < n_cells:
            fixed[n:] = True
        grain_mask = fixed
    n_grain = max(int(np.sum(grain_mask)), 1)

    if burning.size and snapshot_times.size:
        burning_fraction = np.sum(burning[:, grain_mask], axis=1) / n_grain
        times = snapshot_times
    else:
        times = np.asarray(result.get("time", []), dtype=float)
        burning_fraction = np.zeros_like(times)

    exact_metrics = False
    first_by_cell = np.full(n_cells, np.nan)
    if ignition_times.size:
        n = min(n_cells, ignition_times.size)
        valid = np.isfinite(ignition_times[:n]) & (ignition_times[:n] < 1.0e9)
        valid_idx = np.flatnonzero(valid)
        first_by_cell[valid_idx] = ignition_times[:n][valid]
        exact_metrics = bool(np.any(valid & grain_mask[:n]))
    elif burning.size and snapshot_times.size:
        for i in range(n_cells):
            idx = np.flatnonzero(burning[:, i])
            if idx.size:
                first_by_cell[i] = snapshot_times[int(idx[0])]

    if exact_metrics:
        finite_grain_times = np.sort(first_by_cell[grain_mask & np.isfinite(first_by_cell)])
        if times.size:
            burning_fraction = np.array([
                np.sum(finite_grain_times <= t) / n_grain for t in times
            ], dtype=float)

        def exact_threshold_time(threshold: float) -> float:
            required = int(np.ceil(threshold * n_grain))
            required = min(max(required, 1), n_grain)
            if finite_grain_times.size < required:
                return float("nan")
            return float(finite_grain_times[required - 1])

        t10 = exact_threshold_time(0.10)
        t50 = exact_threshold_time(0.50)
        t90 = exact_threshold_time(0.90)
        t100 = exact_threshold_time(0.999)
        first_ignition = float(finite_grain_times[0]) if finite_grain_times.size else float("nan")
        metric_source = "ignition_time_by_cell"
    else:
        t10 = _first_time_at_fraction(times, burning_fraction, 0.10)
        t50 = _first_time_at_fraction(times, burning_fraction, 0.50)
        t90 = _first_time_at_fraction(times, burning_fraction, 0.90)
        t100 = _first_time_at_fraction(times, burning_fraction, 0.999)
        ignited = np.flatnonzero(burning_fraction > 0.0)
        first_ignition = float(times[int(ignited[0])]) if ignited.size else float("nan")
        metric_source = "snapshots" if times.size else "none"

    finite_cells = np.flatnonzero(np.isfinite(first_by_cell))
    order = finite_cells[np.argsort(first_by_cell[finite_cells])]

    if snapshot_times.size > 1:
        sample_dt = float(np.median(np.diff(snapshot_times)))
    elif times.size > 1:
        sample_dt = float(np.median(np.diff(times)))
    else:
        sample_dt = float("inf")
    spread_10_90 = t90 - t10 if np.isfinite(t10) and np.isfinite(t90) else float("nan")
    instant = bool(
        np.isfinite(spread_10_90)
        and np.isfinite(sample_dt)
        and spread_10_90 <= 2.0 * sample_dt
    )

    return {
        "times_s": times,
        "burning_fraction": burning_fraction,
        "first_ignition_time_s": first_ignition,
        "t10_s": t10,
        "t50_s": t50,
        "t90_s": t90,
        "t100_s": t100,
        "spread_10_90_s": spread_10_90,
        "instant_ignition_collapse": instant,
        "first_burning_time_by_cell_s": first_by_cell,
        "axial_ignition_order": order,
        "spread_metric_source": metric_source,
        "exact_spread_metrics": exact_metrics,
    }


def source_timeseries(result: dict, geo=None, propellant=None) -> dict[str, np.ndarray]:
    """Estimate source-family mass-flow rates at snapshot times."""
    times = _snapshot_times(result)
    if times.size == 0:
        empty = np.array([])
        return {
            "times_s": empty,
            "normal_sidewall_kg_s": empty,
            "erosive_sidewall_kg_s": empty,
            "total_sidewall_kg_s": empty,
            "endface_kg_s": empty,
            "pyrogen_kg_s": empty,
            "pyrogen_surface_heat_power_w": empty,
            "pyrogen_surface_heat_flux_w_m2": empty,
            "radiation_heat_power_w": empty,
            "radiation_heat_flux_w_m2": empty,
            "total_estimated_kg_s": empty,
            "normal_fraction": empty,
            "erosive_fraction": empty,
            "endface_fraction": empty,
            "pyrogen_fraction": empty,
        }

    dx = _infer_dx(result, geo)
    rho_p = _rho_propellant(propellant)
    c_burn = stack_snapshots(result, "C_burn")
    r_total = stack_snapshots(result, "r_total")
    r_erosive = stack_snapshots(result, "r_erosive")
    endface = stack_snapshots(result, "endface_msource")
    pyrogen_surface_heat_flux = stack_snapshots(result, "pyrogen_surface_heat_flux")
    radiation_heat_flux = stack_snapshots(result, "radiation_heat_flux")

    normal_rate = np.maximum(r_total - r_erosive, 0.0)
    normal = rho_p * np.sum(normal_rate * c_burn, axis=1) * dx
    erosive = rho_p * np.sum(np.maximum(r_erosive, 0.0) * c_burn, axis=1) * dx
    sidewall = rho_p * np.sum(np.maximum(r_total, 0.0) * c_burn, axis=1) * dx
    endface_total = np.sum(np.maximum(endface, 0.0), axis=1) * dx if endface.size else np.zeros_like(times)
    if pyrogen_surface_heat_flux.size:
        heat_flux = np.max(np.maximum(pyrogen_surface_heat_flux, 0.0), axis=1)
        heat_power = np.sum(
            np.maximum(pyrogen_surface_heat_flux, 0.0) * c_burn,
            axis=1,
        ) * dx
    else:
        heat_flux = np.zeros_like(times)
        heat_power = np.zeros_like(times)
    if radiation_heat_flux.size:
        rad_flux = np.max(np.maximum(radiation_heat_flux, 0.0), axis=1)
        rad_power = np.sum(
            np.maximum(radiation_heat_flux, 0.0) * c_burn,
            axis=1,
        ) * dx
    else:
        rad_flux = np.zeros_like(times)
        rad_power = np.zeros_like(times)

    hist_t = np.asarray(result.get("time", []), dtype=float)
    mdot = np.asarray(result.get("mdot_ig", []), dtype=float)
    if hist_t.size and mdot.size:
        pyrogen = np.interp(times, hist_t, mdot)
    else:
        pyrogen = np.zeros_like(times)

    total = normal + erosive + endface_total + pyrogen
    denom = np.where(total > 0.0, total, 1.0)
    return {
        "times_s": times,
        "normal_sidewall_kg_s": normal,
        "erosive_sidewall_kg_s": erosive,
        "total_sidewall_kg_s": sidewall,
        "endface_kg_s": endface_total,
        "pyrogen_kg_s": pyrogen,
        "pyrogen_surface_heat_power_w": heat_power,
        "pyrogen_surface_heat_flux_w_m2": heat_flux,
        "radiation_heat_power_w": rad_power,
        "radiation_heat_flux_w_m2": rad_flux,
        "total_estimated_kg_s": total,
        "normal_fraction": normal / denom,
        "erosive_fraction": erosive / denom,
        "endface_fraction": endface_total / denom,
        "pyrogen_fraction": pyrogen / denom,
    }


def energy_momentum_timeseries(result: dict) -> dict[str, np.ndarray]:
    """Return per-step energy and pyrogen momentum audit histories.

    Energy residual convention:
    ``gas_sensible_dE_dt - convective_scalar_flux_power
    - thermal_source_power - clipping_correction_power``.
    Positive convective power enters the gas control volume.
    """
    times = np.asarray(result.get("time", []), dtype=float)
    keys = (
        "gas_sensible_energy_before",
        "gas_sensible_energy",
        "gas_sensible_dE_dt",
        "normal_sidewall_thermal_power",
        "erosive_sidewall_thermal_power",
        "endface_thermal_power",
        "pyrogen_gas_thermal_power",
        "pyrogen_enthalpy_power",
        "pyrogen_surface_heat_power",
        "gas_surface_heat_sink_power",
        "radiation_heat_power",
        "radiation_sink_power",
        "convective_wall_heat_sink_power",
        "convective_scalar_flux_power",
        "nozzle_scalar_flux_power",
        "nozzle_enthalpy_power",
        "thermal_source_power",
        "clipping_correction_power",
        "energy_residual",
        "pyrogen_momentum_expected",
        "pyrogen_momentum_deposited",
        "pyrogen_momentum_residual",
    )
    out = {"times_s": times}
    for key in keys:
        values = np.asarray(result.get(key, np.zeros_like(times)), dtype=float)
        if values.size != times.size:
            values = np.zeros_like(times)
        out[key] = values
    return out


def sensible_enthalpy_audit(result: dict) -> dict[str, np.ndarray]:
    """Reconstruct the model's transported ``Cp*T`` ledger.

    The solver does not carry thermodynamic total energy. The between-step
    storage jump is separate because mixture properties and geometry may be
    refreshed between one PISO storage sample and the next.
    """
    energy = energy_momentum_timeseries(result)
    source_reconstructed = (
        energy["normal_sidewall_thermal_power"]
        + energy["erosive_sidewall_thermal_power"]
        + energy["endface_thermal_power"]
        + energy["pyrogen_enthalpy_power"]
        - energy["gas_surface_heat_sink_power"]
        - energy["radiation_sink_power"]
        - energy["convective_wall_heat_sink_power"]
    )
    storage_jump = np.zeros_like(energy["times_s"])
    if storage_jump.size > 1:
        storage_jump[1:] = (
            energy["gas_sensible_energy_before"][1:]
            - energy["gas_sensible_energy"][:-1]
        )
    return {
        **energy,
        "source_reconstructed_power": source_reconstructed,
        "source_ledger_mismatch_power": (
            energy["thermal_source_power"] - source_reconstructed
        ),
        "nozzle_boundary_mismatch_power": (
            energy["convective_scalar_flux_power"]
            + energy["nozzle_enthalpy_power"]
        ),
        "between_step_storage_jump": storage_jump,
    }


def step_diagnostics_timeseries(result: dict) -> dict[str, np.ndarray]:
    """Return per-step diagnostic histories used for early-failure probes."""
    times = np.asarray(result.get("time", []), dtype=float)
    keys = (
        "dt",
        "P_head",
        "P_exit",
        "massflow",
        "n_burning",
        "n_ignited",
        "radiation_emitter_count",
        "radiation_receiver_count",
        "min_gas_temperature",
        "max_gas_temperature",
        "min_surface_temperature",
        "max_surface_temperature",
        "min_pressure",
        "max_pressure",
        "max_mach",
    )
    out = {"times_s": times}
    for key in keys:
        values = np.asarray(result.get(key, np.zeros_like(times)), dtype=float)
        if values.size != times.size:
            values = np.zeros_like(times)
        out[key] = values
    return out


def _time_integral(values: np.ndarray, dt: np.ndarray, mask: np.ndarray) -> float:
    if values.size == 0 or dt.size == 0 or mask.size == 0:
        return 0.0
    n = min(values.size, dt.size, mask.size)
    if n <= 0:
        return 0.0
    return float(np.sum(values[:n][mask[:n]] * dt[:n][mask[:n]]))


def _max_abs(values: np.ndarray, mask: np.ndarray | None = None) -> float:
    if values.size == 0:
        return float("nan")
    if mask is None:
        selected = values
    else:
        n = min(values.size, mask.size)
        selected = values[:n][mask[:n]]
    if selected.size == 0:
        return float("nan")
    return float(np.max(np.abs(selected)))


def _max_value(values: np.ndarray, mask: np.ndarray | None = None) -> float:
    if values.size == 0:
        return float("nan")
    if mask is None:
        selected = values
    else:
        n = min(values.size, mask.size)
        selected = values[:n][mask[:n]]
    if selected.size == 0:
        return float("nan")
    return float(np.max(selected))


def _min_value(values: np.ndarray, mask: np.ndarray | None = None) -> float:
    if values.size == 0:
        return float("nan")
    if mask is None:
        selected = values
    else:
        n = min(values.size, mask.size)
        selected = values[:n][mask[:n]]
    if selected.size == 0:
        return float("nan")
    return float(np.min(selected))


def _first_pressure_time(times: np.ndarray, pressure: np.ndarray,
                         threshold_pa: float) -> float:
    if times.size == 0 or pressure.size == 0:
        return float("nan")
    n = min(times.size, pressure.size)
    idx = np.flatnonzero(pressure[:n] >= threshold_pa)
    if idx.size == 0:
        return float("nan")
    return float(times[int(idx[0])])


def _first_positive_time(times: np.ndarray, values: np.ndarray,
                         threshold: float) -> float:
    if times.size == 0 or values.size == 0:
        return float("nan")
    n = min(times.size, values.size)
    idx = np.flatnonzero(values[:n] > threshold)
    if idx.size == 0:
        return float("nan")
    return float(times[int(idx[0])])


def _first_sustained_time(times: np.ndarray, values: np.ndarray,
                          threshold: float, count: int = 3) -> float:
    if times.size == 0 or values.size == 0 or count <= 0:
        return float("nan")
    n = min(times.size, values.size)
    run = 0
    for i in range(n):
        if values[i] > threshold:
            run += 1
            if run >= count:
                return float(times[i - count + 1])
        else:
            run = 0
    return float("nan")


def _first_threshold_time(times: np.ndarray, values: np.ndarray,
                          threshold: float, less_than: bool = False) -> float:
    if times.size == 0 or values.size == 0:
        return float("nan")
    n = min(times.size, values.size)
    if less_than:
        idx = np.flatnonzero(values[:n] < threshold)
    else:
        idx = np.flatnonzero(values[:n] > threshold)
    if idx.size == 0:
        return float("nan")
    return float(times[int(idx[0])])


def _nearest_index(times: np.ndarray, event_time: float) -> int:
    if times.size == 0 or not np.isfinite(event_time):
        return -1
    return int(np.argmin(np.abs(times - event_time)))


def _nearest_snapshot(result: dict, event_time: float) -> tuple[int, dict[str, Any] | None]:
    snapshots = result.get("snapshots", [])
    if not snapshots or not np.isfinite(event_time):
        return -1, None
    snap_times = _snapshot_times(result)
    idx = _nearest_index(snap_times, event_time)
    if idx < 0:
        return -1, None
    return idx, snapshots[idx]


def _step_flag_values(result: dict, event_time: float,
                      clipping_step: np.ndarray) -> dict[str, Any]:
    times = np.asarray(result.get("time", []), dtype=float)
    idx = _nearest_index(times, event_time)
    if idx < 0:
        return {
            "history_time_s": float("nan"),
            "dt_below_1e_8": False,
            "pressure_above_100mpa": False,
            "mach_above_1e3": False,
            "clipping_dominated_step": False,
        }
    dt = np.asarray(result.get("dt", []), dtype=float)
    max_pressure = np.asarray(result.get("max_pressure", []), dtype=float)
    max_mach = np.asarray(result.get("max_mach", []), dtype=float)
    return {
        "history_time_s": float(times[idx]),
        "dt_below_1e_8": bool(dt.size > idx and dt[idx] < 1.0e-8),
        "pressure_above_100mpa": bool(max_pressure.size > idx and max_pressure[idx] > 100.0e6),
        "mach_above_1e3": bool(max_mach.size > idx and max_mach[idx] > 1.0e3),
        "clipping_dominated_step": bool(clipping_step.size > idx and clipping_step[idx]),
    }


def early_time_diagnostics(result: dict, window_s: float = 0.030,
                           spread: dict[str, Any] | None = None) -> dict[str, Any]:
    """Summarize early-time termination, radiation, clipping, and pressure rise."""
    times = np.asarray(result.get("time", []), dtype=float)
    pressure = np.asarray(result.get("P_head", []), dtype=float)
    dt = np.asarray(result.get("dt", []), dtype=float)
    if dt.size != times.size:
        if times.size > 1:
            dt = np.diff(np.r_[times[0], times])
            if dt.size:
                dt[0] = dt[1] if dt.size > 1 else 0.0
        else:
            dt = np.zeros_like(times)

    if times.size:
        mask = times <= min(float(window_s), float(times[-1]))
    else:
        mask = np.array([], dtype=bool)

    summary = result.get("summary", {})
    termination = str(summary.get("termination", "unknown"))
    history_cap_reached = bool(
        summary.get("history_cap_reached", False)
        or termination == "history array full"
    )

    radiation_enabled = float(summary.get("active_radiation_emissivity", 0.0)) > 0.0
    emitter_count = np.asarray(result.get("radiation_emitter_count", []), dtype=float)
    receiver_count = np.asarray(result.get("radiation_receiver_count", []), dtype=float)
    radiation_heat_power = np.asarray(result.get("radiation_heat_power", []), dtype=float)
    radiation_sink_power = np.asarray(result.get("radiation_sink_power", []), dtype=float)

    max_emitters = _max_value(emitter_count, mask)
    max_receivers = _max_value(receiver_count, mask)
    max_radiation_power = _max_value(radiation_heat_power, mask)
    if not radiation_enabled:
        max_emitters = 0.0
        max_receivers = 0.0
        max_radiation_power = 0.0
    first_radiation_time = _first_positive_time(
        times, np.maximum(radiation_heat_power, receiver_count), 0.0
    )
    radiation_enabled_zero_activity = bool(
        radiation_enabled
        and (not np.isfinite(max_emitters) or max_emitters <= 0.0)
        and (not np.isfinite(max_receivers) or max_receivers <= 0.0)
        and (not np.isfinite(max_radiation_power) or max_radiation_power <= 0.0)
    )

    clipping_power = np.asarray(result.get("clipping_correction_power", []), dtype=float)
    thermal_power = np.asarray(result.get("thermal_source_power", []), dtype=float)
    convective_power = np.asarray(result.get("convective_scalar_flux_power", []), dtype=float)
    clipping_energy = _time_integral(clipping_power, dt, mask)
    thermal_energy = _time_integral(thermal_power, dt, mask)
    convective_energy = _time_integral(convective_power, dt, mask)
    radiation_heat_energy = _time_integral(radiation_heat_power, dt, mask)
    radiation_sink_energy = _time_integral(radiation_sink_power, dt, mask)
    energy_scale = max(abs(thermal_energy), abs(convective_energy), 1.0e-12)
    clipping_dominated = bool(abs(clipping_energy) >= 0.5 * energy_scale)
    n_step = min(clipping_power.size, thermal_power.size, convective_power.size)
    if n_step:
        step_scale = np.maximum.reduce([
            np.abs(thermal_power[:n_step]),
            np.abs(convective_power[:n_step]),
            np.ones(n_step) * 1.0e-12,
        ])
        clipping_step = np.abs(clipping_power[:n_step]) >= 0.5 * step_scale
    else:
        clipping_step = np.array([], dtype=bool)
    radiation_energy_mismatch = bool(
        max(abs(radiation_heat_energy), abs(radiation_sink_energy)) > 1.0e-9
        and abs(radiation_heat_energy - radiation_sink_energy)
        > 1.0e-3 * max(abs(radiation_heat_energy), abs(radiation_sink_energy), 1.0)
    )

    ignition_times = np.asarray(result.get("ignition_time_by_cell", []), dtype=float)
    finite_ignition = np.flatnonzero(ignition_times < 1.0e9)
    if finite_ignition.size:
        local = int(np.argmin(ignition_times[finite_ignition]))
        first_ignition_cell = int(finite_ignition[local])
        first_ignition_time = float(ignition_times[first_ignition_cell])
    else:
        first_ignition_cell = int(summary.get("first_ignition_cell", -1))
        first_ignition_time = float(summary.get("first_ignition_time_s", float("nan")))

    erosive_power = np.asarray(result.get("erosive_sidewall_thermal_power", []), dtype=float)
    erosive_threshold = 0.01 * _max_value(erosive_power)
    if not np.isfinite(erosive_threshold) or erosive_threshold <= 0.0:
        erosive_threshold = 0.0
    first_erosive_time = _first_sustained_time(times, erosive_power, erosive_threshold)

    if spread is None:
        spread = ignition_spread_metrics(result)

    max_pressure_series = np.asarray(result.get("max_pressure", []), dtype=float)
    max_mach_series = np.asarray(result.get("max_mach", []), dtype=float)
    first_dt_collapse_time = _first_threshold_time(times, dt, 1.0e-8, less_than=True)
    first_pressure_collapse_time = _first_threshold_time(
        times, max_pressure_series, 100.0e6
    )
    first_mach_collapse_time = _first_threshold_time(times, max_mach_series, 1.0e3)
    first_clipping_dominated_time = _first_threshold_time(
        times[:clipping_step.size], clipping_step.astype(float), 0.5
    )
    first_history_cap_time = float(times[-1]) if history_cap_reached and times.size else float("nan")

    energy_residual = np.asarray(result.get("energy_residual", []), dtype=float)
    momentum_residual = np.asarray(result.get("pyrogen_momentum_residual", []), dtype=float)
    max_energy_residual = _max_abs(energy_residual, mask)
    max_momentum_residual = _max_abs(momentum_residual, mask)
    residual_scale = max(
        _max_abs(thermal_power, mask),
        _max_abs(convective_power, mask),
        _max_abs(clipping_power, mask),
        1.0,
    )
    energy_residual_relative = (
        max_energy_residual / residual_scale if np.isfinite(max_energy_residual) else float("nan")
    )

    # A runtime termination_code == 4 (numerical-collapse abort) is
    # authoritative: the trip's per-step thresholds (Mach > 100,
    # dt < 1e-9, P > 1 GPa) are stricter than the in-window classifier
    # thresholds (Mach > 1000, dt < 1e-8, P > 100 MPa), so a clean
    # abort may not trip any in-window signal. Promote to collapse so
    # collapse_class agrees with diagnostic_failure_mode.
    trip_aborted = int(summary.get("termination_code", -1)) == 4
    collapse_detected = bool(
        history_cap_reached
        or trip_aborted
        or np.isfinite(first_dt_collapse_time)
        or np.isfinite(first_pressure_collapse_time)
        or np.isfinite(first_mach_collapse_time)
        or clipping_dominated
    )

    exact_t10 = float(spread.get("t10_s", float("nan")))
    exact_t90 = float(spread.get("t90_s", float("nan")))
    residuals_closed = bool(
        np.isfinite(energy_residual_relative)
        and energy_residual_relative < 1.0e-6
        and (not np.isfinite(max_momentum_residual) or max_momentum_residual < 1.0e-6)
    )
    if collapse_detected:
        collapse_class = "collapse"
    elif np.isfinite(exact_t10) and np.isfinite(exact_t90) and residuals_closed:
        collapse_class = "stable"
    else:
        collapse_class = "borderline"

    first_collapse_times = [
        value for value in (
            first_dt_collapse_time,
            first_pressure_collapse_time,
            first_mach_collapse_time,
            first_clipping_dominated_time,
            first_history_cap_time,
        )
        if np.isfinite(value)
    ]
    first_collapse_time = min(first_collapse_times) if first_collapse_times else float("nan")

    collapse_branch_suspect = "none"
    if radiation_energy_mismatch:
        collapse_branch_suspect = "radiation_source_sink_accounting"
    elif collapse_detected and radiation_enabled and np.isfinite(first_erosive_time):
        collapse_branch_suspect = "piso_nozzle_front_numerical_instability"
    elif collapse_detected:
        collapse_branch_suspect = "grid_or_timestep_sensitivity"

    if trip_aborted:
        failure_mode = "numerical_collapse_aborted"
    elif history_cap_reached and clipping_dominated:
        failure_mode = "timestep_front_numerical_pathology"
    elif history_cap_reached:
        failure_mode = "history_cap_or_step_limit"
    elif radiation_energy_mismatch:
        failure_mode = "radiation_energy_accounting_suspect"
    elif radiation_enabled_zero_activity:
        failure_mode = "radiation_enabled_zero_activity"
    elif not np.isfinite(first_ignition_time):
        failure_mode = "physically_weak_ignition"
    elif clipping_dominated:
        failure_mode = "timestep_front_numerical_pathology"
    else:
        failure_mode = "normal_completed"

    return {
        "window_s": float(window_s),
        "termination": termination,
        "termination_code": int(summary.get("termination_code", -1)),
        "history_cap_reached": history_cap_reached,
        "history_capacity": int(summary.get("history_capacity", 0)),
        "steps": int(summary.get("steps", len(times))),
        "final_time_s": float(times[-1]) if times.size else float("nan"),
        "dt_min_s": _min_value(dt, mask),
        "dt_median_s": float(np.median(dt[mask])) if dt.size and np.any(mask) else float("nan"),
        "dt_final_s": float(dt[-1]) if dt.size else float("nan"),
        "first_dt_collapse_time_s": first_dt_collapse_time,
        "first_pressure_collapse_time_s": first_pressure_collapse_time,
        "first_mach_collapse_time_s": first_mach_collapse_time,
        "first_clipping_dominated_time_s": first_clipping_dominated_time,
        "first_history_cap_time_s": first_history_cap_time,
        "first_collapse_time_s": first_collapse_time,
        "pressure_time_1_mpa_s": _first_pressure_time(times, pressure, 1.0e6),
        "pressure_time_2_mpa_s": _first_pressure_time(times, pressure, 2.0e6),
        "pressure_time_5_mpa_s": _first_pressure_time(times, pressure, 5.0e6),
        "pressure_time_10_mpa_s": _first_pressure_time(times, pressure, 10.0e6),
        "first_ignition_time_s": first_ignition_time,
        "first_ignition_cell": first_ignition_cell,
        "first_sustained_erosive_time_s": first_erosive_time,
        "first_full_grain_ignition_time_s": float(spread.get("t100_s", float("nan"))),
        "max_burning_cells": _max_value(np.asarray(result.get("n_burning", []), dtype=float), mask),
        "final_burning_cells": float(np.asarray(result.get("n_burning", [np.nan]), dtype=float)[-1])
        if len(result.get("n_burning", [])) else float("nan"),
        "max_ignited_cells": _max_value(np.asarray(result.get("n_ignited", []), dtype=float), mask),
        "final_ignited_cells": float(np.asarray(result.get("n_ignited", [np.nan]), dtype=float)[-1])
        if len(result.get("n_ignited", [])) else float("nan"),
        "radiation_enabled": radiation_enabled,
        "active_radiation_emissivity": float(summary.get("active_radiation_emissivity", 0.0)),
        "max_radiation_emitters": max_emitters,
        "max_radiation_receivers": max_receivers,
        "first_radiation_time_s": first_radiation_time,
        "radiation_enabled_zero_activity": radiation_enabled_zero_activity,
        "radiation_energy_mismatch": radiation_energy_mismatch,
        "max_nozzle_massflow_kg_s": _max_value(np.asarray(result.get("massflow", []), dtype=float), mask),
        "max_clipping_power_abs_w": _max_abs(clipping_power, mask),
        "max_convective_scalar_flux_power_abs_w": _max_abs(convective_power, mask),
        "max_thermal_source_power_abs_w": _max_abs(thermal_power, mask),
        "max_energy_residual_abs_w": max_energy_residual,
        "energy_residual_relative": energy_residual_relative,
        "max_momentum_residual_abs_n": max_momentum_residual,
        "clipping_energy_j": clipping_energy,
        "thermal_source_energy_j": thermal_energy,
        "convective_scalar_flux_energy_j": convective_energy,
        "radiation_heat_energy_j": radiation_heat_energy,
        "radiation_sink_energy_j": radiation_sink_energy,
        "clipping_dominated_energy": clipping_dominated,
        "min_gas_temperature_k": _min_value(np.asarray(result.get("min_gas_temperature", []), dtype=float), mask),
        "max_gas_temperature_k": _max_value(np.asarray(result.get("max_gas_temperature", []), dtype=float), mask),
        "min_surface_temperature_k": _min_value(np.asarray(result.get("min_surface_temperature", []), dtype=float), mask),
        "max_surface_temperature_k": _max_value(np.asarray(result.get("max_surface_temperature", []), dtype=float), mask),
        "min_pressure_pa": _min_value(np.asarray(result.get("min_pressure", []), dtype=float), mask),
        "max_pressure_pa": _max_value(np.asarray(result.get("max_pressure", []), dtype=float), mask),
        "max_mach": _max_value(np.asarray(result.get("max_mach", []), dtype=float), mask),
        "collapse_detected": collapse_detected,
        "collapse_class": collapse_class,
        "collapse_branch_suspect": collapse_branch_suspect,
        "residuals_closed": residuals_closed,
        "diagnostic_failure_mode": failure_mode,
    }


def _event_cell_from_snapshot(snapshot: dict[str, Any], preference: str,
                              fallback_cell: int = -1) -> int:
    if snapshot is None:
        return fallback_cell
    if fallback_cell >= 0:
        return int(fallback_cell)
    if preference == "radiation":
        values = np.asarray(snapshot.get("radiation_heat_flux", []), dtype=float)
        if values.size and np.max(values) > 0.0:
            return int(np.argmax(values))
    if preference == "pressure":
        values = np.asarray(snapshot.get("P", []), dtype=float)
        if values.size:
            return int(np.argmax(values))
    values = np.asarray(snapshot.get("Mach", []), dtype=float)
    if values.size:
        return int(np.argmax(np.abs(values)))
    return fallback_cell


def collapse_event_trace(result: dict, early: dict[str, Any],
                         spread: dict[str, Any]) -> list[dict[str, Any]]:
    """Compact per-cell snapshot trace around ignition/radiation/collapse events."""
    snapshots = result.get("snapshots", [])
    if not snapshots:
        return []

    times = np.asarray(result.get("time", []), dtype=float)
    clipping_power = np.asarray(result.get("clipping_correction_power", []), dtype=float)
    thermal_power = np.asarray(result.get("thermal_source_power", []), dtype=float)
    convective_power = np.asarray(result.get("convective_scalar_flux_power", []), dtype=float)
    n_step = min(clipping_power.size, thermal_power.size, convective_power.size)
    if n_step:
        step_scale = np.maximum.reduce([
            np.abs(thermal_power[:n_step]),
            np.abs(convective_power[:n_step]),
            np.ones(n_step) * 1.0e-12,
        ])
        clipping_step = np.abs(clipping_power[:n_step]) >= 0.5 * step_scale
    else:
        clipping_step = np.array([], dtype=bool)

    events: list[tuple[str, float, str, int]] = []
    first_cell = int(early.get("first_ignition_cell", -1))
    events.append((
        "first_ignition",
        float(early.get("first_ignition_time_s", float("nan"))),
        "pressure",
        first_cell,
    ))
    events.append((
        "first_radiation",
        float(early.get("first_radiation_time_s", float("nan"))),
        "radiation",
        -1,
    ))
    events.append((
        "first_dt_below_1e-8",
        float(early.get("first_dt_collapse_time_s", float("nan"))),
        "mach",
        -1,
    ))
    events.append((
        "first_pressure_above_100mpa",
        float(early.get("first_pressure_collapse_time_s", float("nan"))),
        "pressure",
        -1,
    ))
    events.append((
        "first_mach_above_1e3",
        float(early.get("first_mach_collapse_time_s", float("nan"))),
        "mach",
        -1,
    ))
    events.append((
        "first_clipping_dominated",
        float(early.get("first_clipping_dominated_time_s", float("nan"))),
        "mach",
        -1,
    ))
    events.append((
        "history_cap",
        float(early.get("first_history_cap_time_s", float("nan"))),
        "mach",
        -1,
    ))

    rows: list[dict[str, Any]] = []
    for event_name, event_time, preference, fallback_cell in events:
        if not np.isfinite(event_time):
            continue
        snap_idx, snapshot = _nearest_snapshot(result, event_time)
        if snapshot is None:
            continue
        event_cell = _event_cell_from_snapshot(snapshot, preference, fallback_cell)
        p = np.asarray(snapshot.get("P", []), dtype=float)
        n_cells = int(p.size)
        if n_cells == 0:
            continue
        if event_cell < 0 or event_cell >= n_cells:
            event_cell = int(np.argmax(p))
        start = max(0, event_cell - 3)
        stop = min(n_cells, event_cell + 4)
        flags = _step_flag_values(result, event_time, clipping_step)
        snap_time = float(snapshot.get("t", float("nan")))
        arrays = {
            "P": p,
            "T": np.asarray(snapshot.get("T", np.full(n_cells, np.nan)), dtype=float),
            "T_surf": np.asarray(snapshot.get("T_surf", np.full(n_cells, np.nan)), dtype=float),
            "Mach": np.asarray(snapshot.get("Mach", np.full(n_cells, np.nan)), dtype=float),
            "r_total": np.asarray(snapshot.get("r_total", np.full(n_cells, np.nan)), dtype=float),
            "r_erosive": np.asarray(snapshot.get("r_erosive", np.full(n_cells, np.nan)), dtype=float),
            "radiation_heat_flux": np.asarray(
                snapshot.get("radiation_heat_flux", np.full(n_cells, np.nan)), dtype=float
            ),
            "mass_source": np.asarray(snapshot.get("mass_source", np.full(n_cells, np.nan)), dtype=float),
            "thermal_source": np.asarray(snapshot.get("thermal_source", np.full(n_cells, np.nan)), dtype=float),
        }
        for cell in range(start, stop):
            rows.append({
                "event": event_name,
                "event_time_s": event_time,
                "snapshot_index": snap_idx,
                "snapshot_time_s": snap_time,
                "history_time_s": flags["history_time_s"],
                "event_cell": event_cell,
                "cell": cell,
                "cell_offset": cell - event_cell,
                "P_pa": arrays["P"][cell],
                "T_k": arrays["T"][cell],
                "T_surf_k": arrays["T_surf"][cell],
                "Mach": arrays["Mach"][cell],
                "r_total_m_s": arrays["r_total"][cell],
                "r_erosive_m_s": arrays["r_erosive"][cell],
                "radiation_heat_flux_w_m2": arrays["radiation_heat_flux"][cell],
                "mass_source_kg_m_s": arrays["mass_source"][cell],
                "thermal_source_k_kg_m_s": arrays["thermal_source"][cell],
                "dt_below_1e_8": flags["dt_below_1e_8"],
                "pressure_above_100mpa": flags["pressure_above_100mpa"],
                "mach_above_1e3": flags["mach_above_1e3"],
                "clipping_dominated_step": flags["clipping_dominated_step"],
                "history_cap_reached": bool(early.get("history_cap_reached", False)),
                "spread_metric_source": spread.get("spread_metric_source", "unknown"),
            })
    return rows


def _nearest_source_at_peak(sources: dict[str, np.ndarray], peak_time_s: float) -> dict[str, float]:
    times = sources.get("times_s", np.array([]))
    if times.size == 0 or not np.isfinite(peak_time_s):
        return {
            "normal_fraction_at_peak": float("nan"),
            "erosive_fraction_at_peak": float("nan"),
            "endface_fraction_at_peak": float("nan"),
            "pyrogen_fraction_at_peak": float("nan"),
            "source_snapshot_time_s": float("nan"),
        }
    idx = int(np.argmin(np.abs(times - peak_time_s)))
    return {
        "normal_fraction_at_peak": float(sources["normal_fraction"][idx]),
        "erosive_fraction_at_peak": float(sources["erosive_fraction"][idx]),
        "endface_fraction_at_peak": float(sources["endface_fraction"][idx]),
        "pyrogen_fraction_at_peak": float(sources["pyrogen_fraction"][idx]),
        "source_snapshot_time_s": float(times[idx]),
    }


def classify_driver(pressure: dict[str, float], pyrogen: dict[str, Any],
                    spread: dict[str, Any], sources_at_peak: dict[str, float]) -> dict[str, Any]:
    """Classify the most likely pressure-spike driver."""
    pyrogen_active = bool(pyrogen.get("pyrogen_active_at_peak", False))
    pyro_frac = float(sources_at_peak.get("pyrogen_fraction_at_peak", np.nan))
    erosive_frac = float(sources_at_peak.get("erosive_fraction_at_peak", np.nan))
    endface_frac = float(sources_at_peak.get("endface_fraction_at_peak", np.nan))
    instant = bool(spread.get("instant_ignition_collapse", False))

    flags = {
        "pyrogen_combustion": pyrogen_active and (not np.isfinite(pyro_frac) or pyro_frac >= 0.25),
        "erosive_snap_on": np.isfinite(erosive_frac) and erosive_frac >= 0.35 and (
            not pyrogen_active or erosive_frac >= pyro_frac
        ),
        "endface_sources": np.isfinite(endface_frac) and endface_frac >= 0.35 and (
            not pyrogen_active or endface_frac >= pyro_frac
        ),
        "immediate_goodman_spread": instant,
    }

    if flags["pyrogen_combustion"]:
        primary = "pyrogen_combustion"
    elif flags["erosive_snap_on"]:
        primary = "erosive_snap_on"
    elif flags["endface_sources"]:
        primary = "endface_sources"
    elif flags["immediate_goodman_spread"]:
        primary = "immediate_goodman_spread"
    else:
        primary = "mixed_or_unclassified"

    return {
        "primary_driver": primary,
        **flags,
        "peak_time_s": pressure.get("peak_time_s", float("nan")),
    }


def analyze_ignition_spike(result: dict, geo=None, propellant=None) -> dict[str, Any]:
    """Run the full diagnostic reduction for one simulation result."""
    pressure = pressure_landmarks(result)
    pyrogen = pyrogen_landmarks(result, pressure["peak_time_s"])
    spread = ignition_spread_metrics(result)
    sources = source_timeseries(result, geo=geo, propellant=propellant)
    energy = energy_momentum_timeseries(result)
    step_diagnostics = step_diagnostics_timeseries(result)
    early = early_time_diagnostics(result, spread=spread)
    collapse_trace = collapse_event_trace(result, early, spread)
    sources_at_peak = _nearest_source_at_peak(sources, pressure["peak_time_s"])
    classification = classify_driver(pressure, pyrogen, spread, sources_at_peak)
    classification["diagnostic_failure_mode"] = early["diagnostic_failure_mode"]
    classification["history_cap_reached"] = early["history_cap_reached"]
    classification["radiation_enabled_zero_activity"] = early["radiation_enabled_zero_activity"]
    classification["clipping_dominated_energy"] = early["clipping_dominated_energy"]
    classification["radiation_energy_mismatch"] = early["radiation_energy_mismatch"]
    return {
        "pressure": pressure,
        "pyrogen": pyrogen,
        "ignition_spread": spread,
        "sources": sources,
        "energy": energy,
        "step_diagnostics": step_diagnostics,
        "early": early,
        "collapse_trace": collapse_trace,
        "sources_at_peak": sources_at_peak,
        "classification": classification,
    }


def _flatten_scalars(prefix: str, data: dict[str, Any], rows: list[tuple[str, Any]]) -> None:
    for key, value in data.items():
        name = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            _flatten_scalars(name, value, rows)
        elif np.isscalar(value) or isinstance(value, (str, bool)):
            rows.append((name, value))


def classification_report(diagnostics: dict[str, Any]) -> str:
    """Compact human-readable report for one diagnostic run."""
    pressure = diagnostics["pressure"]
    pyrogen = diagnostics["pyrogen"]
    spread = diagnostics["ignition_spread"]
    sources = diagnostics["sources_at_peak"]
    source_series = diagnostics["sources"]
    early = diagnostics.get("early", {})
    cls = diagnostics["classification"]
    heat_flux_peak = float("nan")
    heat_power_peak = float("nan")
    rad_flux_peak = float("nan")
    rad_power_peak = float("nan")
    if source_series["times_s"].size:
        heat_flux_peak = float(np.max(source_series["pyrogen_surface_heat_flux_w_m2"]))
        heat_power_peak = float(np.max(source_series["pyrogen_surface_heat_power_w"]))
        rad_flux_peak = float(np.max(source_series["radiation_heat_flux_w_m2"]))
        rad_power_peak = float(np.max(source_series["radiation_heat_power_w"]))
    lines = [
        f"primary_driver: {cls['primary_driver']}",
        f"startup_window_peak: {pressure['startup_window_peak_pressure_mpa']:.3g} MPa "
        f"at {pressure['startup_window_peak_time_s']:.6g} s "
        f"(window_end={pressure['startup_window_end_s']:.6g} s)",
        f"global_peak: {pressure['global_peak_pressure_mpa']:.3g} MPa "
        f"at {pressure['global_peak_time_s']:.6g} s",
        f"pyrogen_active_at_peak: {pyrogen['pyrogen_active_at_peak']} "
        f"(mdot={pyrogen['mdot_at_pressure_peak_g_s']:.3g} g/s)",
        f"ignition_spread_10_90: {spread['spread_10_90_s']:.6g} s "
        f"(instant={spread['instant_ignition_collapse']})",
        "source_fractions_at_peak: "
        f"normal={sources['normal_fraction_at_peak']:.3g}, "
        f"erosive={sources['erosive_fraction_at_peak']:.3g}, "
        f"endface={sources['endface_fraction_at_peak']:.3g}, "
        f"pyrogen={sources['pyrogen_fraction_at_peak']:.3g}",
        f"pyrogen_surface_heat_peak: {heat_flux_peak / 1.0e6:.3g} MW/m^2, "
        f"{heat_power_peak / 1000.0:.3g} kW",
        f"adjacent_radiation_heat_peak: {rad_flux_peak / 1.0e6:.3g} MW/m^2, "
        f"{rad_power_peak / 1000.0:.3g} kW",
        f"diagnostic_failure_mode: {cls.get('diagnostic_failure_mode', 'unknown')}",
        f"termination: {early.get('termination', 'unknown')} "
        f"(steps={early.get('steps', float('nan'))}, "
        f"t_end={early.get('final_time_s', float('nan')):.6g} s, "
        f"history_cap={early.get('history_cap_reached', False)})",
        f"dt_min_median: {early.get('dt_min_s', float('nan')):.3g} s, "
        f"{early.get('dt_median_s', float('nan')):.3g} s",
        f"first_ignition: cell={early.get('first_ignition_cell', -1)}, "
        f"t={early.get('first_ignition_time_s', float('nan')):.6g} s; "
        f"first_erosive={early.get('first_sustained_erosive_time_s', float('nan')):.6g} s",
        f"radiation_activity: enabled={early.get('radiation_enabled', False)}, "
        f"emitters={early.get('max_radiation_emitters', float('nan')):.3g}, "
        f"receivers={early.get('max_radiation_receivers', float('nan')):.3g}, "
        f"first_nonzero={early.get('first_radiation_time_s', float('nan')):.6g} s",
        f"early_extrema_30ms: maxP={early.get('max_pressure_pa', float('nan')) / 1.0e6:.3g} MPa, "
        f"Tgas=[{early.get('min_gas_temperature_k', float('nan')):.3g}, "
        f"{early.get('max_gas_temperature_k', float('nan')):.3g}] K, "
        f"Tsurf_max={early.get('max_surface_temperature_k', float('nan')):.3g} K, "
        f"Mach_max={early.get('max_mach', float('nan')):.3g}",
        f"early_energy_flags: clipping_dominated={early.get('clipping_dominated_energy', False)}, "
        f"radiation_zero_activity={early.get('radiation_enabled_zero_activity', False)}, "
        f"radiation_energy_mismatch={early.get('radiation_energy_mismatch', False)}",
        f"collapse: class={early.get('collapse_class', 'unknown')}, "
        f"first={early.get('first_collapse_time_s', float('nan')):.6g} s, "
        f"branch={early.get('collapse_branch_suspect', 'unknown')}",
    ]
    return "\n".join(lines) + "\n"


def literature_evaluation_report(case_name: str, rows: list[dict[str, Any]]) -> str:
    """Write a compact literature-framed interpretation for probe matrices."""
    by_variant = {str(row.get("variant")): row for row in rows}
    ambient = (
        by_variant.get("ambient_emissivity_0p45")
        or by_variant.get("ambient_emissivity_0.45")
        or by_variant.get("ambient_nominal_T850_dt1e-4")
        or by_variant.get("ambient_initial_gas")
        or {}
    )
    no_rad = (
        by_variant.get("ambient_no_radiation_T850_dt1e-4")
        or by_variant.get("ambient_no_radiation")
        or {}
    )
    baseline = (
        by_variant.get("baseline_hotfill")
        or by_variant.get("baseline_hotfill_T850_dt1e-4")
        or by_variant.get("baseline")
        or {}
    )
    no_erosive = (
        by_variant.get("no_erosive_hotfill")
        or by_variant.get("no_erosive_hotfill_T850_dt1e-4")
        or by_variant.get("no_erosive")
        or {}
    )

    def _fmt_peak(row: dict[str, Any]) -> str:
        if not row:
            return "not run"
        return (
            f"{float(row.get('startup_window_peak_pressure_mpa', float('nan'))):.3g} MPa "
            f"at {float(row.get('startup_window_peak_time_s', float('nan'))):.6g} s; "
            f"mode={row.get('diagnostic_failure_mode', 'unknown')}"
        )

    ambient_mode = str(ambient.get("diagnostic_failure_mode", "unknown"))
    no_rad_mode = str(no_rad.get("diagnostic_failure_mode", "unknown"))
    baseline_mode = str(baseline.get("diagnostic_failure_mode", "unknown"))
    no_erosive_mode = str(no_erosive.get("diagnostic_failure_mode", "unknown"))
    class_counts: dict[str, int] = {}
    for row in rows:
        cls = str(row.get("collapse_class", "unknown"))
        class_counts[cls] = class_counts.get(cls, 0) + 1

    def _is_stable(name: str) -> bool:
        return str(by_variant.get(name, {}).get("collapse_class", "")) == "stable"

    def _is_collapse(name: str) -> bool:
        return str(by_variant.get(name, {}).get("collapse_class", "")) == "collapse"

    sweep_names = [
        "ambient_emissivity_0", "ambient_emissivity_0p05",
        "ambient_emissivity_0p10", "ambient_emissivity_0p20",
        "ambient_emissivity_0p30", "ambient_emissivity_0p40",
        "ambient_emissivity_0p45", "ambient_emissivity_0p50",
        "ambient_emissivity_0p60", "ambient_emissivity_0p75",
        "ambient_emissivity_0p90",
    ]
    stable_sweep = [name for name in sweep_names if _is_stable(name)]
    collapse_sweep = [name for name in sweep_names if _is_collapse(name)]

    if _is_stable("ambient_rad045_receiver_heat_no_sink"):
        dominant_branch = "radiation source/sink accounting"
    elif _is_stable("ambient_rad045_no_erosive"):
        dominant_branch = "erosive feedback after radiation-assisted spread"
    elif any(_is_stable(name) for name in (
        "ambient_rad045_cfl025", "ambient_rad045_cfl010",
        "ambient_rad045_dt2e-5", "ambient_rad045_cells50",
        "ambient_rad045_cells200",
    )):
        dominant_branch = "grid/timestep sensitivity"
    elif str(ambient.get("collapse_class", "")) == "collapse":
        dominant_branch = "PISO/nozzle/front numerical instability"
    else:
        dominant_branch = "not isolated"

    if ambient.get("radiation_enabled_zero_activity", False):
        conclusion = (
            "Nominal ambient+radiation is not yet a physical radiation result: "
            "radiation is enabled but no emitter/receiver activity appears in "
            "the early diagnostic window."
        )
    elif ambient.get("history_cap_reached", False) or ambient.get("clipping_dominated_energy", False):
        conclusion = (
            "Nominal ambient+radiation is dominated by numerical/front behavior "
            "before it can be used as a literature comparison."
        )
    elif ambient_mode == "physically_weak_ignition":
        conclusion = (
            "Nominal ambient+radiation currently behaves like weak pyrogen "
            "pressurization without completing flame spread."
        )
    else:
        conclusion = (
            "Nominal ambient+radiation completed the early transient without "
            "the diagnostic failure flags targeted by this probe."
        )

    if ambient.get("radiation_enabled_zero_activity", False):
        salita_probe_note = (
            "- In this probe, zero early emitter/receiver activity would make "
            "the radiation path a code-path or gating suspect."
        )
    else:
        salita_probe_note = (
            "- In this probe, radiation is active after first ignition; the "
            "suspect behavior is the numerical/front response once radiation "
            "participates."
        )

    lines = [
        f"{case_name} radiation-probe literature evaluation",
        "",
        "Data anchors:",
        f"- ambient nominal: {_fmt_peak(ambient)}",
        f"- ambient no radiation: {_fmt_peak(no_rad)}",
        f"- hot-fill baseline: {_fmt_peak(baseline)}",
        f"- no erosive: {_fmt_peak(no_erosive)}",
        f"- collapse classes: {class_counts}",
        f"- emissivity sweep stable: {stable_sweep}",
        f"- emissivity sweep collapse: {collapse_sweep}",
        f"- dominant branch suspect: {dominant_branch}",
        "",
        "Sutton/DeMar:",
        "- Pyrogen pressurization around the 1-2 MPa range is a plausible chamber-filling/ignition-pressure level.",
        "- Failure to transition from pyrogen heating into flame spread and pressure establishment is not the desired Hasegawa A match.",
        "- Pyrogen heating should be primarily convective; radiation should not dominate before a burning surface exists.",
        "",
        "Salita:",
        "- Radiation should be explicit and physically sourced, not a hidden calibration multiplier.",
        salita_probe_note,
        "",
        "d'Agostino:",
        "- The rapid spike in baseline/no-radiation runs and the weaker no-erosive run support erosive feedback as the correct pressure-establishment driver.",
        "- Do not tune an empirical C_hc-style multiplier until the nominal ambient path is internally consistent.",
        "",
        "Peretz/Pardue/Cavallini:",
        "- Critical-surface-temperature ignition and per-cell spread remain literature-consistent.",
        "- History-cap termination, timestep collapse, clipping dominance, or front ringing should be resolved as numerical pathologies before model calibration.",
        "",
        "Diagnostic conclusion:",
        f"- {conclusion}",
        f"- Comparison modes: ambient={ambient_mode}, no_radiation={no_rad_mode}, "
        f"baseline={baseline_mode}, no_erosive={no_erosive_mode}.",
    ]
    return "\n".join(lines) + "\n"


def write_diagnostic_outputs(diagnostics: dict[str, Any], output_dir,
                             case_name: str) -> None:
    """Write summary and time-series CSVs plus a text report."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[tuple[str, Any]] = []
    _flatten_scalars("", {
        "pressure": diagnostics["pressure"],
        "pyrogen": diagnostics["pyrogen"],
        "sources_at_peak": diagnostics["sources_at_peak"],
        "classification": diagnostics["classification"],
        "early": diagnostics.get("early", {}),
    }, rows)
    spread = diagnostics["ignition_spread"]
    for key in ("first_ignition_time_s", "t10_s", "t50_s", "t90_s",
                "t100_s", "spread_10_90_s", "instant_ignition_collapse",
                "spread_metric_source", "exact_spread_metrics"):
        rows.append((f"ignition_spread.{key}", spread[key]))

    with open(output_dir / f"{case_name}_summary.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        writer.writerows(rows)

    sources = diagnostics["sources"]
    source_keys = [k for k in sources if k != "times_s"]
    with open(output_dir / f"{case_name}_sources.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s"] + source_keys)
        for i, t in enumerate(sources["times_s"]):
            writer.writerow([t] + [sources[k][i] for k in source_keys])

    energy = diagnostics.get("energy", {"times_s": np.array([])})
    energy_keys = [k for k in energy if k != "times_s"]
    with open(output_dir / f"{case_name}_energy_momentum.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s"] + energy_keys)
        for i, t in enumerate(energy["times_s"]):
            writer.writerow([t] + [energy[k][i] for k in energy_keys])

    step_diag = diagnostics.get("step_diagnostics", {"times_s": np.array([])})
    step_keys = [k for k in step_diag if k != "times_s"]
    with open(output_dir / f"{case_name}_step_diagnostics.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s"] + step_keys)
        for i, t in enumerate(step_diag["times_s"]):
            writer.writerow([t] + [step_diag[k][i] for k in step_keys])

    early = diagnostics.get("early", {})
    with open(output_dir / f"{case_name}_early_diagnostics.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for key, value in early.items():
            writer.writerow([key, value])

    collapse_trace = diagnostics.get("collapse_trace", [])
    trace_keys = [
        "event", "event_time_s", "snapshot_index", "snapshot_time_s",
        "history_time_s", "event_cell", "cell", "cell_offset",
        "P_pa", "T_k", "T_surf_k", "Mach", "r_total_m_s",
        "r_erosive_m_s", "radiation_heat_flux_w_m2",
        "mass_source_kg_m_s", "thermal_source_k_kg_m_s",
        "dt_below_1e_8", "pressure_above_100mpa", "mach_above_1e3",
        "clipping_dominated_step", "history_cap_reached",
        "spread_metric_source",
    ]
    with open(output_dir / f"{case_name}_collapse_trace.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=trace_keys)
        writer.writeheader()
        for row in collapse_trace:
            writer.writerow(row)

    first_by_cell = spread["first_burning_time_by_cell_s"]
    with open(output_dir / f"{case_name}_ignition_times.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["cell", "first_burning_time_s"])
        for i, value in enumerate(first_by_cell):
            writer.writerow([i, value])

    with open(output_dir / f"{case_name}_classification.txt", "w") as f:
        f.write(classification_report(diagnostics))


def plot_diagnostic_figures(result: dict, diagnostics: dict[str, Any],
                            output_dir, case_name: str) -> None:
    """Write overview and x-t diagnostic plots."""
    if not HAS_MATPLOTLIB:
        return
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    t = np.asarray(result.get("time", []), dtype=float)
    p = np.asarray(result.get("P_head", []), dtype=float) / 1.0e6
    sources = diagnostics["sources"]
    spread = diagnostics["ignition_spread"]

    fig, axes = plt.subplots(4, 1, figsize=(11, 11), sharex=False)
    axes[0].plot(t, p, "b-", linewidth=1.8)
    axes[0].axvline(
        diagnostics["pressure"]["startup_window_peak_time_s"],
        color="k", linestyle=":", linewidth=1.0, label="startup peak",
    )
    axes[0].axvline(
        diagnostics["pressure"]["global_peak_time_s"],
        color="0.5", linestyle="--", linewidth=0.8, label="global peak",
    )
    axes[0].set_ylabel("P_head [MPa]")
    axes[0].set_title(f"{case_name}: pressure and startup drivers")
    axes[0].legend(loc="best", fontsize=8)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, np.asarray(result.get("mdot_ig", [])) * 1000.0,
                 color="tab:orange", linewidth=1.4, label="mdot_ig")
    axes[1].set_ylabel("mdot_ig [g/s]")
    axes[1].grid(True, alpha=0.3)

    if sources["times_s"].size:
        st = sources["times_s"]
        axes[2].plot(st, sources["normal_fraction"], label="normal")
        axes[2].plot(st, sources["erosive_fraction"], label="erosive")
        axes[2].plot(st, sources["endface_fraction"], label="endface")
        axes[2].plot(st, sources["pyrogen_fraction"], label="pyrogen")
    axes[2].set_ylabel("source fraction")
    axes[2].legend(loc="best", fontsize=8)
    axes[2].grid(True, alpha=0.3)

    if spread["times_s"].size:
        axes[3].plot(spread["times_s"], spread["burning_fraction"],
                     "g-", linewidth=1.6)
    axes[3].set_ylabel("burning fraction")
    axes[3].set_xlabel("time [s]")
    axes[3].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_dir / f"{case_name}_overview.png", dpi=160)
    plt.close(fig)

    snapshots = result.get("snapshots", [])
    if not snapshots:
        return
    snap_t = _snapshot_times(result)
    x_mm = np.asarray(snapshots[0]["x"], dtype=float) * 1000.0
    maps = [
        ("T", "T_gas [K]", stack_snapshots(result, "T")),
        ("T_surf", "T_surf [K]", stack_snapshots(result, "T_surf")),
        ("is_burning", "is_burning [-]", stack_snapshots(result, "is_burning").astype(float)),
        ("r_erosive", "r_erosive [mm/s]", stack_snapshots(result, "r_erosive") * 1000.0),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    for ax, (_key, title, data) in zip(axes.flat, maps):
        pcm = ax.pcolormesh(snap_t, x_mm, data.T, shading="nearest")
        fig.colorbar(pcm, ax=ax)
        ax.set_title(title)
        ax.set_xlabel("time [s]")
        ax.set_ylabel("x [mm]")
    fig.tight_layout()
    fig.savefig(output_dir / f"{case_name}_xt.png", dpi=160)
    plt.close(fig)
