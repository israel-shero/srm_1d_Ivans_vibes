import numpy as np
import pytest

from scripts.run_numerical_verification import (
    FULL_CFL_LEVELS,
    STARTUP_CFL_LEVELS,
    _adjacent_changes,
    _extract_metrics,
)


def test_extract_metrics_counts_clipping_and_fill_window():
    result = {
        "time": np.array([0.0, 0.002, 0.004]),
        "P_head": np.array([1e5, 3e6, 2e6]),
        "max_mach": np.array([0.0, 4.0, 9.0]),
        "min_pressure": np.array([1e5, 9e4, 8e4]),
        "max_gas_temperature": np.array([300.0, 1000.0, 900.0]),
        "clipping_correction_power": np.array([0.0, -2.0, 0.0]),
        "energy_residual": np.array([0.0, 1.0, -1.0]),
        "thermal_source_power": np.array([0.0, 10.0, 10.0]),
        "summary": {
            "cells": 12, "steps": 3, "termination": "t_max reached",
            "mass_balance_error": 0.01, "dt_min": 1e-6,
            "dt_median": 2e-6, "wall_time": 0.1,
        },
    }

    metrics = _extract_metrics(result)
    assert metrics["peak_pressure_mpa"] == pytest.approx(3.0)
    assert metrics["peak_time_s"] == pytest.approx(0.002)
    assert metrics["max_fill_mach"] == pytest.approx(4.0)
    assert metrics["clipping_active_steps"] == 1
    assert metrics["clipping_active_fraction"] == pytest.approx(1 / 3)
    assert metrics["integrated_abs_clipping_energy_j"] > 0.0


def test_extract_metrics_includes_full_burn_performance():
    result = {
        "time": np.array([0.0, 0.002]),
        "P_head": np.array([1e5, 2e6]),
        "max_mach": np.array([0.0, 0.2]),
        "min_pressure": np.array([1e5, 9e4]),
        "max_gas_temperature": np.array([300.0, 1000.0]),
        "clipping_correction_power": np.zeros(2),
        "energy_residual": np.zeros(2),
        "thermal_source_power": np.ones(2),
        "summary": {
            "cells": 12, "steps": 2, "termination": "t_max reached",
            "mass_balance_error": 0.01, "dt_min": 1e-6,
            "dt_median": 2e-6, "wall_time": 0.1,
            "t_first_burnout": 1.2, "D_throat_final": 0.02,
        },
    }
    performance = {
        "total_impulse": 100.0, "average_thrust": 50.0,
        "peak_thrust": 80.0, "burn_time": 2.0, "average_Isp": 200.0,
    }

    metrics = _extract_metrics(result, performance)
    assert metrics["total_impulse_ns"] == pytest.approx(100.0)
    assert metrics["performance_burn_time_s"] == pytest.approx(2.0)
    assert metrics["first_burnout_time_s"] == pytest.approx(1.2)


def test_adjacent_changes_use_finer_value_as_denominator():
    rows = [
        {"label": "coarse", "peak_pressure_mpa": 9.8},
        {"label": "fine", "peak_pressure_mpa": 10.0},
    ]
    change = _adjacent_changes(rows, "peak_pressure_mpa")[0]
    assert change["relative_change_percent"] == pytest.approx(2.0)


def test_full_profile_avoids_known_history_capacity_limit():
    assert STARTUP_CFL_LEVELS[-1] == pytest.approx(0.075)
    assert FULL_CFL_LEVELS == (0.3, 0.2, 0.1)
