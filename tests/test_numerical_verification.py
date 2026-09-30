import csv

import numpy as np
import pytest

from scripts.run_numerical_verification import (
    FULL_CFL_LEVELS,
    STARTUP_CFL_LEVELS,
    _adjacent_changes,
    _extract_metrics,
    _write_limit_csv,
    _write_mach_limit_csv,
)


def _localized_result():
    return {
        "time": np.array([0.0, 0.002]),
        "P_head": np.array([1e5, 2e6]),
        "max_mach": np.array([0.0, 0.2]),
        "min_pressure": np.array([1e5, 9e4]),
        "max_gas_temperature": np.array([300.0, 1000.0]),
        "clipping_correction_power": np.zeros(2),
        "energy_residual": np.zeros(2),
        "thermal_source_power": np.ones(2),
        "summary": {
            "cells": 2, "steps": 2, "termination": "t_max reached",
            "mass_balance_error": 0.0, "dt_min": 1e-6,
            "dt_median": 2e-6, "wall_time": 0.1,
        },
        "numerical_limits": {
            "x_m": np.array([0.1, 0.2]),
            "temperature_floor": {
                "activation_count_by_cell": np.array([0, 0]),
                "duration_s_by_cell": np.array([0.0, 0.0]),
                "first_activation_time_s_by_cell": np.array([np.nan, np.nan]),
                "last_activation_time_s_by_cell": np.array([np.nan, np.nan]),
                "absolute_correction_energy_j_by_cell": np.array([0.0, 0.0]),
            },
            "temperature_ceiling": {
                "activation_count_by_cell": np.array([2, 1]),
                "duration_s_by_cell": np.array([0.002, 0.001]),
                "first_activation_time_s_by_cell": np.array([0.01, 0.02]),
                "last_activation_time_s_by_cell": np.array([0.04, 0.03]),
                "absolute_correction_energy_j_by_cell": np.array([3.0, 5.0]),
            },
            "pressure_floor": {
                "activation_count_by_cell": np.array([0, 0]),
                "duration_s_by_cell": np.array([0.0, 0.0]),
                "first_activation_time_s_by_cell": np.array([np.nan, np.nan]),
                "last_activation_time_s_by_cell": np.array([np.nan, np.nan]),
            },
            "port_mach_cap": {
                "threshold_mach": 1.0,
                "x_m": np.array([0.15]),
                "activation_count_by_face": np.array([3]),
                "duration_s_by_face": np.array([0.003]),
                "first_activation_time_s_by_face": np.array([0.01]),
                "last_activation_time_s_by_face": np.array([0.04]),
            },
        },
    }


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


def test_extract_metrics_localizes_numerical_limits():
    result = _localized_result()

    metrics = _extract_metrics(result)
    assert metrics["limit_temperature_ceiling_active_cells"] == 2
    assert metrics["limit_temperature_ceiling_total_cell_duration_s"] == pytest.approx(0.003)
    assert metrics["limit_temperature_ceiling_max_duration_x_m"] == pytest.approx(0.1)
    assert metrics["limit_temperature_ceiling_max_energy_x_m"] == pytest.approx(0.2)
    assert metrics["limit_temperature_ceiling_first_activation_time_s"] == 0.01
    assert metrics["limit_temperature_ceiling_last_activation_time_s"] == 0.04
    assert metrics["limit_port_mach_cap_active_faces"] == 1
    assert metrics["limit_port_mach_cap_total_activations"] == 3


def test_limit_csv_preserves_each_cell(tmp_path):
    result = _localized_result()
    path = tmp_path / "limits.csv"
    _write_limit_csv(path, result)

    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 2
    assert rows[1]["x_m"] == "0.2"
    assert rows[1]["temperature_ceiling_activations"] == "1"
    assert rows[1]["temperature_ceiling_absolute_correction_energy_j"] == "5.0"


def test_mach_limit_csv_preserves_each_interior_face(tmp_path):
    result = _localized_result()
    path = tmp_path / "mach_limits.csv"
    _write_mach_limit_csv(path, result)

    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 1
    assert rows[0]["interior_face_index"] == "1"
    assert rows[0]["x_m"] == "0.15"
    assert rows[0]["activations"] == "3"


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
