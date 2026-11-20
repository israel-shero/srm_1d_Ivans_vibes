import csv

import numpy as np
import pytest

from scripts.run_numerical_verification import (
    FULL_CFL_LEVELS,
    STARTUP_CFL_LEVELS,
    _adjacent_changes,
    _extract_metrics,
    _ignition_refresh_delays,
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
                "threshold_pa": 1000.0,
                "activation_count_by_cell": np.array([0, 0]),
                "duration_s_by_cell": np.array([0.0, 0.0]),
                "first_activation_time_s_by_cell": np.array([np.nan, np.nan]),
                "last_activation_time_s_by_cell": np.array([np.nan, np.nan]),
                "maximum_raw_deficit_pa_by_cell": np.array([0.0, 2500.0]),
                "maximum_deficit_event": {
                    "deficit_pa": np.array([0.0, 2500.0]),
                    "step_start_time_s": np.array([0.0, 0.012]),
                    "nozzle_state_before_first": np.array([0, 2]),
                },
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


def test_ignition_refresh_delays_use_next_strict_cadence_boundary():
    result = {
        "time": np.array([0.0, 0.1, 0.2, 0.3, 0.4]),
        "ignition_time_by_cell": np.array([0.0, 0.1, 0.2, 1.0e10]),
    }

    delays = _ignition_refresh_delays(result, interval=2)

    np.testing.assert_allclose(delays, [0.2, 0.1, 0.2])


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


def test_extract_metrics_includes_history_storage_provenance():
    result = _localized_result()
    result["summary"].update({
        "history_capacity": 1234,
        "history_storage": "memmap",
        "history_storage_bytes": 414624,
    })

    metrics = _extract_metrics(result)

    assert metrics["history_capacity_rows"] == 1234
    assert metrics["history_storage"] == "memmap"
    assert metrics["history_storage_bytes"] == 414624


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
    assert metrics["limit_pressure_floor_threshold_pa"] == 1000.0
    assert metrics["limit_pressure_floor_maximum_raw_deficit_pa"] == 2500.0
    assert metrics["limit_pressure_floor_maximum_raw_deficit_x_m"] == 0.2
    assert metrics[
        "limit_pressure_floor_maximum_deficit_event_deficit_pa"
    ] == 2500.0
    assert metrics[
        "limit_pressure_floor_maximum_deficit_event_step_start_time_s"
    ] == 0.012
    assert metrics[
        "limit_pressure_floor_maximum_deficit_event_nozzle_state_before_first"
    ] == 2


def test_extract_metrics_omits_pressure_deficit_location_when_inactive():
    result = _localized_result()
    result["numerical_limits"]["pressure_floor"][
        "maximum_raw_deficit_pa_by_cell"
    ] = np.zeros(2)

    metrics = _extract_metrics(result)

    assert metrics["limit_pressure_floor_maximum_raw_deficit_pa"] == 0.0
    assert metrics["limit_pressure_floor_maximum_raw_deficit_x_m"] is None
    assert metrics[
        "limit_pressure_floor_maximum_deficit_event_deficit_pa"
    ] is None
    assert metrics[
        "limit_pressure_floor_maximum_deficit_event_step_start_time_s"
    ] is None
    assert metrics[
        "limit_pressure_floor_maximum_deficit_event_nozzle_state_before_first"
    ] is None


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
    assert rows[1]["pressure_floor_maximum_raw_deficit_pa"] == "2500.0"
    assert rows[1]["pressure_floor_maximum_deficit_event_deficit_pa"] == "2500.0"
    assert rows[1]["pressure_floor_maximum_deficit_event_step_start_time_s"] == "0.012"
    assert rows[1]["pressure_floor_maximum_deficit_event_nozzle_state_before_first"] == "2"


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


@pytest.mark.parametrize('termination', [
    'history array full', 'numerical collapse aborted', 'canceled by user',
    'pressure cutoff', 'complete burnout',
])
def test_full_duration_point_rejects_premature_termination(termination):
    from scripts.run_numerical_verification import _require_requested_duration

    with pytest.raises(RuntimeError, match='did not reach configured t_max'):
        _require_requested_duration({'summary': {'termination': termination}})


def test_full_duration_point_accepts_requested_stop():
    from scripts.run_numerical_verification import _require_requested_duration

    _require_requested_duration({'summary': {'termination': 't_max reached'}})


def test_generated_geometry_fingerprint_detects_changed_fmm_table():
    from types import SimpleNamespace
    from scripts.run_numerical_verification import _geometry_fingerprints

    tables = {'N': 2, 'fmm_port_flat': np.array([0.01, 0.02])}
    geometry = SimpleNamespace(compile_geometry_arrays=lambda: tables)
    before = _geometry_fingerprints(geometry)
    tables['fmm_port_flat'][0] = 0.015
    after = _geometry_fingerprints(geometry)

    assert before['N'] == after['N'] == 2
    assert before['fmm_port_flat']['sha256'] != after['fmm_port_flat']['sha256']


@pytest.mark.parametrize('termination, expected_status', [
    ('t_max reached', 'complete'), ('history array full', 'failed'),
])
def test_single_point_records_selected_closure_and_rejects_truncation(
    tmp_path, monkeypatch, termination, expected_status,
):
    import json
    from types import SimpleNamespace
    import scripts.run_numerical_verification as runner

    result = _localized_result()
    result['summary']['termination'] = termination
    seen_options = {}

    def fake_run(path, **options):
        seen_options.update(options)
        return result, {}, None, SimpleNamespace(
            compile_geometry_arrays=lambda: {'N': 2},
        ), None

    monkeypatch.setattr(runner, 'artifact_dir', lambda *a, **kw: tmp_path)
    monkeypatch.setattr(runner, '_source_fingerprints', lambda: {'solver': 'abc'})
    monkeypatch.setattr(runner, '_runtime_versions', lambda: {'python': 'test'})
    monkeypatch.setattr(runner, 'run_from_ric', fake_run)
    monkeypatch.setattr(runner, 'verify_run_health', lambda *a, **kw: True)
    call = lambda: runner.run_single_point(
        'full', 100, 0.075, erosive_closure='mukunda_paul',
    )
    if expected_status == 'failed':
        with pytest.raises(RuntimeError, match='history array full'):
            call()
    else:
        call()

    payload = json.loads((tmp_path / 'point.json').read_text())
    assert seen_options['erosive_closure'] == 'mukunda_paul'
    assert payload['configuration']['cfl_target'] == 0.075
    assert payload['status'] == expected_status
    assert payload['metrics']['termination'] == termination
    assert payload['metrics']['recorded_end_time_s'] == 0.002
    assert payload['provenance']['source_sha256'] == {'solver': 'abc'}
