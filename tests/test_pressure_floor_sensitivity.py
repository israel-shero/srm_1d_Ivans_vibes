import pytest

from scripts.run_pressure_floor_sensitivity import (
    PRESSURE_FLOORS_PA,
    SENSITIVITY_METRICS,
    SENSITIVITY_SCREEN_PERCENT,
    TARGET_CELLS,
    _build_comparisons,
    _build_run_matrix,
    _relative_change_percent,
    _summarize_run_statuses,
)


def test_run_matrix_covers_each_requested_grid_and_floor():
    matrix = _build_run_matrix()

    assert len(matrix) == 6
    assert [row["target_propellant_cells"] for row in matrix] == [
        800, 800, 800, 1600, 1600, 1600,
    ]
    assert {row["pressure_floor_pa"] for row in matrix} == set(PRESSURE_FLOORS_PA)
    assert TARGET_CELLS == (800, 1600)


def test_relative_change_uses_reference_magnitude_and_handles_missing_values():
    assert _relative_change_percent(10.2, 10.0) == pytest.approx(2.0)
    assert _relative_change_percent(9.8, 10.0) == pytest.approx(2.0)
    assert _relative_change_percent(1.0, 0.0) is None
    assert _relative_change_percent(None, 10.0) is None


def test_comparisons_screen_metrics_against_same_grid_reference():
    metrics = {metric: 10.0 for metric in SENSITIVITY_METRICS}
    rows = [
        {
            "target_propellant_cells": 800,
            "pressure_floor_pa": 1000.0,
            "status": "complete",
            "metrics": metrics,
        },
        {
            "target_propellant_cells": 800,
            "pressure_floor_pa": 100.0,
            "status": "complete",
            "metrics": {
                "peak_pressure_mpa": 10.1,
                "peak_time_s": 10.2,
                "max_fill_mach": 9.9,
            },
        },
        {
            "target_propellant_cells": 1600,
            "pressure_floor_pa": 1000.0,
            "status": "complete",
            "metrics": metrics,
        },
        {
            "target_propellant_cells": 1600,
            "pressure_floor_pa": 10000.0,
            "status": "failed",
            "metrics": None,
        },
    ]

    comparisons = _build_comparisons(rows)
    comparison = next(
        row for row in comparisons
        if row["target_propellant_cells"] == 800
        and row["pressure_floor_pa"] == 100.0
    )
    assert comparison["status"] == "compared"
    assert comparison["metrics"]["peak_pressure_mpa"]["exceeds_screen"] is False
    assert comparison["metrics"]["peak_time_s"]["relative_change_percent"] == pytest.approx(2.0)
    assert comparison["metrics"]["peak_time_s"]["exceeds_screen"] is False
    assert comparison["metrics"]["max_fill_mach"]["relative_change_percent"] == pytest.approx(1.0)
    failed = next(
        row for row in comparisons
        if row["target_propellant_cells"] == 1600
        and row["pressure_floor_pa"] == 10000.0
    )
    assert failed["status"] == "unavailable"
    assert all(
        metric["relative_change_percent"] is None
        for metric in failed["metrics"].values()
    )
    assert SENSITIVITY_SCREEN_PERCENT == 2.0


def test_study_status_surfaces_unhealthy_or_failed_rows():
    complete = [
        {**config, "status": "complete"} for config in _build_run_matrix()
    ]
    assert _summarize_run_statuses(complete) == {
        "all_runs_complete": True,
        "counts": {"complete": 6},
        "study_status": "complete",
    }

    complete[-1]["status"] = "unhealthy"
    summary = _summarize_run_statuses(complete)
    assert summary["all_runs_complete"] is False
    assert summary["study_status"] == "complete_with_failures"
    assert summary["counts"] == {"complete": 5, "unhealthy": 1}
