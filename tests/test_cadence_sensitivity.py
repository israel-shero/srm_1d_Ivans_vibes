import copy

import pytest

from scripts import run_cadence_sensitivity as cadence


def test_run_matrix_is_the_explicit_fixed_grid_cadence_matrix():
    assert cadence.TARGET_CELLS == 800
    assert cadence.CFL_TARGET == 0.075
    assert cadence.T_MAX_S == 0.05
    assert cadence.CADENCE_MATRIX == (
        (161, 161), (80, 161), (161, 80), (80, 80),
    )
    assert cadence._build_run_matrix() == [
        {
            "target_propellant_cells": 800,
            "burn_update_interval": burn,
            "geometry_update_interval": geometry,
        }
        for burn, geometry in cadence.CADENCE_MATRIX
    ]


def test_comparisons_calculate_changes_against_baseline_and_mark_missing_rows():
    baseline_metrics = {metric: 10.0 for metric in cadence.COMPARISON_METRICS}
    baseline = {
        **cadence._build_run_matrix()[0], "status": "complete",
        "metrics": baseline_metrics,
    }
    burn_only = {
        **cadence._build_run_matrix()[1], "status": "complete",
        "metrics": {**baseline_metrics, "peak_pressure_mpa": 10.2},
    }
    geometry_only = {
        **cadence._build_run_matrix()[2], "status": "failed", "metrics": None,
    }
    both = {
        **cadence._build_run_matrix()[3], "status": "complete",
        "metrics": {**baseline_metrics, "peak_time_s": 9.0},
    }

    comparisons = cadence._build_comparisons(
        [baseline, burn_only, geometry_only, both]
    )
    assert [row["status"] for row in comparisons] == [
        "compared", "unavailable", "compared",
    ]
    pressure = comparisons[0]["metrics"]["peak_pressure_mpa"]
    assert pressure["absolute_change"] == pytest.approx(0.2)
    assert pressure["relative_change_percent"] == pytest.approx(2.0)
    assert pressure["exceeds_screen"] is False
    missing = comparisons[1]["metrics"]["minimum_pressure_pa"]
    assert missing["value"] is None
    assert missing["reference_value"] == 10.0
    assert missing["relative_change_percent"] is None
    assert missing["exceeds_screen"] is None
    assert comparisons[2]["metrics"]["peak_time_s"]["relative_change_percent"] == 10.0
    assert comparisons[2]["metrics"]["peak_time_s"]["exceeds_screen"] is True


def test_comparison_relative_change_handles_none_and_zero_reference():
    assert cadence._relative_change_percent(None, 2.0) is None
    assert cadence._relative_change_percent(1.0, 0.0) is None
    assert cadence._relative_change_percent(12.0, 10.0) == pytest.approx(20.0)


def test_summary_and_claims_preserve_screening_limits():
    rows = [
        {**config, "status": "complete"}
        for config in cadence._build_run_matrix()
    ]
    assert cadence._summarize_run_statuses(rows) == {
        "all_runs_complete": True,
        "counts": {"complete": 4},
        "study_status": "complete",
    }
    rows[-1]["status"] = "failed"
    summary = cadence._summarize_run_statuses(rows)
    assert summary["all_runs_complete"] is False
    assert summary["study_status"] == "complete_with_failures"
    payload = cadence._new_payload(
        cadence.PROJECT_ROOT / "motors/machbusterNew.ric", {}
    )
    assert payload["claims"] == {
        "convergence_established": False,
        "physics_validated": False,
        "research_note_results_reproduced": False,
    }


def test_resume_validation_checks_inputs_configuration_and_matrix():
    motor_path = cadence.PROJECT_ROOT / "motors/machbusterNew.ric"
    payload = {
        "provenance": {
            "motor_sha256": cadence._sha256(motor_path),
            "baseline_config_sha256": cadence._sha256(cadence.CONFIG_PATH),
            "motor_path": "motors/machbusterNew.ric",
        },
        "configuration": {"t_max": 0.05},
        "scope": {"cadence_matrix": [list(pair) for pair in cadence.CADENCE_MATRIX]},
    }
    cadence._validate_resume(payload, motor_path, {"t_max": 0.05})

    bad_config = copy.deepcopy(payload)
    bad_config["configuration"] = {"t_max": 0.06}
    with pytest.raises(ValueError, match="configuration"):
        cadence._validate_resume(bad_config, motor_path, {"t_max": 0.05})

    bad_matrix = copy.deepcopy(payload)
    bad_matrix["scope"]["cadence_matrix"] = [[1, 1]]
    with pytest.raises(ValueError, match="cadence matrix"):
        cadence._validate_resume(bad_matrix, motor_path, {"t_max": 0.05})

    bad_checksum = copy.deepcopy(payload)
    bad_checksum["provenance"]["baseline_config_sha256"] = "invalid"
    with pytest.raises(ValueError, match="baseline-config checksum"):
        cadence._validate_resume(bad_checksum, motor_path, {"t_max": 0.05})
