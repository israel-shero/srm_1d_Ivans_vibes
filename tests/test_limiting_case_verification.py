import json

import numpy as np
import pytest

from scripts.run_limiting_case_verification import (
    _closure_limiting_cases,
    _relative_change,
    _steady_interval_metrics,
    run_limiting_cases,
)
from srm_1d.openmotor_adapter import convert_propellant, load_ric


def test_steady_interval_metrics_identifies_stationary_pressure():
    result = {
        "time": np.linspace(0.0, 2.0, 21),
        "P_head": np.full(21, 8.0e6),
    }

    metrics = _steady_interval_metrics(result, 0.6, 1.8)

    assert metrics["status"] == "evaluated"
    assert metrics["mean_pressure_mpa"] == pytest.approx(8.0)
    assert metrics["pressure_cv"] == 0.0
    assert metrics["pressure_slope_percent_per_s"] == pytest.approx(0.0)
    assert metrics["screen"]["passed"]


def test_steady_interval_metrics_rejects_large_pressure_slope():
    time = np.linspace(0.0, 2.0, 21)
    result = {
        "time": time,
        "P_head": 8.0e6 + 1.0e6 * time,
    }

    metrics = _steady_interval_metrics(result, 0.6, 1.8)

    assert metrics["pressure_slope_percent_per_s"] > 5.0
    assert not metrics["screen"]["passed"]


def test_steady_interval_requires_two_samples():
    result = {
        "time": np.array([0.0, 1.0]),
        "P_head": np.array([1.0e5, 1.0e5]),
    }

    metrics = _steady_interval_metrics(result, 0.2, 0.8)

    assert metrics["status"] == "not_evaluated"


def test_relative_change_uses_reference_magnitude():
    assert _relative_change(10.0, 8.0) == pytest.approx(-20.0)
    assert _relative_change(0.0, 1.0) is None


def test_closure_limiting_cases_use_converted_propellant():
    motor = load_ric("motors/machbusterNew.ric")
    propellant = convert_propellant(motor["propellant"])

    checks = _closure_limiting_cases(propellant)

    assert checks["zero_crossflow"]["normal_rate_preserved"]
    assert checks["finite_crossflow"]["normal_rate_preserved"]
    assert checks["finite_crossflow"]["erosive_rate_m_s"] > 0.0


def test_resume_rejects_changed_motor_checksum(tmp_path):
    manifest = tmp_path / "limiting_cases.json"
    manifest.write_text(json.dumps({
        "provenance": {
            "motor_sha256": "wrong",
            "baseline_config_sha256": "wrong",
        },
        "runs": {},
    }))

    with pytest.raises(ValueError, match="motor checksum"):
        run_limiting_cases(resume_path=manifest)
