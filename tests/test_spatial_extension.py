import pytest

from scripts.run_spatial_extension import (
    _approximate_observed_order,
    _successive_changes,
)


def test_successive_changes_use_fine_value_as_denominator():
    rows = [
        {"target_propellant_cells": 100, "metric": 9.0},
        {"target_propellant_cells": 200, "metric": 10.0},
    ]

    changes = _successive_changes(rows, "metric")

    assert changes[0]["relative_change_percent"] == pytest.approx(10.0)


def test_observed_order_and_extrapolation_for_monotone_factor_two_series():
    analysis = _approximate_observed_order([8.0, 9.0, 9.5, 9.75])

    assert analysis["triplets"][0]["approximate_order"] == pytest.approx(1.0)
    assert analysis["triplets"][0]["approximate_extrapolated_value"] == 10.0
    assert analysis["triplets"][1]["approximate_order"] == pytest.approx(1.0)


def test_observed_order_is_withheld_for_oscillatory_series():
    analysis = _approximate_observed_order([8.0, 9.0, 8.5])

    assert analysis["triplets"][0]["monotone"] is False
    assert analysis["triplets"][0]["approximate_order"] is None
