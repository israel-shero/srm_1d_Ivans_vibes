import numpy as np
import pytest

from cases import (
    CaseNotFoundError,
    MeasurementUnavailableError,
    get_case,
    load_measurement,
    load_registry,
)


def test_registry_has_unique_expected_cases():
    registry = load_registry()
    case_ids = [case["id"] for case in registry["cases"]]
    assert len(case_ids) == len(set(case_ids)) == 8
    assert {"hasegawa_a", "zerox", "chunc", "ballsstick_3in"} <= set(case_ids)


def test_get_case_returns_copy():
    first = get_case("chunc")
    first["label"] = "changed"
    assert get_case("chunc")["label"] == "Chunc / machbusterNew"


def test_unknown_case_raises():
    with pytest.raises(CaseNotFoundError):
        get_case("not-a-case")


def test_hasegawa_symbol_trace_loads_without_mutating_source():
    first = load_measurement("hasegawa_a")
    second = load_measurement("hasegawa_a")
    assert len(first["time"]) == len(first["pressure"]) == 36
    first["time"][0] = 99.0
    assert second["time"][0] == 0.0


def test_zerox_reversed_columns_load_correctly():
    trace = load_measurement("zerox")
    assert len(trace["time"]) == 58
    np.testing.assert_allclose(trace["time"][:2], [0.0, 0.155])
    np.testing.assert_allclose(trace["pressure"][:2], [0.0, 0.0239])


def test_chunc_csv_loads_in_mpa():
    trace = load_measurement("chunc")
    assert len(trace["time"]) == 346
    np.testing.assert_allclose(trace["time"][0], 0.01)
    np.testing.assert_allclose(trace["pressure"][0], 8.466000756)


def test_pathfinder_pa_converts_to_mpa():
    trace = load_measurement("pathfinder54_overlay")
    assert len(trace["time"]) == 409
    np.testing.assert_allclose(trace["pressure"][0], 0.0901485969)


@pytest.mark.parametrize("case_id", ["ballsstick_3in", "st2_54_2800", "hlbl_test"])
def test_unsafe_or_absent_measurements_are_refused(case_id):
    with pytest.raises(MeasurementUnavailableError):
        load_measurement(case_id)
