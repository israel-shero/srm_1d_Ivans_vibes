"""Tests for openMotor propellant-catalog imports."""

import pytest
import yaml

from srm_1d.propellant_catalog import (
    get_propellant, list_propellants, load_propellant_catalog,
    missing_transport_fields, with_catalog_propellant,
)


def test_library_has_77_distinct_entries():
    names = list_propellants()
    assert len(names) == 77
    assert len(names) == len(set(names))
    assert "RCS - Super Thunder (H13ST Version)" in names


def test_selected_propellant_is_a_detached_copy():
    first = get_propellant("RCS - Super Thunder (H13ST Version)")
    first["tabs"][0]["a"] = -999
    assert get_propellant("RCS - Super Thunder (H13ST Version)")["tabs"][0]["a"] != -999


def test_missing_transport_is_reported_without_fabrication():
    item = get_propellant("RCS - Super Thunder (H13ST Version)")
    assert missing_transport_fields(item) == {
        0: ["mu", "kThermalFrozen", "cpFrozen"]
    }
    assert item["tabs"][0]["mu"] == 0.0


def test_replace_only_flat_motor_propellant():
    original = {
        "propellant": {"name": "before"},
        "grains": [1],
        "nozzle": {"throat": 1},
    }
    updated = with_catalog_propellant(
        original, "RCS - Super Thunder (H13ST Version)"
    )
    assert updated["propellant"]["name"] == "RCS - Super Thunder (H13ST Version)"
    assert updated["grains"] == original["grains"]
    assert original["propellant"] == {"name": "before"}


def test_reject_duplicate_names(tmp_path):
    path = tmp_path / "duplicate.yaml"
    entry = {
        "name": "duplicate",
        "density": 1500,
        "tabs": [{
            "a": 0.1, "n": 0.3, "k": 1.2, "m": 25, "t": 3000,
            "minPressure": 0, "maxPressure": 1e7,
        }],
    }
    path.write_text(
        yaml.safe_dump({"data": [entry, entry]}), encoding="utf-8"
    )
    with pytest.raises(ValueError, match="Duplicate"):
        load_propellant_catalog(path)


def test_exact_name_failure_is_informative():
    with pytest.raises(KeyError, match="Closest names"):
        get_propellant("RCS - Super Thundr (H13ST Version)")
