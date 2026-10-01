import pytest

from scripts.run_fill_configuration_comparison import (
    HISTORICAL_RECORDED,
    LEGACY_BPNV_MARKER,
    _relative_change,
    _resolve_options,
    _variants,
)


def test_variants_separate_pyrogen_from_historical_nonpyrogen_knobs():
    current = {
        "pyrogen": "mtv",
        "roughness": 35.0e-6,
        "kappa": 0.45,
        "T_ignition": 900.0,
        "P_cutoff": 50_000.0,
        "t_max": 3.0,
    }

    variants = _variants(current)

    assert variants["current_mtv"]["pyrogen"] == "mtv"
    assert variants["current_bpnv"]["pyrogen"] == "bpnv"
    assert variants["current_bpnv"]["T_ignition"] == 900.0
    assert variants["current_legacy_bpnv"]["pyrogen"] == LEGACY_BPNV_MARKER
    assert variants["historical_mtv"]["pyrogen"] == "mtv"
    assert variants["historical_mtv"]["T_ignition"] == 756.0
    assert variants["historical_bpnv"]["roughness"] == pytest.approx(32.0e-6)
    assert all(item["t_max"] == 0.05 for item in variants.values())


def test_legacy_bpnv_marker_resolves_pre_fix_gas_properties():
    options = _resolve_options({"pyrogen": LEGACY_BPNV_MARKER, "kappa": 0.44})

    assert options["kappa"] == 0.44
    assert options["pyrogen"].T_flame == 2800.0
    assert options["pyrogen"].M == 0.030
    assert options["pyrogen"].gamma == 1.25
    assert options["pyrogen"].gas_mass_fraction == 1.0


def test_relative_change_uses_named_reference():
    assert _relative_change(4.0, 5.0) == pytest.approx(25.0)
    assert _relative_change(0.0, 1.0) is None
    assert HISTORICAL_RECORDED[200]["max_fill_mach"] == 12.04
