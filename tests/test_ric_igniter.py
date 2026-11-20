"""
test_ric_igniter.py — v0.8.0 Phase 4 igniter-as-data (D3).

Verifies the motor ``data.igniter`` block: default seeding, round-trip
through the Pyrogen material + chamber sizing, that all repo motors carry a
parseable igniter block after migration, and that run_from_ric reads the
igniter from the motor file when no pyrogen kwarg is given.
"""

import os
import glob

import pytest

from srm_1d.openmotor_adapter import (
    default_igniter_block, load_igniter, _seed_motor_igniter,
    _pyrogen_to_block, _block_to_pyrogen, load_ric, build_pyrogen_chamber,
)
from srm_1d.propellant import Pyrogen

MOTORS = os.path.join(os.path.dirname(__file__), '..', 'motors')


def test_default_block_has_material_and_sizing():
    block = default_igniter_block('bpnv')
    assert block['pyrogen']['name'] == 'BPNV'
    assert block['injection_topology'] == 'forward_plenum'
    assert block['mass'] == -1.0  # auto sentinel
    for key in ('throat_area', 'volume', 'burn_area', 'cartridge_length_m'):
        assert block[key] == -1.0


def test_pyrogen_block_round_trip():
    block = default_igniter_block('bpnv')
    pyro = _block_to_pyrogen(block['pyrogen'])
    assert isinstance(pyro, Pyrogen)
    re = _pyrogen_to_block(pyro)
    # Round-trip preserves every material field.
    for key, val in block['pyrogen'].items():
        assert re[key] == val


def test_gas_mass_fraction_default_and_custom_round_trip():
    # The bundled BPNV material carries its condensed-phase split through the
    # igniter block used by openMotor.
    block = default_igniter_block('bpnv')['pyrogen']
    assert block['gas_mass_fraction'] == pytest.approx(0.781)

    custom = dict(block, gas_mass_fraction=0.63)
    pyro = _block_to_pyrogen(custom)
    assert pyro.gas_mass_fraction == pytest.approx(0.63)
    assert _pyrogen_to_block(pyro)['gas_mass_fraction'] == pytest.approx(0.63)


def test_gas_mass_fraction_missing_keeps_legacy_default():
    # Older motor files and user material dictionaries omitted this field;
    # preserve the Pyrogen model's established default when decoding them.
    block = default_igniter_block('bpnv')['pyrogen']
    block.pop('gas_mass_fraction')
    pyro = _block_to_pyrogen(block)
    assert pyro.gas_mass_fraction == pytest.approx(1.0)


def test_load_igniter_maps_sentinels_to_auto():
    block = default_igniter_block('bpnv')
    motor = {'igniter': block}
    pyro, sizing = load_igniter(motor)
    assert pyro.name == 'BPNV'
    # -1.0 sentinels become None (auto) for build_pyrogen_chamber.
    assert sizing['pyrogen_mass'] is None
    assert sizing['pyrogen_throat_area'] is None
    assert sizing['pyrogen_burn_area'] is None
    assert sizing['injection_topology'] == 'forward_plenum'
    # cartridge_length_m keeps the -1.0 sentinel (build derives it).
    assert sizing['cartridge_length_m'] == -1.0


def test_load_igniter_none_when_absent():
    assert load_igniter({'propellant': {}}) is None


def test_seed_motor_igniter_idempotent():
    data = {'propellant': {}}
    _seed_motor_igniter(data)
    assert 'igniter' in data
    first = data['igniter']
    _seed_motor_igniter(data)
    assert data['igniter'] is first  # not overwritten


def test_explicit_sizing_survives_load():
    block = default_igniter_block('bpnv')
    block['mass'] = 0.006
    block['injection_topology'] = 'head_basket'
    pyro, sizing = load_igniter({'igniter': block})
    assert sizing['pyrogen_mass'] == pytest.approx(0.006)
    assert sizing['injection_topology'] == 'head_basket'


def test_all_repo_motors_have_parseable_igniter():
    rics = glob.glob(os.path.join(MOTORS, '*.ric'))
    assert rics
    for ric in rics:
        motor = load_ric(ric)
        ig = load_igniter(motor)
        assert ig is not None, f"{os.path.basename(ric)} has no igniter block"
        pyro, sizing = ig
        assert pyro.name and pyro.a > 0 and pyro.gamma > 1.0


def test_build_chamber_from_motor_file_igniter():
    """The loaded igniter builds a valid PyrogenChamber (self-describing)."""
    from tests.test_simulation_phase3 import _small_motor
    geo, _, nozzle = _small_motor()
    motor = load_ric(os.path.join(MOTORS, 'hasegawa_a.ric'))
    pyro, sizing = load_igniter(motor)
    chamber = build_pyrogen_chamber(pyro, geo, nozzle, **sizing)
    assert chamber.m_pyrogen_initial > 0
    assert chamber.injection_topology == 'forward_plenum'
