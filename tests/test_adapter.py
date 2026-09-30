"""Tests for the openMotor adapter."""
import pytest
import numpy as np

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

from srm_1d.openmotor_adapter import (
    convert_propellant, convert_geometry, convert_nozzle,
    ric_to_sim_args, load_pyrogen, run_from_ric,
    _INHIBIT_MAP,
)


# Sample .ric data (from the BATES test motor)
SAMPLE_RIC_PROPELLANT = {
    'density': 1700.0,
    'name': 'Hasegawa A Prop',
    'tabs': [{
        'a': 4.821e-05,
        'k': 1.19,
        'm': 25.4,
        'maxPressure': 10000000.0,
        'minPressure': 0.0,
        'n': 0.3,
        't': 3041.0,
    }],
}

SAMPLE_RIC_GRAINS = [
    {'type': 'BATES', 'properties': {
        'coreDiameter': 0.038, 'diameter': 0.07,
        'inhibitedEnds': 'Neither', 'length': 0.12}},
    {'type': 'BATES', 'properties': {
        'coreDiameter': 0.038, 'diameter': 0.07,
        'inhibitedEnds': 'Neither', 'length': 0.12}},
    {'type': 'BATES', 'properties': {
        'coreDiameter': 0.038, 'diameter': 0.07,
        'inhibitedEnds': 'Neither', 'length': 0.12}},
    {'type': 'BATES', 'properties': {
        'coreDiameter': 0.038, 'diameter': 0.07,
        'inhibitedEnds': 'Neither', 'length': 0.12}},
]

SAMPLE_RIC_NOZZLE = {
    'throat': 0.02, 'exit': 0.04, 'divAngle': 15.0,
    'efficiency': 0.95, 'erosionCoeff': 1.5e-10,
    'slagCoeff': 0.0, 'convAngle': 60.0, 'throatLength': 0.0,
}

SAMPLE_RIC_CONFIG = {
    'ambPressure': 101356.5, 'timestep': 0.025,
    'burnoutWebThres': 2.54e-05, 'burnoutThrustThres': 0.1,
    'maxPressure': 10342500.0,
}

SAMPLE_GAS_PROPS = {'mu': 8.842e-5, 'k': 0.3685, 'Cp': 2060.0}


def _write_minimal_ric(path):
    """Write a small safe_load-compatible .ric for adapter tests."""
    if not HAS_YAML:
        pytest.skip("PyYAML not installed")
    data = {
        'version': [0, 0, 0],
        'data': {
            'propellant': SAMPLE_RIC_PROPELLANT,
            'grains': [SAMPLE_RIC_GRAINS[0]],
            'nozzle': SAMPLE_RIC_NOZZLE,
            'config': SAMPLE_RIC_CONFIG,
        },
    }
    path.write_text(yaml.safe_dump(data), encoding='utf-8')


class TestPyrogenLoading:
    def test_load_builtin_pyrogen(self):
        pyro = load_pyrogen('bpnv')
        assert pyro.name == 'BPNV'
        assert pyro.a > 0.0
        assert pyro.impetus_W == pytest.approx(5000.0)
        assert pyro.heat_flux_cal_cm2_s == pytest.approx(69.4)
        mtv = load_pyrogen('mtv')
        assert mtv.name == 'MTV'
        assert mtv.heat_flux_cal_cm2_s == pytest.approx(110.0)

    def test_missing_pyrogen_raises_before_run(self):
        with pytest.raises(ValueError, match="Unknown pyrogen"):
            load_pyrogen('not-a-real-pyrogen')

    def test_run_from_ric_discovers_sibling_pyrogen(self, tmp_path, monkeypatch):
        ric_path = tmp_path / 'motor.ric'
        _write_minimal_ric(ric_path)
        (tmp_path / 'motor.pyrogen.yaml').write_text(
            "name: sibling\n"
            "a: 2.0e-5\n"
            "n: 0.5\n"
            "rho: 1700.0\n"
            "T_flame: 2800.0\n"
            "M: 0.030\n"
            "gamma: 1.25\n",
            encoding='utf-8',
        )

        captured = {}

        def fake_run_simulation(geo, propellant, **kwargs):
            captured.update(kwargs)
            return {'time': np.array([0.0]), 'P_head': np.array([101325.0])}

        monkeypatch.setattr(
            'srm_1d.openmotor_adapter.run_simulation', fake_run_simulation
        )
        monkeypatch.setattr(
            'srm_1d.openmotor_adapter.compute_motor_performance',
            lambda result, nozzle, prop, P_ambient=101325.0: {},
        )
        monkeypatch.setattr(
            'srm_1d.openmotor_adapter.print_performance_summary',
            lambda perf, nozzle: None,
        )

        run_from_ric(str(ric_path), gas_props=SAMPLE_GAS_PROPS)

        chamber = captured['pyrogen_chamber']
        assert chamber.pyrogen.name == 'sibling'
        assert chamber.pyrogen.heat_flux_cal_cm2_s is None
        assert chamber.m_pyrogen_initial > 0.0
        assert captured['T_ignition'] == pytest.approx(756.0)  # v0.7.5 re-LHS default

    def test_run_from_ric_missing_pyrogen_is_informative(self, tmp_path):
        ric_path = tmp_path / 'motor.ric'
        _write_minimal_ric(ric_path)
        with pytest.raises(ValueError, match="Pass pyrogen='bpnv'"):
            run_from_ric(str(ric_path), gas_props=SAMPLE_GAS_PROPS)


class TestPropellantConversion:
    def test_burn_rate_coefficient(self):
        """a should pass through unchanged (same units)."""
        prop = convert_propellant(SAMPLE_RIC_PROPELLANT, SAMPLE_GAS_PROPS)
        assert prop.tabs[0].a == pytest.approx(4.821e-05, rel=1e-6)

    def test_molecular_weight_conversion(self):
        """MW should convert from g/mol to kg/mol."""
        prop = convert_propellant(SAMPLE_RIC_PROPELLANT, SAMPLE_GAS_PROPS)
        assert prop.tabs[0].molecular_weight == pytest.approx(0.0254, rel=1e-6)

    def test_density(self):
        prop = convert_propellant(SAMPLE_RIC_PROPELLANT, SAMPLE_GAS_PROPS)
        assert prop.rho_propellant == pytest.approx(1700.0)

    def test_gamma(self):
        prop = convert_propellant(SAMPLE_RIC_PROPELLANT, SAMPLE_GAS_PROPS)
        assert prop.tabs[0].gamma == pytest.approx(1.19)

    def test_radiation_emissivity_default_is_opt_in(self):
        prop = convert_propellant(SAMPLE_RIC_PROPELLANT, SAMPLE_GAS_PROPS)
        assert prop.radiation_emissivity == pytest.approx(0.0)

        aluminized = dict(SAMPLE_RIC_PROPELLANT)
        aluminized['name'] = 'Hasegawa Propellant 1 (69AP/17HTPB/14Al)'
        prop_al = convert_propellant(aluminized, SAMPLE_GAS_PROPS)
        assert prop_al.radiation_emissivity == pytest.approx(0.0)

    def test_radiation_emissivity_explicit_override(self):
        ric_prop = dict(SAMPLE_RIC_PROPELLANT)
        ric_prop['radiation_emissivity'] = 0.32
        prop = convert_propellant(ric_prop, SAMPLE_GAS_PROPS)
        assert prop.radiation_emissivity == pytest.approx(0.32)

    def test_pressure_range(self):
        """min/max pressures should pass through from .ric tab."""
        prop = convert_propellant(SAMPLE_RIC_PROPELLANT, SAMPLE_GAS_PROPS)
        assert prop.tabs[0].min_pressure == pytest.approx(0.0)
        assert prop.tabs[0].max_pressure == pytest.approx(1e7)

    def test_multi_tab_preserved(self):
        """All tabs should be carried through (no longer dropped)."""
        ric_prop = {
            'density': 1700.0,
            'name': 'two-tab',
            'tabs': [
                {'a': 1e-5, 'n': 0.3, 'k': 1.2, 'm': 25.0, 't': 3000.0,
                 'minPressure': 0.0, 'maxPressure': 5e6},
                {'a': 2e-5, 'n': 0.4, 'k': 1.2, 'm': 25.0, 't': 3000.0,
                 'minPressure': 5e6, 'maxPressure': 1e7},
            ],
        }
        prop = convert_propellant(ric_prop, SAMPLE_GAS_PROPS)
        assert len(prop.tabs) == 2
        assert prop.tabs[0].a == pytest.approx(1e-5)
        assert prop.tabs[1].a == pytest.approx(2e-5)

    def test_transport_from_gas_props(self):
        """Explicit gas props should be used when provided."""
        prop = convert_propellant(SAMPLE_RIC_PROPELLANT, SAMPLE_GAS_PROPS)
        assert prop.mu_gas == pytest.approx(8.842e-5)
        assert prop.k_gas == pytest.approx(0.3685)
        assert prop.Cp_gas == pytest.approx(2060.0)

    def test_missing_transport_hard_faults(self):
        """v0.8.0 D7: no gas props AND no per-tab transport in the .ric
        must hard-fault — srm_1d never fabricates transport."""
        with pytest.raises(ValueError, match="transport"):
            convert_propellant(SAMPLE_RIC_PROPELLANT, None)

    def test_transport_read_from_ric_tabs(self):
        """With per-tab transport present, gas_props=None reads the active
        variant from the .ric (frozen by default)."""
        ric_prop = {
            'name': 'X', 'density': 1700,
            'transportVariant': 'frozen',
            'tabs': [{
                'minPressure': 0, 'maxPressure': 1e7,
                'a': 1e-5, 'n': 0.4, 'k': 1.2, 't': 2800, 'm': 25,
                'mu': 8.842e-5,
                'kThermalFrozen': 0.3685, 'cpFrozen': 2060.0,
                'kThermalEffective': 0.6517, 'cpEffective': 2764.0,
            }],
        }
        prop = convert_propellant(ric_prop, None)
        assert prop.mu_gas == pytest.approx(8.842e-5)
        assert prop.k_gas == pytest.approx(0.3685)   # frozen
        assert prop.Cp_gas == pytest.approx(2060.0)
        # Effective variant selects the other slot
        ric_prop['transportVariant'] = 'effective'
        prop_e = convert_propellant(ric_prop, None)
        assert prop_e.k_gas == pytest.approx(0.6517)
        assert prop_e.Cp_gas == pytest.approx(2764.0)


class TestGeometryConversion:
    def test_segment_count(self):
        geo = convert_geometry(SAMPLE_RIC_GRAINS)
        assert len(geo.segments) == 4

    def test_bore_diameter(self):
        geo = convert_geometry(SAMPLE_RIC_GRAINS)
        assert geo.segments[0].D_bore_fwd == pytest.approx(0.038)

    def test_outer_diameter(self):
        geo = convert_geometry(SAMPLE_RIC_GRAINS)
        assert geo.D_outer == pytest.approx(0.07)

    def test_uninhibited_ends(self):
        geo = convert_geometry(SAMPLE_RIC_GRAINS)
        for seg in geo.segments:
            assert seg.inhibit_fwd is False
            assert seg.inhibit_aft is False

    def test_inhibited_top(self):
        """Top = forward face inhibited."""
        grains = [{'type': 'BATES', 'properties': {
            'coreDiameter': 0.038, 'diameter': 0.07,
            'inhibitedEnds': 'Top', 'length': 0.12}}]
        geo = convert_geometry(grains)
        assert geo.segments[0].inhibit_fwd is True
        assert geo.segments[0].inhibit_aft is False

    def test_inhibited_both(self):
        grains = [{'type': 'BATES', 'properties': {
            'coreDiameter': 0.038, 'diameter': 0.07,
            'inhibitedEnds': 'Both', 'length': 0.12}}]
        geo = convert_geometry(grains)
        assert geo.segments[0].inhibit_fwd is True
        assert geo.segments[0].inhibit_aft is True

    def test_unsupported_grain_raises(self):
        # EndBurningGrain isn't BATES, isn't Conical, and isn't an FMM
        # type — should raise from the adapter dispatch.
        grains = [{'type': 'EndBurningGrain', 'properties': {
            'diameter': 0.080, 'length': 0.100,
            'inhibitedEnds': 'Neither'}}]
        with pytest.raises(ValueError):
            convert_geometry(grains)

    def test_conical_bore_diameters(self):
        grains = [{'type': 'Conical', 'properties': {
            'forwardCoreDiameter': 0.030, 'aftCoreDiameter': 0.050,
            'diameter': 0.080, 'inhibitedEnds': 'Both', 'length': 0.500}}]
        geo = convert_geometry(grains)
        assert geo.segments[0].D_bore_fwd == pytest.approx(0.030)
        assert geo.segments[0].D_bore_aft == pytest.approx(0.050)
        assert geo.D_outer == pytest.approx(0.080)
        assert geo.segments[0].inhibit_fwd is True
        assert geo.segments[0].inhibit_aft is True
        assert geo.segments[0].fmm_table is None

    def test_default_spacing(self):
        """Default gap auto-computed as max(3mm, 5%·D_outer); v0.6.0 snapping
        rounds segment/gap lengths to integer dx, so L_motor matches the
        analytical expectation only within a few-cell tolerance."""
        geo = convert_geometry(SAMPLE_RIC_GRAINS)
        expected_gap = max(0.003, 0.07 * 0.05)
        expected_L = 4 * 0.12 + 5 * expected_gap
        assert abs(geo.L_motor - expected_L) < 5 * geo.dx

    def test_bonded_inhibited_interface_gets_no_default_gap(self):
        """Both inhibited interface faces means adjacent grains are bonded."""
        grains = [
            {'type': 'BATES', 'properties': {
                'coreDiameter': 0.020, 'diameter': 0.050,
                'inhibitedEnds': 'Both', 'length': 0.100}},
            {'type': 'Conical', 'properties': {
                'forwardCoreDiameter': 0.020, 'aftCoreDiameter': 0.030,
                'diameter': 0.050, 'inhibitedEnds': 'Both',
                'length': 0.100}},
        ]
        geo = convert_geometry(grains)
        first_end = geo.segments[0].x_start + geo.segments[0].length
        second_start = geo.segments[1].x_start
        assert first_end == pytest.approx(second_start)

    def test_bonded_interface_right_cell_maps_to_downstream_segment(self):
        """Snapped touching interfaces must not create an epsilon gap/overlap."""
        grains = [
            {'type': 'BATES', 'properties': {
                'coreDiameter': 0.0190500381,
                'diameter': 0.04762509525,
                'inhibitedEnds': 'Both',
                'length': 0.381000762,
            }},
            {'type': 'Conical', 'properties': {
                'forwardCoreDiameter': 0.0190500381,
                'aftCoreDiameter': 0.02857505715,
                'diameter': 0.04762509525,
                'inhibitedEnds': 'Both',
                'length': 0.381000762,
            }},
        ]
        geo = convert_geometry(grains)
        ga = geo.compile_geometry_arrays()
        interface = geo.segments[0].x_start + geo.segments[0].length
        right_idx = int(np.flatnonzero(ga['x_centers'] > interface)[0])
        left_idx = right_idx - 1
        assert ga['cell_segment_id'][left_idx] == 0
        assert ga['cell_segment_id'][right_idx] == 1

    def test_uninhibited_interface_keeps_default_gap(self):
        """A burning interface face needs gas access, so keep a gap."""
        grains = [
            {'type': 'BATES', 'properties': {
                'coreDiameter': 0.020, 'diameter': 0.050,
                'inhibitedEnds': 'Bottom', 'length': 0.100}},
            {'type': 'BATES', 'properties': {
                'coreDiameter': 0.020, 'diameter': 0.050,
                'inhibitedEnds': 'Bottom', 'length': 0.100}},
        ]
        geo = convert_geometry(grains)
        first_end = geo.segments[0].x_start + geo.segments[0].length
        second_start = geo.segments[1].x_start
        assert second_start - first_end > 0.0

    def test_target_propellant_cells(self):
        geo = convert_geometry(SAMPLE_RIC_GRAINS, target_propellant_cells=200)
        # Snapping may shave a few cells from the target depending on
        # gap clamping; must be within ~10% of requested.
        propellant_cells = sum(int(round(s.length / geo.dx))
                               for s in geo.segments)
        assert abs(propellant_cells - 200) < 20


class TestNozzleConversion:
    def test_throat_diameter(self):
        noz = convert_nozzle(SAMPLE_RIC_NOZZLE)
        assert noz.D_throat == pytest.approx(0.02)

    def test_exit_diameter(self):
        noz = convert_nozzle(SAMPLE_RIC_NOZZLE)
        assert noz.D_exit == pytest.approx(0.04)

    def test_erosion_unit_conversion(self):
        """1.5e-10 m/(s·Pa) should become 150 μm/(s·MPa)."""
        noz = convert_nozzle(SAMPLE_RIC_NOZZLE)
        assert noz.erosion_coeff == pytest.approx(150.0, rel=1e-6)

    def test_slag_zero(self):
        noz = convert_nozzle(SAMPLE_RIC_NOZZLE)
        assert noz.slag_coeff == 0.0

    def test_slag_unit_conversion(self):
        """2e-3 (m*Pa)/s should become 2e-9 (m*MPa)/s."""
        ric_nozzle = dict(SAMPLE_RIC_NOZZLE, slagCoeff=2e-3)
        noz = convert_nozzle(ric_nozzle)
        assert noz.slag_coeff == pytest.approx(2e-9)

    def test_divergence_angle(self):
        noz = convert_nozzle(SAMPLE_RIC_NOZZLE)
        assert noz.div_angle == pytest.approx(15.0)

    def test_convergence_angle(self):
        noz = convert_nozzle(SAMPLE_RIC_NOZZLE)
        assert noz.conv_angle == pytest.approx(60.0)

    def test_throat_length(self):
        noz = convert_nozzle(SAMPLE_RIC_NOZZLE)
        assert noz.throat_length == pytest.approx(0.0)

    def test_efficiency(self):
        noz = convert_nozzle(SAMPLE_RIC_NOZZLE)
        assert noz.efficiency == pytest.approx(0.95)


class TestAmbientPressureRouting:
    def test_ambient_pressure_from_ric_config(self):
        """ambPressure should be plumbed through ric_to_sim_args."""
        motor = {
            'propellant': SAMPLE_RIC_PROPELLANT,
            'grains': SAMPLE_RIC_GRAINS,
            'nozzle': SAMPLE_RIC_NOZZLE,
            'config': SAMPLE_RIC_CONFIG,
        }
        args = ric_to_sim_args(motor, gas_props=SAMPLE_GAS_PROPS)
        assert args['P_ambient'] == pytest.approx(101356.5)
        assert 'nozzle' in args
        assert args['nozzle'].D_throat == pytest.approx(0.02)


class TestTabSelection:
    """Hard-switchover tab lookup on Propellant (mirrors openMotor's
    Propellant.getCombustionProperties)."""

    def _two_tab_prop(self):
        ric_prop = {
            'density': 1700.0, 'name': 'two-tab',
            'tabs': [
                {'a': 1e-5, 'n': 0.3, 'k': 1.2, 'm': 25.0, 't': 3000.0,
                 'minPressure': 0.0, 'maxPressure': 5e6},
                {'a': 2e-5, 'n': 0.4, 'k': 1.2, 'm': 25.0, 't': 3000.0,
                 'minPressure': 5e6, 'maxPressure': 1e7},
            ],
        }
        return convert_propellant(ric_prop, SAMPLE_GAS_PROPS)

    def test_select_tab_picks_low_range(self):
        prop = self._two_tab_prop()
        tab = prop.select_tab(3e6)
        assert tab.a == pytest.approx(1e-5)

    def test_select_tab_picks_high_range(self):
        prop = self._two_tab_prop()
        tab = prop.select_tab(7e6)
        assert tab.a == pytest.approx(2e-5)

    def test_select_tab_out_of_range_falls_back(self):
        """Pressure above all tabs → closest-boundary fallback."""
        prop = self._two_tab_prop()
        tab = prop.select_tab(2e7)
        assert tab.a == pytest.approx(2e-5)  # closest to 1e7 boundary

    def test_burn_rate_normal_tab_dependent(self):
        """burn_rate_normal should switch tabs as pressure crosses boundary."""
        prop = self._two_tab_prop()
        r_low = prop.burn_rate_normal(3e6)
        r_high = prop.burn_rate_normal(7e6)
        assert r_low == pytest.approx(1e-5 * 3e6 ** 0.3, rel=1e-6)
        assert r_high == pytest.approx(2e-5 * 7e6 ** 0.4, rel=1e-6)

    def test_representative_tab_widest_range(self):
        prop = self._two_tab_prop()
        rep = prop.representative_tab()
        # Both tabs are 5 MPa wide — first one returned
        assert rep.a == pytest.approx(1e-5)
