import numpy as np
import pytest

from srm_1d import run_simulation
from srm_1d.igniter_plenum import PyrogenChamber
from srm_1d.nozzle import Nozzle
from srm_1d.propellant import Pyrogen
from srm_1d.simulation import (
    _goodman_ignition_sources_and_mass,
    _cal_cm2_s_to_w_m2,
    _pyrogen_surface_heat_power,
    _pyrogen_surface_thermal_sink,
    _thermal_source_power,
)
from tests._motor_fixtures import (
    hasegawa_propellant_1,
    single_cylinder_geo,
)


def _test_chamber():
    pyro = Pyrogen(
        name="phase3-test",
        a=3.0e-5,
        n=0.5,
        rho=1700.0,
        T_flame=2800.0,
        M=0.030,
        gamma=1.25,
        impetus_W=5000.0,
        heat_flux_cal_cm2_s=69.4,
    )
    return PyrogenChamber(
        pyrogen=pyro,
        m_pyrogen_initial=0.003,
        A_burn_initial=5.0e-4,
        A_throat=2.0e-5,
        V_plenum=3.0e-6,
        burn_law="end_burning",
    )


def _small_motor():
    geo = single_cylinder_geo(
        D_bore=0.030, D_outer=0.060, length=0.120,
        target_propellant_cells=12,
    )
    prop = hasegawa_propellant_1()
    nozzle = Nozzle(D_throat=0.010, D_exit=0.020, efficiency=0.95)
    return geo, prop, nozzle


def test_legacy_igniter_kwargs_are_not_accepted():
    geo, prop, nozzle = _small_motor()
    with pytest.raises(TypeError):
        run_simulation(
            geo, prop, nozzle, _test_chamber(),
            igniter_tau=0.1,
            t_max=0.001,
        )


def test_pyrogen_heat_flux_unit_conversion():
    assert _cal_cm2_s_to_w_m2(1.0) == pytest.approx(41840.0)
    assert _cal_cm2_s_to_w_m2(69.4) == pytest.approx(69.4 * 41840.0)


def test_pyrogen_surface_heat_power_uses_measured_flux_when_uncapped():
    power, flux = _pyrogen_surface_heat_power(
        0.1, 2800.0, 300.0, 0.1, 0.01, 2000.0, 1.0e6,
    )
    assert power == pytest.approx(1000.0)
    assert flux == pytest.approx(1.0e6)


def test_pyrogen_surface_heat_power_is_sensible_enthalpy_capped():
    power, flux = _pyrogen_surface_heat_power(
        0.001, 800.0, 300.0, 0.1, 0.01, 2000.0, 1.0e9,
    )
    assert power == pytest.approx(1000.0)
    assert flux == pytest.approx(1.0e6)


def test_pyrogen_surface_thermal_sink_conserves_energy_units():
    """v0.7.1 (Phase 3): sink and pyrogen_enthalpy_source are now in W/m
    (enthalpy injection per unit length). Conservation check: sink * dx
    equals the input power (when the available pyrogen enthalpy exceeds
    the requested sink, no clamping kicks in)."""
    power_w = 1200.0
    sink = _pyrogen_surface_thermal_sink(
        power_w, dx=0.02, pyrogen_enthalpy_source_w_per_m=1.0e9,
    )
    assert sink * 0.02 == pytest.approx(power_w)


def test_thermal_source_power_matches_solver_units():
    """v0.7.1 (Phase 3): thermal_source is now W/m so the power is the
    sum times dx, no Cp factor."""
    thermal_source = np.array([10.0, -2.5])
    assert _thermal_source_power(thermal_source, 0.01, 2) == pytest.approx(0.075)


def test_adjacent_radiation_heats_only_neighbors_and_conserves_sink():
    N = 3
    dx = 0.01
    P = np.full(N, 101325.0)
    T = np.full(N, 300.0)
    T_surf = np.full(N, 293.0)
    delta = np.full(N, 1.0e-6)
    has_ignited = np.array([False, True, False])
    is_burning = np.array([False, True, False])
    is_grain = np.array([True, True, True])
    ignition_time = np.full(N, 1.0e10)
    r_total = np.array([0.0, 0.01, 0.0])
    r_erosive = np.zeros(N)
    mass_source = np.zeros(N)
    thermal_source = np.zeros(N)
    C_burn = np.ones(N)
    endface = np.zeros(N)
    pyrogen_flux = np.zeros(N)
    radiation_flux = np.zeros(N)
    radiation_sink_power = np.zeros(N)
    radiation_emitter = np.zeros(N, dtype=np.bool_)
    x = np.array([0.005, 0.015, 0.025])
    Re = np.zeros(N)
    D_hyd = np.full(N, 0.03)
    f = np.zeros(N)

    mass_source_by_species = np.zeros((N, 3))  # v0.7.1: 3-species
    # v0.7.2 Phase B: spatial coupling args. flame_spread_enabled=False
    # keeps this radiation-focused test independent of the Phase B
    # augmentation (zeros for G_cum + disabled flag => no behavior change
    # vs Phase A).
    flame_spread_augment = np.ones(N)
    # v0.7.1 Phase 3.5: per-species Cp args. Cp_propellant = 2060 (matches
    # the prior Cp_gas value), Cp_pyrogen = 1385 (BPNV-like derived Cp);
    # this test sets mdot_igniter=0 so Cp_pyrogen is unused — value picked
    # for documentation, not for the assertion.
    # v0.7.3 Phase B.2: Y_species (all-ambient → no pyrogen-hot cells).
    # v0.7.3 Phase B.4: zero per-cell pyrogen heat flux array (no DeMar).
    Y_species_test = np.zeros((N, 3))
    Y_species_test[:, 2] = 1.0  # ambient
    pyrogen_heat_flux_arr_test = np.zeros(N)
    out = _goodman_ignition_sources_and_mass(
        P, T, T_surf, delta, has_ignited, is_burning, is_grain,
        ignition_time, r_total, r_erosive,
        mass_source, thermal_source,
        C_burn, endface, pyrogen_flux,
        radiation_flux, radiation_sink_power, radiation_emitter,
        x, Re, D_hyd, f,
        0.0, 1.0e-6, 1700.0, 3041.0, 293.0,
        0.5, 0.3685, 50e-6, 0.45, 1.0e-7, 0.3,
        10000.0, N, dx, 0.0, 300.0,
        2060.0, 1385.0,
        0.0, 0.45, False, 0.0,
        mass_source_by_species,
        flame_spread_augment, False,
        Y_species_test, pyrogen_heat_flux_arr_test,
        False, np.ones(N, dtype=np.bool_), 0,  # v0.7.4 Phase F: gate off, fwd_plenum
    )

    radiation_heat_power = out[4]
    radiation_sink_total_power = out[5]
    assert radiation_flux[0] > 0.0
    assert radiation_flux[1] == pytest.approx(0.0)
    assert radiation_flux[2] > 0.0
    assert radiation_sink_power[1] == pytest.approx(radiation_heat_power)
    assert radiation_sink_total_power == pytest.approx(radiation_heat_power)
    assert radiation_sink_power[0] == pytest.approx(0.0)
    assert radiation_sink_power[2] == pytest.approx(0.0)
    # v0.7.1 Phase 3.5: thermal_source is W/m. Bare combustion contribution
    # is prop_source * T_flame * Cp_propellant; the radiation gas sink debits.
    assert thermal_source[1] < 1700.0 * r_total[1] * C_burn[1] * 3041.0 * 2060.0


def test_adjacent_radiation_sink_can_be_disabled_for_diagnostics():
    N = 3
    dx = 0.01
    P = np.full(N, 101325.0)
    T = np.full(N, 300.0)
    T_surf = np.full(N, 293.0)
    delta = np.full(N, 1.0e-6)
    has_ignited = np.array([False, True, False])
    is_burning = np.array([False, True, False])
    is_grain = np.array([True, True, True])
    ignition_time = np.full(N, 1.0e10)
    r_total = np.array([0.0, 0.01, 0.0])
    r_erosive = np.zeros(N)
    mass_source = np.zeros(N)
    thermal_source = np.zeros(N)
    C_burn = np.ones(N)
    endface = np.zeros(N)
    pyrogen_flux = np.zeros(N)
    radiation_flux = np.zeros(N)
    radiation_sink_power = np.zeros(N)
    radiation_emitter = np.zeros(N, dtype=np.bool_)
    x = np.array([0.005, 0.015, 0.025])
    Re = np.zeros(N)
    D_hyd = np.full(N, 0.03)
    f = np.zeros(N)

    mass_source_by_species = np.zeros((N, 3))  # v0.7.1: 3-species
    # v0.7.2 Phase B: flame_spread_enabled=False keeps this radiation-
    # focused test independent of the Phase B augmentation.
    flame_spread_augment = np.ones(N)
    # v0.7.1 Phase 3.5: per-species Cp args (see sibling test for notes).
    # v0.7.3 Phase B.2/B.4: pass-through args for newly required signature.
    Y_species_test = np.zeros((N, 3))
    Y_species_test[:, 2] = 1.0
    pyrogen_heat_flux_arr_test = np.zeros(N)
    out = _goodman_ignition_sources_and_mass(
        P, T, T_surf, delta, has_ignited, is_burning, is_grain,
        ignition_time, r_total, r_erosive,
        mass_source, thermal_source,
        C_burn, endface, pyrogen_flux,
        radiation_flux, radiation_sink_power, radiation_emitter,
        x, Re, D_hyd, f,
        0.0, 1.0e-6, 1700.0, 3041.0, 293.0,
        0.5, 0.3685, 50e-6, 0.45, 1.0e-7, 0.3,
        10000.0, N, dx, 0.0, 300.0,
        2060.0, 1385.0,
        0.0, 0.45, True, 0.0,
        mass_source_by_species,
        flame_spread_augment, False,
        Y_species_test, pyrogen_heat_flux_arr_test,
        False, np.ones(N, dtype=np.bool_), 0,  # v0.7.4 Phase F: gate off, fwd_plenum
    )

    assert out[4] > 0.0
    assert out[5] == pytest.approx(0.0)
    assert np.max(radiation_flux) > 0.0
    assert np.max(radiation_sink_power) == pytest.approx(0.0)


def test_numerical_collapse_trip_aborts_with_termination_code_4():
    """A forced dt-floor regime trips the classified collapse abort.

    Set dt_max below the trip threshold (1e-9 s) so every step satisfies
    ``dt < 1e-9``. After three consecutive collapsed steps the loop must
    exit with termination_code == 4 instead of running until ``t_max``,
    history-cap, or pressure-cutoff. The classifier should then label
    the run ``numerical_collapse_aborted``.
    """
    from srm_1d.tools.ignition_diagnostics import early_time_diagnostics

    geo, prop, nozzle = _small_motor()
    result = run_simulation(
        geo, prop, nozzle, _test_chamber(),
        T_ignition=850.0,
        t_max=0.001,
        P_cutoff=1.0,
        dt_max=5.0e-10,
        burn_update_interval=1,
        snapshot_interval=0.001,
        cfl_target=0.5,
        verbose=False,
    )

    summary = result['summary']
    assert summary['termination_code'] == 4
    assert summary['termination'] == 'numerical collapse aborted'
    assert summary['steps'] >= 3

    early = early_time_diagnostics(result, window_s=summary['steps'] * 1.0e-9)
    assert early['diagnostic_failure_mode'] == 'numerical_collapse_aborted'
    assert early['termination_code'] == 4


def test_explicit_history_capacity_is_reported_without_changing_termination():
    geo, prop, nozzle = _small_motor()
    result = run_simulation(
        geo, prop, nozzle, _test_chamber(),
        T_ignition=294.0,
        t_max=0.001,
        P_cutoff=1.0,
        history_capacity=20_000,
        snapshot_interval=0.001,
        verbose=False,
    )

    assert result['summary']['history_capacity'] == 20_000
    assert result['summary']['termination'] != 'history array full'


def test_explicit_history_capacity_must_be_positive():
    geo, prop, nozzle = _small_motor()
    with pytest.raises(ValueError, match="history_capacity must be positive"):
        run_simulation(
            geo, prop, nozzle, _test_chamber(),
            history_capacity=0,
            verbose=False,
        )


def test_burn_update_interval_must_be_positive():
    geo, prop, nozzle = _small_motor()
    with pytest.raises(ValueError, match="burn_update_interval must be positive"):
        run_simulation(
            geo, prop, nozzle, _test_chamber(),
            burn_update_interval=0,
            verbose=False,
        )


def test_geometry_update_interval_must_be_positive():
    geo, prop, nozzle = _small_motor()
    with pytest.raises(
        ValueError, match="geometry_update_interval must be positive"
    ):
        run_simulation(
            geo, prop, nozzle, _test_chamber(),
            geometry_update_interval=0,
            verbose=False,
        )


def test_run_exposes_per_cell_numerical_limit_diagnostics():
    geo, prop, nozzle = _small_motor()
    result = run_simulation(
        geo, prop, nozzle, _test_chamber(),
        T_ignition=294.0,
        t_max=0.001,
        P_cutoff=1.0,
        snapshot_interval=0.001,
        verbose=False,
    )

    limits = result['numerical_limits']
    assert result['summary']['burn_update_interval'] == max(
        10, geo.N_cells // 5
    )
    assert result['summary']['geometry_update_interval'] == (
        result['summary']['burn_update_interval']
    )
    assert limits['x_m'].shape == (geo.N_cells,)
    for name in ('temperature_floor', 'temperature_ceiling', 'pressure_floor'):
        assert limits[name]['activation_count_by_cell'].shape == (geo.N_cells,)
        assert limits[name]['duration_s_by_cell'].shape == (geo.N_cells,)
        assert limits[name]['first_activation_time_s_by_cell'].shape == (
            geo.N_cells,
        )
        assert limits[name]['last_activation_time_s_by_cell'].shape == (
            geo.N_cells,
        )
        assert np.all(limits[name]['activation_count_by_cell'] >= 0)
    mach_limit = limits['port_mach_cap']
    assert mach_limit['threshold_mach'] == 0.0
    assert mach_limit['x_m'].shape == (geo.N_cells - 1,)
    assert np.count_nonzero(mach_limit['activation_count_by_face']) == 0
    assert np.all(np.isnan(mach_limit['first_activation_time_s_by_face']))


def test_enabled_port_mach_cap_reports_face_activations():
    geo, prop, nozzle = _small_motor()
    result = run_simulation(
        geo, prop, nozzle, _test_chamber(),
        T_ignition=294.0,
        port_mach_cap=1.0e-6,
        t_max=0.001,
        P_cutoff=1.0,
        snapshot_interval=0.001,
        verbose=False,
    )

    diagnostic = result['numerical_limits']['port_mach_cap']
    active = diagnostic['activation_count_by_face'] > 0
    assert np.any(active)
    assert np.all(np.isfinite(diagnostic['first_activation_time_s_by_face'][active]))
    assert np.all(
        diagnostic['last_activation_time_s_by_face'][active]
        >= diagnostic['first_activation_time_s_by_face'][active]
    )


def test_pyrogen_driven_run_reports_ignition_and_pyrogen_state():
    geo, prop, nozzle = _small_motor()
    result = run_simulation(
        geo, prop, nozzle, _test_chamber(),
        T_ignition=294.0,
        t_max=0.010,
        P_cutoff=1.0,
        dt_max=2.0e-5,
        burn_update_interval=1,
        snapshot_interval=0.002,
        cfl_target=0.5,
    )

    summary = result['summary']
    assert summary['pyrogen_mass_burned'] > 0.0
    assert summary['pyrogen_peak_P'] > 101325.0
    assert np.max(result['P_head']) > 101325.0
    assert np.max(result['mdot_ig']) > 0.0
    assert len(result['P_ig']) == len(result['time'])

    assert result['snapshots']
    last = result['snapshots'][-1]
    assert 'T_surf' in last
    assert 'pyrogen_surface_heat_flux' in last
    assert 'radiation_heat_flux' in last
    assert np.any(last['is_burning'])
    assert np.max(last['T_surf']) >= 294.0
    assert np.max(last['pyrogen_surface_heat_flux']) >= 0.0
    assert result['summary']['radiation_emissivity'] == pytest.approx(prop.radiation_emissivity)
    assert 'pyrogen_enthalpy_power' in result
    assert 'normal_sidewall_thermal_power' in result
    assert 'erosive_sidewall_thermal_power' in result
    assert 'endface_thermal_power' in result
    assert 'pyrogen_gas_thermal_power' in result
    assert 'gas_sensible_energy' in result
    assert 'gas_sensible_dE_dt' in result
    assert 'convective_scalar_flux_power' in result
    assert 'nozzle_scalar_flux_power' in result
    assert 'clipping_correction_power' in result
    assert 'gas_surface_heat_sink_power' in result
    assert 'convective_wall_heat_sink_power' in result
    assert 'energy_residual' in result
    assert 'pyrogen_momentum_residual' in result
    assert 'dt' in result
    assert 'n_burning' in result
    assert 'n_ignited' in result
    assert 'radiation_emitter_count' in result
    assert 'radiation_receiver_count' in result
    assert 'max_gas_temperature' in result
    assert 'max_surface_temperature' in result
    assert 'max_mach' in result
    assert len(result['dt']) == len(result['time'])
    assert summary['termination_code'] >= 0
    assert summary['history_capacity'] >= summary['steps']
    assert 'first_ignition_time_s' in summary
    assert 'first_ignition_cell' in summary
    assert 'energy_residual_convention' in summary
    assert summary['energy_convention']['classification'] == (
        'sensible-enthalpy scalar transport balance'
    )
