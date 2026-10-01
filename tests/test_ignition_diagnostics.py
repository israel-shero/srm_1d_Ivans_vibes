import copy
from types import SimpleNamespace

import numpy as np
import pytest

from srm_1d import run_simulation
from examples.ignition_spike_diagnostic import (
    DEFAULT_VARIANTS,
    RADIATION_COLLAPSE_DEFAULT_VARIANTS,
    RADIATION_COLLAPSE_VARIANTS,
    RADIATION_PROBE_DEFAULT_VARIANTS,
    RADIATION_PROBE_VARIANTS,
    VARIANTS,
    build_parser,
)
from srm_1d.igniter_plenum import PyrogenChamber
from srm_1d.nozzle import Nozzle
from srm_1d.propellant import Pyrogen
from tests._motor_fixtures import (
    example_bates_geo,
    hasegawa_propellant_1,
    single_cylinder_geo,
)
from srm_1d.tools.ignition_diagnostics import (
    analyze_ignition_spike,
    collapse_event_trace,
    early_time_diagnostics,
    energy_momentum_timeseries,
    ignition_spread_metrics,
    pressure_landmarks,
    sensible_enthalpy_audit,
    step_diagnostics_timeseries,
    source_timeseries,
)


def _snap(t, burning, r_total=None, r_erosive=None, endface=None,
          pyrogen_surface_heat_flux=None, radiation_heat_flux=None,
          mach=None, mass_source=None, thermal_source=None):
    burning = np.asarray(burning, dtype=bool)
    n = len(burning)
    if r_total is None:
        r_total = np.zeros(n)
    if r_erosive is None:
        r_erosive = np.zeros(n)
    if endface is None:
        endface = np.zeros(n)
    if pyrogen_surface_heat_flux is None:
        pyrogen_surface_heat_flux = np.zeros(n)
    if radiation_heat_flux is None:
        radiation_heat_flux = np.zeros(n)
    if mach is None:
        mach = np.zeros(n)
    if mass_source is None:
        mass_source = np.zeros(n)
    if thermal_source is None:
        thermal_source = np.zeros(n)
    return {
        "t": float(t),
        "x": np.arange(n, dtype=float),
        "P": np.full(n, 101325.0),
        "T": np.full(n, 300.0),
        "T_surf": np.full(n, 293.0),
        "is_burning": burning,
        "is_grain": np.ones(n, dtype=bool),
        "C_burn": np.ones(n),
        "r_total": np.asarray(r_total, dtype=float),
        "r_erosive": np.asarray(r_erosive, dtype=float),
        "endface_msource": np.asarray(endface, dtype=float),
        "pyrogen_surface_heat_flux": np.asarray(pyrogen_surface_heat_flux, dtype=float),
        "radiation_heat_flux": np.asarray(radiation_heat_flux, dtype=float),
        "Mach": np.asarray(mach, dtype=float),
        "mass_source": np.asarray(mass_source, dtype=float),
        "thermal_source": np.asarray(thermal_source, dtype=float),
    }


def _synthetic_result(time, pressure, mdot, snapshots=()):
    time = np.asarray(time, dtype=float)
    return {
        "time": time,
        "P_head": np.asarray(pressure, dtype=float),
        "P_exit": np.asarray(pressure, dtype=float),
        "mdot_ig": np.asarray(mdot, dtype=float),
        "P_ig": np.asarray(mdot, dtype=float) * 1.0e6 + 101325.0,
        "m_pyrogen": np.maximum(0.0, 0.001 - np.cumsum(mdot) * 0.01),
        "snapshots": list(snapshots),
        "P_ambient": 101325.0,
        "summary": {},
    }


def _synthetic_probe_result(**overrides):
    time = np.asarray(overrides.pop("time", [0.0, 0.001, 0.002]), dtype=float)
    n = time.size
    result = _synthetic_result(
        time,
        overrides.pop("pressure", np.linspace(101325.0, 2.0e6, n)),
        overrides.pop("mdot", np.zeros(n)),
    )
    result.update({
        "dt": np.asarray(overrides.pop("dt", np.full(n, 0.001)), dtype=float),
        "massflow": np.asarray(overrides.pop("massflow", np.zeros(n)), dtype=float),
        "n_burning": np.asarray(overrides.pop("n_burning", np.ones(n)), dtype=float),
        "n_ignited": np.asarray(overrides.pop("n_ignited", np.ones(n)), dtype=float),
        "radiation_emitter_count": np.asarray(
            overrides.pop("radiation_emitter_count", np.zeros(n)), dtype=float
        ),
        "radiation_receiver_count": np.asarray(
            overrides.pop("radiation_receiver_count", np.zeros(n)), dtype=float
        ),
        "radiation_heat_power": np.asarray(
            overrides.pop("radiation_heat_power", np.zeros(n)), dtype=float
        ),
        "radiation_sink_power": np.asarray(
            overrides.pop("radiation_sink_power", np.zeros(n)), dtype=float
        ),
        "clipping_correction_power": np.asarray(
            overrides.pop("clipping_correction_power", np.zeros(n)), dtype=float
        ),
        "thermal_source_power": np.asarray(
            overrides.pop("thermal_source_power", np.ones(n)), dtype=float
        ),
        "convective_scalar_flux_power": np.asarray(
            overrides.pop("convective_scalar_flux_power", np.ones(n)), dtype=float
        ),
        "erosive_sidewall_thermal_power": np.asarray(
            overrides.pop("erosive_sidewall_thermal_power", np.zeros(n)), dtype=float
        ),
        "min_gas_temperature": np.asarray(
            overrides.pop("min_gas_temperature", np.full(n, 293.0)), dtype=float
        ),
        "max_gas_temperature": np.asarray(
            overrides.pop("max_gas_temperature", np.full(n, 900.0)), dtype=float
        ),
        "min_surface_temperature": np.asarray(
            overrides.pop("min_surface_temperature", np.full(n, 293.0)), dtype=float
        ),
        "max_surface_temperature": np.asarray(
            overrides.pop("max_surface_temperature", np.full(n, 500.0)), dtype=float
        ),
        "min_pressure": np.asarray(
            overrides.pop("min_pressure", np.full(n, 101325.0)), dtype=float
        ),
        "max_pressure": np.asarray(
            overrides.pop("max_pressure", np.linspace(101325.0, 2.0e6, n)),
            dtype=float,
        ),
        "max_mach": np.asarray(overrides.pop("max_mach", np.zeros(n)), dtype=float),
        "ignition_time_by_cell": np.asarray(
            overrides.pop("ignition_time_by_cell", [0.001, 1.0e10]), dtype=float
        ),
        "summary": {
            "termination": overrides.pop("termination", "t_max reached"),
            "termination_code": overrides.pop("termination_code", 0),
            "history_cap_reached": overrides.pop("history_cap_reached", False),
            "history_capacity": overrides.pop("history_capacity", 100),
            "steps": n,
            "active_radiation_emissivity": overrides.pop(
                "active_radiation_emissivity", 0.0
            ),
        },
    })
    result.update(overrides)
    return result


def _test_chamber(**pyrogen_overrides):
    pyro = Pyrogen(
        name="diagnostic-test",
        a=3.0e-5,
        n=0.5,
        rho=1700.0,
        T_flame=2800.0,
        M=0.030,
        gamma=1.25,
        impetus_W=5000.0,
        heat_flux_cal_cm2_s=69.4,
    )
    for key, value in pyrogen_overrides.items():
        setattr(pyro, key, value)
    return PyrogenChamber(
        pyrogen=pyro,
        m_pyrogen_initial=0.003,
        A_burn_initial=5.0e-4,
        A_throat=2.0e-5,
        V_plenum=3.0e-6,
        burn_law="end_burning",
    )


def _small_motor(inhibited=True):
    if inhibited:
        geo = single_cylinder_geo(
            D_bore=0.030, D_outer=0.060, length=0.120,
            target_propellant_cells=12,
        )
    else:
        geo = example_bates_geo(target_propellant_cells=32)
    prop = hasegawa_propellant_1()
    nozzle = Nozzle(D_throat=0.010, D_exit=0.020, efficiency=0.95)
    return geo, prop, nozzle


def _run_small(**kwargs):
    geo, prop, nozzle = _small_motor(kwargs.pop("inhibited", True))
    result = run_simulation(
        geo, prop, nozzle, _test_chamber(),
        T_ignition=kwargs.pop("T_ignition", 294.0),
        t_max=kwargs.pop("t_max", 0.006),
        P_cutoff=kwargs.pop("P_cutoff", 1.0),
        dt_max=kwargs.pop("dt_max", 2.0e-5),
        burn_update_interval=kwargs.pop("burn_update_interval", 1),
        snapshot_interval=kwargs.pop("snapshot_interval", 0.001),
        cfl_target=kwargs.pop("cfl_target", 0.5),
        verbose=False,
        **kwargs,
    )
    return result, geo, prop


def test_classifies_pyrogen_driven_spike_when_peak_overlaps_mdot():
    result = _synthetic_result(
        [0.0, 0.1, 0.2],
        [101325.0, 3.0e6, 1.5e6],
        [0.0, 0.02, 0.01],
    )
    diagnostics = analyze_ignition_spike(result)
    assert diagnostics["classification"]["primary_driver"] == "pyrogen_combustion"
    assert diagnostics["classification"]["pyrogen_combustion"]
    assert "energy" in diagnostics


def test_diagnostic_cli_default_variants_are_registered():
    args = build_parser().parse_args([])
    assert args.variants == DEFAULT_VARIANTS
    assert set(args.variants).issubset(VARIANTS)


def test_radiation_probe_default_variants_are_registered():
    args = build_parser().parse_args(["--mode", "radiation-probe"])
    assert args.variants == DEFAULT_VARIANTS
    assert set(RADIATION_PROBE_DEFAULT_VARIANTS).issubset(RADIATION_PROBE_VARIANTS)


def test_radiation_collapse_default_variants_are_registered():
    args = build_parser().parse_args(["--mode", "radiation-collapse"])
    assert args.variants == DEFAULT_VARIANTS
    assert set(RADIATION_COLLAPSE_DEFAULT_VARIANTS).issubset(
        RADIATION_COLLAPSE_VARIANTS
    )


def test_step_diagnostics_timeseries_preserves_probe_histories():
    result = _synthetic_probe_result(
        dt=[1.0e-4, 2.0e-4],
        time=[0.0, 1.0e-4],
        n_burning=[0, 1],
        radiation_receiver_count=[0, 1],
        max_mach=[0.0, 0.25],
    )
    step = step_diagnostics_timeseries(result)
    np.testing.assert_allclose(step["dt"], [1.0e-4, 2.0e-4])
    np.testing.assert_allclose(step["n_burning"], [0, 1])
    np.testing.assert_allclose(step["radiation_receiver_count"], [0, 1])
    np.testing.assert_allclose(step["max_mach"], [0.0, 0.25])


def test_early_diagnostics_classifies_normal_completed_run():
    result = _synthetic_probe_result()
    early = early_time_diagnostics(result)
    assert early["diagnostic_failure_mode"] == "normal_completed"
    assert not early["history_cap_reached"]
    assert not early["radiation_enabled_zero_activity"]


def test_early_diagnostics_flags_history_cap_or_step_limit():
    result = _synthetic_probe_result(
        termination="history array full",
        termination_code=3,
        history_cap_reached=True,
        history_capacity=3,
    )
    early = early_time_diagnostics(result)
    assert early["history_cap_reached"]
    assert early["diagnostic_failure_mode"] == "history_cap_or_step_limit"


def test_early_diagnostics_flags_radiation_enabled_zero_activity():
    result = _synthetic_probe_result(
        active_radiation_emissivity=0.45,
        n_burning=[0, 0, 0],
        n_ignited=[0, 0, 0],
        ignition_time_by_cell=[1.0e10, 1.0e10],
    )
    early = early_time_diagnostics(result)
    assert early["radiation_enabled_zero_activity"]
    assert early["diagnostic_failure_mode"] == "radiation_enabled_zero_activity"


def test_early_diagnostics_flags_clipping_dominated_energy():
    result = _synthetic_probe_result(
        clipping_correction_power=[-100.0, -100.0, -100.0],
        thermal_source_power=[1.0, 1.0, 1.0],
        convective_scalar_flux_power=[1.0, 1.0, 1.0],
    )
    early = early_time_diagnostics(result)
    assert early["clipping_dominated_energy"]
    assert early["diagnostic_failure_mode"] == "timestep_front_numerical_pathology"
    assert early["collapse_class"] == "collapse"


def test_exact_ignition_spread_uses_ignition_time_by_cell():
    result = _synthetic_result(
        [0.0, 0.001, 0.002],
        [101325.0, 2.0e6, 3.0e6],
        [0.0, 0.0, 0.0],
        snapshots=[_snap(0.0, [False] * 10)],
    )
    result["ignition_time_by_cell"] = np.array([
        0.001, 0.002, 0.003, 0.004, 0.005,
        0.006, 0.007, 0.008, 0.009, 0.010,
    ])

    spread = ignition_spread_metrics(result)

    assert spread["exact_spread_metrics"]
    assert spread["spread_metric_source"] == "ignition_time_by_cell"
    assert spread["first_ignition_time_s"] == pytest.approx(0.001)
    assert spread["t10_s"] == pytest.approx(0.001)
    assert spread["t50_s"] == pytest.approx(0.005)
    assert spread["t90_s"] == pytest.approx(0.009)
    assert spread["t100_s"] == pytest.approx(0.010)


def test_collapse_detector_reports_threshold_events():
    result = _synthetic_probe_result(
        time=[0.0, 1.0e-4, 2.0e-4, 3.0e-4],
        dt=[1.0e-4, 5.0e-9, 5.0e-9, 5.0e-9],
        max_pressure=[101325.0, 2.0e6, 150.0e6, 150.0e6],
        max_mach=[0.0, 10.0, 20.0, 2000.0],
        clipping_correction_power=[0.0, 0.0, -100.0, -100.0],
        thermal_source_power=[1.0, 1.0, 1.0, 1.0],
        convective_scalar_flux_power=[1.0, 1.0, 1.0, 1.0],
        active_radiation_emissivity=0.45,
        erosive_sidewall_thermal_power=[0.0, 1.0, 1.0, 1.0],
    )

    early = early_time_diagnostics(result)

    assert early["collapse_detected"]
    assert early["collapse_class"] == "collapse"
    assert early["first_dt_collapse_time_s"] == pytest.approx(1.0e-4)
    assert early["first_pressure_collapse_time_s"] == pytest.approx(2.0e-4)
    assert early["first_mach_collapse_time_s"] == pytest.approx(3.0e-4)
    assert early["first_clipping_dominated_time_s"] == pytest.approx(2.0e-4)
    assert early["first_collapse_time_s"] == pytest.approx(1.0e-4)
    assert (
        early["collapse_branch_suspect"]
        == "piso_nozzle_front_numerical_instability"
    )


def test_collapse_event_trace_records_local_cells_and_limiter_flags():
    snap = _snap(
        0.001,
        [True, True, True, True, True],
        r_total=[0.1, 0.2, 0.3, 0.4, 0.5],
        r_erosive=[0.0, 0.1, 0.2, 0.3, 0.4],
        radiation_heat_flux=[0.0, 0.0, 50.0, 0.0, 0.0],
        mach=[1.0, 2.0, 5.0, 3.0, 2.0],
        mass_source=[10.0, 20.0, 30.0, 40.0, 50.0],
        thermal_source=[100.0, 200.0, 300.0, 400.0, 500.0],
    )
    result = _synthetic_probe_result(
        time=[0.0, 0.001, 0.002],
        dt=[1.0e-4, 1.0e-9, 1.0e-9],
        snapshots=[snap],
        max_pressure=[101325.0, 150.0e6, 150.0e6],
        max_mach=[0.0, 2000.0, 2000.0],
        radiation_heat_power=[0.0, 1.0, 1.0],
        radiation_receiver_count=[0.0, 1.0, 1.0],
        clipping_correction_power=[0.0, -100.0, -100.0],
        thermal_source_power=[1.0, 1.0, 1.0],
        convective_scalar_flux_power=[1.0, 1.0, 1.0],
        ignition_time_by_cell=[0.001, 0.001, 0.001, 0.001, 0.001],
        active_radiation_emissivity=0.45,
    )
    spread = ignition_spread_metrics(result)
    early = early_time_diagnostics(result, spread=spread)

    rows = collapse_event_trace(result, early, spread)

    radiation_rows = [row for row in rows if row["event"] == "first_radiation"]
    assert radiation_rows
    assert {row["cell"] for row in radiation_rows} == {0, 1, 2, 3, 4}
    assert all(row["event_cell"] == 2 for row in radiation_rows)
    center = next(row for row in radiation_rows if row["cell"] == 2)
    assert center["radiation_heat_flux_w_m2"] == pytest.approx(50.0)
    assert center["mass_source_kg_m_s"] == pytest.approx(30.0)
    assert center["thermal_source_k_kg_m_s"] == pytest.approx(300.0)
    assert center["dt_below_1e_8"]
    assert center["pressure_above_100mpa"]
    assert center["mach_above_1e3"]
    assert center["clipping_dominated_step"]


def test_classifies_erosive_snap_on_after_pyrogen_burnout():
    result = _synthetic_result(
        [0.0, 0.1, 0.2],
        [101325.0, 3.0e6, 1.0e6],
        [0.02, 0.0, 0.0],
        snapshots=[
            _snap(0.1, [True, True], r_total=[10.0, 10.0], r_erosive=[9.0, 9.0]),
            _snap(0.2, [True, True], r_total=[1.0, 1.0], r_erosive=[0.0, 0.0]),
        ],
    )
    diagnostics = analyze_ignition_spike(
        result, propellant=SimpleNamespace(rho_propellant=1.0)
    )
    assert diagnostics["classification"]["primary_driver"] == "erosive_snap_on"
    assert diagnostics["classification"]["erosive_snap_on"]


def test_pressure_landmarks_label_startup_peak_separately_from_global_peak():
    result = _synthetic_result(
        [0.0, 0.05, 0.2, 2.0],
        [101325.0, 1.0e6, 0.8e6, 5.0e6],
        [0.02, 0.0, 0.0, 0.0],
    )
    pressure = pressure_landmarks(result, startup_margin_s=0.10)
    assert pressure["startup_window_peak_time_s"] == pytest.approx(0.05)
    assert pressure["startup_window_peak_pressure_mpa"] == pytest.approx(1.0)
    assert pressure["global_peak_time_s"] == pytest.approx(2.0)
    assert pressure["global_peak_pressure_mpa"] == pytest.approx(5.0)
    assert pressure["peak_time_s"] == pressure["startup_window_peak_time_s"]


def test_detects_instant_ignition_collapse():
    result = _synthetic_result(
        [0.0, 0.1, 0.2],
        [101325.0, 1.0e6, 0.8e6],
        [0.0, 0.0, 0.0],
        snapshots=[
            _snap(0.0, [False, False, False, False]),
            _snap(0.1, [True, True, True, True]),
            _snap(0.2, [True, True, True, True]),
        ],
    )
    spread = ignition_spread_metrics(result)
    assert spread["instant_ignition_collapse"]
    assert spread["spread_10_90_s"] == pytest.approx(0.0)


def test_source_timeseries_estimates_sources_without_mutating_result():
    result = _synthetic_result(
        [0.0, 0.1],
        [101325.0, 2.0e6],
        [0.0, 1.0],
        snapshots=[
            _snap(0.1, [True, True], r_total=[3.0, 3.0],
                  r_erosive=[1.0, 1.0], endface=[0.5, 0.5]),
        ],
    )
    original = copy.deepcopy(result)
    sources = source_timeseries(
        result, propellant=SimpleNamespace(rho_propellant=2.0)
    )
    assert sources["normal_sidewall_kg_s"][0] == pytest.approx(8.0)
    assert sources["erosive_sidewall_kg_s"][0] == pytest.approx(4.0)
    assert sources["endface_kg_s"][0] == pytest.approx(1.0)
    assert sources["pyrogen_kg_s"][0] == pytest.approx(1.0)
    assert sources["pyrogen_surface_heat_power_w"][0] == pytest.approx(0.0)
    assert sources["pyrogen_surface_heat_flux_w_m2"][0] == pytest.approx(0.0)
    assert sources["radiation_heat_power_w"][0] == pytest.approx(0.0)
    assert sources["radiation_heat_flux_w_m2"][0] == pytest.approx(0.0)
    np.testing.assert_array_equal(result["P_head"], original["P_head"])
    np.testing.assert_array_equal(
        result["snapshots"][0]["r_total"],
        original["snapshots"][0]["r_total"],
    )


def test_default_diagnostic_controls_preserve_behavior():
    base, _geo, _prop = _run_small()
    explicit, _geo, _prop = _run_small(
        initial_gas_temperature=None,
        diagnostic_disable_erosive=False,
        diagnostic_disable_endfaces=False,
        diagnostic_disable_momentum=False,
        diagnostic_disable_pyrogen_surface_heating=False,
        diagnostic_disable_adjacent_radiation=False,
        igniter_axial_momentum_fraction=1.0,
    )
    np.testing.assert_allclose(base["P_head"], explicit["P_head"])
    np.testing.assert_allclose(base["mdot_ig"], explicit["mdot_ig"])


def test_source_timeseries_reports_pyrogen_surface_heating_power():
    result = _synthetic_result(
        [0.0, 0.1],
        [101325.0, 2.0e6],
        [0.0, 1.0],
        snapshots=[
            _snap(
                0.1, [False, False],
                pyrogen_surface_heat_flux=[1000.0, 500.0],
            ),
        ],
    )
    sources = source_timeseries(result)
    assert sources["pyrogen_surface_heat_power_w"][0] == pytest.approx(1500.0)
    assert sources["pyrogen_surface_heat_flux_w_m2"][0] == pytest.approx(1000.0)


def test_source_timeseries_reports_adjacent_radiation_heating_power():
    result = _synthetic_result(
        [0.0, 0.1],
        [101325.0, 2.0e6],
        [0.0, 1.0],
        snapshots=[
            _snap(
                0.1, [False, False],
                radiation_heat_flux=[2000.0, 1000.0],
            ),
        ],
    )
    sources = source_timeseries(result)
    assert sources["radiation_heat_power_w"][0] == pytest.approx(3000.0)
    assert sources["radiation_heat_flux_w_m2"][0] == pytest.approx(2000.0)


def test_energy_momentum_timeseries_preserves_result_ledgers():
    result = _synthetic_result(
        [0.0, 0.1],
        [101325.0, 2.0e6],
        [0.0, 1.0],
    )
    result["pyrogen_enthalpy_power"] = np.array([10.0, 20.0])
    result["normal_sidewall_thermal_power"] = np.array([1.0, 2.0])
    result["convective_scalar_flux_power"] = np.array([-3.0, -4.0])
    result["clipping_correction_power"] = np.array([0.0, -5.0])
    result["gas_sensible_energy"] = np.array([100.0, 140.0])
    result["gas_sensible_dE_dt"] = np.array([0.0, 400.0])
    result["pyrogen_momentum_residual"] = np.array([0.0, 1.0e-12])
    audit = energy_momentum_timeseries(result)
    np.testing.assert_allclose(audit["pyrogen_enthalpy_power"], [10.0, 20.0])
    np.testing.assert_allclose(audit["normal_sidewall_thermal_power"], [1.0, 2.0])
    np.testing.assert_allclose(audit["convective_scalar_flux_power"], [-3.0, -4.0])
    np.testing.assert_allclose(audit["clipping_correction_power"], [0.0, -5.0])
    np.testing.assert_allclose(audit["gas_sensible_energy"], [100.0, 140.0])
    np.testing.assert_allclose(audit["gas_sensible_dE_dt"], [0.0, 400.0])
    np.testing.assert_allclose(audit["pyrogen_momentum_residual"], [0.0, 1.0e-12])
    np.testing.assert_allclose(audit["energy_residual"], [0.0, 0.0])


def test_sensible_enthalpy_audit_reconstructs_source_and_storage_jumps():
    result = _synthetic_result([0.0, 0.1], [101325.0, 2.0e6], [0.0, 1.0])
    result.update({
        "normal_sidewall_thermal_power": np.array([10.0, 20.0]),
        "erosive_sidewall_thermal_power": np.array([1.0, 2.0]),
        "endface_thermal_power": np.array([3.0, 4.0]),
        "pyrogen_enthalpy_power": np.array([5.0, 6.0]),
        "gas_surface_heat_sink_power": np.array([2.0, 3.0]),
        "radiation_sink_power": np.array([1.0, 1.0]),
        "convective_wall_heat_sink_power": np.array([4.0, 5.0]),
        "thermal_source_power": np.array([12.0, 23.0]),
        "convective_scalar_flux_power": np.array([-7.0, -8.0]),
        "nozzle_enthalpy_power": np.array([7.0, 8.0]),
        "gas_sensible_energy_before": np.array([100.0, 131.0]),
        "gas_sensible_energy": np.array([130.0, 150.0]),
    })

    audit = sensible_enthalpy_audit(result)

    np.testing.assert_allclose(audit["source_reconstructed_power"], [12.0, 23.0])
    np.testing.assert_allclose(audit["source_ledger_mismatch_power"], [0.0, 0.0])
    np.testing.assert_allclose(audit["nozzle_boundary_mismatch_power"], [0.0, 0.0])
    np.testing.assert_allclose(audit["between_step_storage_jump"], [0.0, 1.0])


def test_ambient_initial_gas_ignites_from_pyrogen_surface_heating():
    ambient, _geo, prop = _run_small(
        T_ignition=500.0,
        initial_gas_temperature=propellant_initial_temperature(),
    )
    disabled, _geo, _prop = _run_small(
        T_ignition=500.0,
        initial_gas_temperature=propellant_initial_temperature(),
        diagnostic_disable_pyrogen_surface_heating=True,
    )
    ambient_spread = ignition_spread_metrics(ambient)
    disabled_spread = ignition_spread_metrics(disabled)
    assert ambient["summary"]["initial_gas_temperature"] == pytest.approx(293.0)
    assert np.isfinite(ambient_spread["first_ignition_time_s"])
    assert max(
        np.max(s["pyrogen_surface_heat_flux"]) for s in ambient["snapshots"]
    ) > 0.0
    assert max(
        np.max(s["pyrogen_surface_heat_flux"]) for s in disabled["snapshots"]
    ) == pytest.approx(0.0)
    heated_cells = set()
    for snap in ambient["snapshots"]:
        heated = np.flatnonzero(snap["pyrogen_surface_heat_flux"] > 0.0)
        if heated.size:
            heated_cells.update(int(i) for i in heated)
    assert len(heated_cells) == 1
    assert (
        np.isnan(disabled_spread["first_ignition_time_s"])
        or disabled_spread["first_ignition_time_s"] >= ambient_spread["first_ignition_time_s"]
    )


def propellant_initial_temperature():
    return hasegawa_propellant_1().T_initial


def test_custom_pyrogen_missing_heat_flux_requires_explicit_disable():
    geo, prop, nozzle = _small_motor()
    chamber = _test_chamber(heat_flux_cal_cm2_s=None)
    with pytest.raises(ValueError, match="heat_flux_cal_cm2_s"):
        run_simulation(
            geo, prop, nozzle, chamber,
            T_ignition=500.0,
            initial_gas_temperature=prop.T_initial,
            t_max=0.001,
            P_cutoff=1.0,
            dt_max=2.0e-5,
            burn_update_interval=1,
            snapshot_interval=0.001,
            verbose=False,
        )

    result = run_simulation(
        geo, prop, nozzle, chamber,
        T_ignition=500.0,
        initial_gas_temperature=prop.T_initial,
        diagnostic_disable_pyrogen_surface_heating=True,
        t_max=0.001,
        P_cutoff=1.0,
        dt_max=2.0e-5,
        burn_update_interval=1,
        snapshot_interval=0.001,
        verbose=False,
    )
    assert result["summary"]["diagnostic_disable_pyrogen_surface_heating"]


def test_disabling_erosive_preserves_normal_burn_and_zeroes_erosive_increment():
    baseline, _geo, _prop = _run_small(roughness=100e-6, T_ignition=294.0)
    disabled, _geo, _prop = _run_small(
        roughness=100e-6,
        T_ignition=294.0,
        diagnostic_disable_erosive=True,
    )
    baseline_r_erosive = max(np.max(s["r_erosive"]) for s in baseline["snapshots"])
    disabled_r_erosive = max(np.max(s["r_erosive"]) for s in disabled["snapshots"])
    disabled_r_total = max(np.max(s["r_total"]) for s in disabled["snapshots"])
    assert baseline_r_erosive > 0.0
    assert disabled_r_erosive == pytest.approx(0.0)
    assert disabled_r_total > 0.0


def test_disabling_endfaces_removes_endface_source_without_changing_sidewall_geometry():
    baseline, _geo, _prop = _run_small(inhibited=False, T_ignition=294.0)
    disabled, _geo, _prop = _run_small(
        inhibited=False,
        T_ignition=294.0,
        diagnostic_disable_endfaces=True,
    )
    assert max(np.max(s["endface_msource"]) for s in baseline["snapshots"]) > 0.0
    assert max(np.max(s["endface_msource"]) for s in disabled["snapshots"]) == pytest.approx(0.0)
    np.testing.assert_allclose(
        baseline["snapshots"][0]["C_burn"],
        disabled["snapshots"][0]["C_burn"],
        atol=1.0e-6,
    )


def test_disabling_momentum_removes_igniter_momentum_source():
    baseline, _geo, _prop = _run_small(T_ignition=294.0)
    disabled, _geo, _prop = _run_small(
        T_ignition=294.0,
        diagnostic_disable_momentum=True,
    )
    baseline_momentum = max(np.max(s["momentum_source"]) for s in baseline["snapshots"])
    disabled_momentum = max(np.max(s["momentum_source"]) for s in disabled["snapshots"])
    assert baseline_momentum > 0.0
    assert disabled_momentum == pytest.approx(0.0)
