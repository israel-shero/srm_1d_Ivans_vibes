"""Tests for the PISO solver numerical building blocks."""
import numpy as np
import pytest
from srm_1d.solver import (
    NOZZLE_STATE_BALANCED,
    NOZZLE_STATE_CHOKED_OUT,
    NOZZLE_STATE_SUBSONIC_IN,
    NOZZLE_STATE_SUBSONIC_OUT,
    _piso_step_with_energy_diagnostics,
    thomas_solve,
    compute_dt_cfl,
    piso_step,
    _nozzle_boundary_flow,
)


class TestThomasSolve:
    def test_identity_system(self):
        """b=1, a=c=0: solution equals RHS."""
        N = 5
        x = thomas_solve(np.zeros(N), np.ones(N), np.zeros(N),
                         np.array([1., 2., 3., 4., 5.]), N)
        np.testing.assert_allclose(x, [1., 2., 3., 4., 5.])

    def test_known_tridiagonal(self):
        """Solve a known 4x4 system and check against numpy.linalg.solve."""
        N = 4
        a = np.array([0.0, -1.0, -1.0, -1.0])
        b = np.array([4.0, 4.0, 4.0, 4.0])
        c = np.array([-1.0, -1.0, -1.0, 0.0])
        d = np.array([1.0, 2.0, 3.0, 4.0])

        # Build full matrix for reference
        A = np.diag(b) + np.diag(a[1:], -1) + np.diag(c[:-1], 1)
        x_ref = np.linalg.solve(A, d)
        x_tdma = thomas_solve(a, b, c, d, N)
        np.testing.assert_allclose(x_tdma, x_ref, atol=1e-12)

    def test_large_diagonally_dominant(self):
        """100-element diagonally dominant system."""
        N = 100
        a = np.full(N, -1.0); a[0] = 0.0
        b = np.full(N, 3.0)
        c = np.full(N, -1.0); c[-1] = 0.0
        d = np.ones(N)

        A = np.diag(b) + np.diag(a[1:], -1) + np.diag(c[:-1], 1)
        x_ref = np.linalg.solve(A, d)
        x_tdma = thomas_solve(a, b, c, d, N)
        np.testing.assert_allclose(x_tdma, x_ref, atol=1e-10)


class TestCFL:
    def test_basic_cfl(self):
        """CFL condition should scale dt inversely with wave speed."""
        u = np.array([10.0, 20.0, 30.0])
        dt = compute_dt_cfl(u, 340.0, 0.01, 3, 0.5, 1.0)
        # max wave speed = 30 + 340 = 370
        expected = 0.5 * 0.01 / 370.0
        assert dt == pytest.approx(expected, rel=1e-6)

    def test_dt_max_cap(self):
        """dt should not exceed dt_max."""
        u = np.array([0.01])
        dt = compute_dt_cfl(u, 340.0, 0.01, 1, 0.5, 1e-6)
        assert dt == pytest.approx(1e-6)

    def test_zero_velocity(self):
        """With zero flow, dt is set by sound speed alone."""
        u = np.zeros(5)
        dt = compute_dt_cfl(u, 340.0, 0.01, 5, 0.5, 1.0)
        expected = 0.5 * 0.01 / 340.0
        assert dt == pytest.approx(expected, rel=1e-6)


class TestPisoSources:
    def _single_cell_source_step(self, A_port_value, mass_rate, source_temperature,
                                 dt=1.0e-4, diagnostics=False,
                                 limit_arrays=None):
        """Closed one-cell source step used by conservative energy tests.

        v0.7.1 (Phase 3): thermal_source carries W/m and PISO takes
        per-cell ``gamma_arr/R_arr/Cp_arr/T_ceiling_arr`` arrays.
        """
        N = 1
        gamma = 1.2
        R_specific = 300.0
        Cp_gas = 2000.0
        T_flame = 3000.0
        T_initial = 300.0
        rho_initial = 1.0
        P_initial = rho_initial * R_specific * T_initial

        rho = np.array([rho_initial])
        u = np.zeros(N + 1)
        P = np.array([P_initial])
        T = np.array([T_initial])
        A_port = np.array([A_port_value])
        D_hyd = np.array([0.035])
        mass_source = np.array([mass_rate])
        thermal_source = mass_source * source_temperature * Cp_gas
        momentum_source = np.zeros(N + 1)
        f_darcy = np.zeros(N)
        gamma_arr = np.full(N, gamma)
        R_arr = np.full(N, R_specific)
        Cp_arr = np.full(N, Cp_gas)
        T_ceiling_arr = np.full(N, T_flame * 1.01)

        step_func = _piso_step_with_energy_diagnostics if diagnostics else piso_step
        args = (
            rho, u, P, T, A_port, D_hyd,
            mass_source, thermal_source, momentum_source, f_darcy,
            0.01, dt, gamma_arr, R_arr, Cp_arr, T_ceiling_arr,
            0.0, P_initial, T_initial, N,
        )
        if diagnostics and limit_arrays is not None:
            return step_func(*args, 0.0, *limit_arrays)
        return step_func(*args)

    def test_single_cell_thermal_source_matches_conservative_temperature(self):
        """No-flow source update should conserve rho*T scalar content."""
        A_port = 1.0e-3
        dx = 0.01
        dt = 1.0e-4
        mass_rate = 0.05
        source_temperature = 1000.0
        old_mass = 1.0 * A_port * dx
        added_mass = mass_rate * dx * dt
        expected = (
            old_mass * 300.0 + added_mass * source_temperature
        ) / (old_mass + added_mass)

        _rho_new, _u_new, _P_new, T_new = self._single_cell_source_step(
            A_port, mass_rate, source_temperature, dt,
        )

        assert T_new[0] == pytest.approx(expected, rel=1.0e-8)

    def test_same_temperature_mass_source_preserves_temperature(self):
        """Adding gas at the cell temperature should not heat the cell."""
        _rho_new, _u_new, _P_new, T_new = self._single_cell_source_step(
            1.0e-3, 0.05, 300.0,
        )

        assert T_new[0] == pytest.approx(300.0, rel=1.0e-10, abs=1.0e-10)

    def test_same_per_length_source_heats_smaller_port_faster(self):
        """For fixed per-length source, lower gas mass heats faster."""
        small = self._single_cell_source_step(1.0e-3, 0.05, 1000.0)
        large = self._single_cell_source_step(4.0e-3, 0.05, 1000.0)

        assert small[3][0] > large[3][0]

    def test_single_cell_energy_diagnostics_close_for_source_update(self):
        """Diagnostic residual should close for an unclipped source update."""
        out = self._single_cell_source_step(
            1.0e-3, 0.05, 1000.0, diagnostics=True,
        )

        assert out[10] == pytest.approx(0.0)
        assert out[11] == pytest.approx(0.0, abs=1.0e-8)

    def test_temperature_ceiling_clip_is_localized(self):
        counts = np.zeros((3, 1), dtype=np.int64)
        durations = np.zeros((3, 1))
        energies = np.zeros((3, 1))
        dt = 1.0e-4

        self._single_cell_source_step(
            1.0e-3, 100.0, 10_000.0, dt,
            diagnostics=True,
            limit_arrays=(counts, durations, energies),
        )

        assert counts[1, 0] == 1
        assert durations[1, 0] == pytest.approx(dt)
        assert energies[1, 0] > 0.0
        assert counts[0, 0] == 0

    def test_thermal_source_controls_injection_temperature(self):
        """Same mass source with hotter thermal source should heat more.

        v0.7.1 (Phase 3): thermal_source carries W/m, and PISO takes
        per-cell ``gamma_arr/R_arr/Cp_arr/T_ceiling_arr`` arrays.
        """
        N = 3
        gamma = 1.2
        R_specific = 300.0
        Cp_gas = 2000.0
        T_flame = 3000.0
        rho = np.full(N, 1.0)
        u = np.zeros(N + 1)
        T = np.full(N, 300.0)
        P = rho * R_specific * T
        A_port = np.full(N, 1.0e-3)
        D_hyd = np.full(N, 0.035)
        mass_source = np.zeros(N)
        mass_source[0] = 0.05
        momentum_source = np.zeros(N + 1)
        f_darcy = np.zeros(N)

        gamma_arr = np.full(N, gamma)
        R_arr = np.full(N, R_specific)
        Cp_arr = np.full(N, Cp_gas)
        T_ceiling_arr = np.full(N, T_flame * 1.01)

        cold_source = mass_source * 500.0 * Cp_gas
        hot_source = mass_source * 2500.0 * Cp_gas

        cold = piso_step(
            rho.copy(), u.copy(), P.copy(), T.copy(), A_port, D_hyd,
            mass_source, cold_source, momentum_source, f_darcy,
            0.01, 1.0e-5, gamma_arr, R_arr, Cp_arr, T_ceiling_arr,
            1.0e-4, 101325.0, 300.0, N,
        )
        hot = piso_step(
            rho.copy(), u.copy(), P.copy(), T.copy(), A_port, D_hyd,
            mass_source, hot_source, momentum_source, f_darcy,
            0.01, 1.0e-5, gamma_arr, R_arr, Cp_arr, T_ceiling_arr,
            1.0e-4, 101325.0, 300.0, N,
        )

        assert hot[3][0] > cold[3][0]

    def test_momentum_source_accelerates_downstream_face(self):
        """Positive face momentum source should push flow downstream."""
        N = 3
        gamma = 1.2
        R_specific = 300.0
        Cp_gas = 2000.0
        T_flame = 3000.0
        rho = np.full(N, 1.0)
        u = np.zeros(N + 1)
        P = np.full(N, 101325.0)
        T = np.full(N, 300.0)
        A_port = np.full(N, 1.0e-3)
        D_hyd = np.full(N, 0.035)
        mass_source = np.zeros(N)
        thermal_source = np.zeros(N)
        f_darcy = np.zeros(N)

        gamma_arr = np.full(N, gamma)
        R_arr = np.full(N, R_specific)
        Cp_arr = np.full(N, Cp_gas)
        T_ceiling_arr = np.full(N, T_flame * 1.01)

        no_momentum = np.zeros(N + 1)
        with_momentum = np.zeros(N + 1)
        with_momentum[1] = 5.0e4

        base = piso_step(
            rho.copy(), u.copy(), P.copy(), T.copy(), A_port, D_hyd,
            mass_source, thermal_source, no_momentum, f_darcy,
            0.01, 1.0e-5, gamma_arr, R_arr, Cp_arr, T_ceiling_arr,
            1.0e-4, 101325.0, 300.0, N,
        )
        driven = piso_step(
            rho.copy(), u.copy(), P.copy(), T.copy(), A_port, D_hyd,
            mass_source, thermal_source, with_momentum, f_darcy,
            0.01, 1.0e-5, gamma_arr, R_arr, Cp_arr, T_ceiling_arr,
            1.0e-4, 101325.0, 300.0, N,
        )

        assert driven[1][1] > base[1][1]

    def test_nozzle_boundary_does_not_drain_ambient_chamber(self):
        """At ambient pressure with no sources, the nozzle should not create vacuum."""
        N = 3
        gamma = 1.2
        R_specific = 300.0
        Cp_gas = 2000.0
        T_flame = 3000.0
        P_ambient = 101325.0
        rho = np.full(N, 1.0)
        u = np.zeros(N + 1)
        P = np.full(N, P_ambient)
        T = np.full(N, 300.0)
        A_port = np.full(N, 1.0e-3)
        D_hyd = np.full(N, 0.035)
        mass_source = np.zeros(N)
        thermal_source = np.zeros(N)
        momentum_source = np.zeros(N + 1)
        f_darcy = np.zeros(N)

        gamma_arr = np.full(N, gamma)
        R_arr = np.full(N, R_specific)
        Cp_arr = np.full(N, Cp_gas)
        T_ceiling_arr = np.full(N, T_flame * 1.01)

        rho_new, u_new, P_new, T_new = piso_step(
            rho.copy(), u.copy(), P.copy(), T.copy(), A_port, D_hyd,
            mass_source, thermal_source, momentum_source, f_darcy,
            0.01, 1.0e-5, gamma_arr, R_arr, Cp_arr, T_ceiling_arr,
            1.0e-4, P_ambient, 300.0, N,
        )

        np.testing.assert_allclose(P_new, P, rtol=1.0e-10, atol=1.0e-6)
        np.testing.assert_allclose(T_new, T)
        assert np.all(np.isfinite(rho_new))
        assert np.allclose(u_new, 0.0)

    def test_subambient_chamber_draws_reverse_nozzle_inflow(self):
        """Below ambient, the open boundary should add reverse mass flow."""
        N = 3
        gamma = 1.2
        R_specific = 300.0
        Cp_gas = 2000.0
        T_flame = 3000.0
        P_ambient = 101325.0
        rho = np.full(N, 0.5)
        u = np.zeros(N + 1)
        P = np.full(N, 0.8 * P_ambient)
        T = np.full(N, 293.0)
        A_port = np.full(N, 1.0e-3)
        D_hyd = np.full(N, 0.035)
        mass_source = np.zeros(N)
        thermal_source = np.zeros(N)
        momentum_source = np.zeros(N + 1)
        f_darcy = np.zeros(N)

        gamma_arr = np.full(N, gamma)
        R_arr = np.full(N, R_specific)
        Cp_arr = np.full(N, Cp_gas)
        T_ceiling_arr = np.full(N, T_flame * 1.01)

        rho_new, u_new, P_new, _T_new = piso_step(
            rho.copy(), u.copy(), P.copy(), T.copy(), A_port, D_hyd,
            mass_source, thermal_source, momentum_source, f_darcy,
            0.01, 1.0e-5, gamma_arr, R_arr, Cp_arr, T_ceiling_arr,
            1.0e-4, P_ambient, 293.0, N,
        )

        assert u_new[-1] < 0.0
        assert P_new[-1] > P[-1]
        assert np.all(np.isfinite(rho_new))


class TestNozzleBoundary:
    def test_choked_outflow_matches_ideal_formula(self):
        gamma = 1.2
        R = 300.0
        T = 2500.0
        P = 2.0e6
        A = 1.0e-4
        mdot, dmdp, upstream_T, state = _nozzle_boundary_flow(
            P, T, A, gamma, R, 101325.0, 293.0,
        )
        gamma_fn = np.sqrt(gamma) * (2.0 / (gamma + 1.0)) ** (
            (gamma + 1.0) / (2.0 * (gamma - 1.0))
        )
        assert mdot == pytest.approx(P * A * gamma_fn / np.sqrt(R * T))
        assert dmdp == pytest.approx(A * gamma_fn / np.sqrt(R * T))
        assert upstream_T == pytest.approx(T)
        assert state == NOZZLE_STATE_CHOKED_OUT

    def test_subsonic_outflow_and_reverse_inflow_are_signed(self):
        gamma = 1.2
        R = 300.0
        A = 1.0e-4
        out = _nozzle_boundary_flow(110000.0, 300.0, A, gamma, R, 101325.0, 293.0)
        bal = _nozzle_boundary_flow(101325.0, 300.0, A, gamma, R, 101325.0, 293.0)
        inflow = _nozzle_boundary_flow(90000.0, 300.0, A, gamma, R, 101325.0, 293.0)

        assert out[0] > 0.0
        assert out[1] > 0.0
        assert out[3] == NOZZLE_STATE_SUBSONIC_OUT
        assert bal[0] == pytest.approx(0.0)
        assert bal[3] == NOZZLE_STATE_BALANCED
        assert inflow[0] < 0.0
        assert inflow[1] > 0.0
        assert inflow[2] == pytest.approx(293.0)
        assert inflow[3] == NOZZLE_STATE_SUBSONIC_IN
