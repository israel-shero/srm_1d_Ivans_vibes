"""
solver.py — 1D Compressible PISO Solver on a Staggered Grid
=============================================================

PURPOSE:
    Solves the 1D unsteady Navier-Stokes equations (mass, momentum,
    energy) for compressible flow using the PISO algorithm. This module
    contains ONLY the numerical building blocks — it knows nothing about
    propellant, burn rates, or grain geometry. Those are handled by
    burn_rate.py and simulation.py.

THE PISO ALGORITHM (per time step):
    1. MOMENTUM PREDICTOR — Solve for u* using the old pressure field.
       Semi-implicit: new velocity, old pressure.

    2. PRESSURE CORRECTION 1 — u* does not satisfy continuity. Derive
       a correction P' from the continuity residual, yielding a
       tridiagonal system solved by the Thomas algorithm (TDMA).
       Update: P = P_old + P',  u = u* - (dt/ρ) · ∇P'

    3. PRESSURE CORRECTION 2 — Repeat correction on the updated field.
       This second pass captures coupling terms that the first correction
       misses, giving better temporal accuracy without outer iterations.
       This is PISO's key advantage over SIMPLE, which instead repeats
       the full predictor-corrector until convergence.

    4. ENERGY EQUATION — Advect temperature with the corrected velocity
       field. Mass source cells inject gas at the flame temperature.

    5. EQUATION OF STATE — Update density: ρ = P / (R·T).

GRID LAYOUT:
    Staggered (MAC) variable arrangement:

        Cell centers: |  0  |  1  |  2  | ... | N-1 |
        Cell faces:  0     1     2     3    ...  N-1    N
                     ^                                  ^
                 head end                            nozzle
                 (wall)                              (choked)

    Scalars (P, ρ, T) live at cell centers — arrays of length N.
    Velocities (u) live at cell faces — array of length N+1.
    Face j sits between cells j-1 and j.

    The staggered arrangement prevents pressure-velocity decoupling
    (checkerboarding). Each face's pressure gradient uses its two
    immediately adjacent cells, so odd and even cells cannot decouple.
    This eliminates the need for Rhie-Chow interpolation.

BOUNDARY CONDITIONS:
    Head end (face 0): Solid wall — u[0] = 0.
    Nozzle end (face N): Signed open throat boundary. Flow is choked
    only when the chamber/ambient pressure ratio requires it; otherwise
    the boundary uses the subsonic isentropic relation and can reverse
    to ambient inflow during failed ignition or taildown.

TIME STEPPING:
    Adaptive CFL-based: dt = CFL · dx / max(|u| + a_sound)

NUMBA JIT:
    All functions are compiled with @njit for C-level speed. This
    requires all data to be passed as arrays/scalars — no Python
    objects inside JIT boundaries.

REFERENCES:
    Issa, R.I. (1986). "Solution of the Implicitly Discretised Fluid
    Flow Equations by Operator-Splitting." J. Comput. Phys., 62.

    Moukalled, F., Mangani, L., Darwish, M. (2016). The Finite Volume
    Method in Computational Fluid Dynamics. Springer. Ch. 15.
"""

import numpy as np

try:
    from numba import njit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        def wrapper(func):
            return func
        return wrapper


# ================================================================
# TDMA (Thomas Algorithm)
# ================================================================

@njit(cache=True, fastmath=True)
def thomas_solve(a, b, c, d, N):
    """
    Solve a tridiagonal system using the Thomas algorithm (TDMA).

    Solves: a[i]·x[i-1] + b[i]·x[i] + c[i]·x[i+1] = d[i]

    This is the workhorse of the pressure correction step. The pressure
    correction equation, when discretized on the staggered grid, produces
    a tridiagonal system because each cell's pressure correction couples
    only to its immediate east and west neighbors.

    Parameters
    ----------
    a : ndarray (N,)
        Sub-diagonal coefficients. a[0] is unused (no left neighbor
        for the first cell).
    b : ndarray (N,)
        Main diagonal coefficients.
    c : ndarray (N,)
        Super-diagonal coefficients. c[N-1] is unused.
    d : ndarray (N,)
        Right-hand side vector.
    N : int
        System size.

    Returns
    -------
    x : ndarray (N,)
        Solution vector.

    Notes
    -----
    TDMA is O(N) — it exploits the tridiagonal structure to avoid the
    O(N³) cost of general Gaussian elimination. For our 150-cell grid,
    this means ~300 operations instead of ~3.4 million.

    The algorithm has two phases:
    1. Forward sweep: eliminate the sub-diagonal by modifying the
       diagonal and RHS, working left to right.
    2. Back substitution: solve for x working right to left.

    Numerically stable as long as the system is diagonally dominant,
    which the pressure correction equation guarantees (the diagonal
    coefficient includes contributions from both faces plus the
    transient density term).
    """
    x = np.zeros(N)
    # Forward sweep
    c_prime = np.zeros(N)
    d_prime = np.zeros(N)
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]
    for i in range(1, N):
        m = a[i] / (b[i] - a[i] * c_prime[i - 1])
        c_prime[i] = c[i] / (b[i] - a[i] * c_prime[i - 1])
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / (b[i] - a[i] * c_prime[i - 1])
    # Back substitution
    x[N - 1] = d_prime[N - 1]
    for i in range(N - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]
    return x


# ================================================================
# PISO Time Step
# ================================================================

NOZZLE_STATE_CHOKED_IN = -2
NOZZLE_STATE_SUBSONIC_IN = -1
NOZZLE_STATE_BALANCED = 0
NOZZLE_STATE_SUBSONIC_OUT = 1
NOZZLE_STATE_CHOKED_OUT = 2

LIMIT_TEMPERATURE_FLOOR = 0
LIMIT_TEMPERATURE_CEILING = 1
LIMIT_PRESSURE_FLOOR = 2
N_LIMIT_DIAGNOSTICS = 3


@njit(cache=True, fastmath=True)
def _critical_pressure_ratio(gamma):
    return (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))


@njit(cache=True, fastmath=True)
def _critical_flow_function(gamma):
    return np.sqrt(gamma) * (2.0 / (gamma + 1.0)) ** (
        (gamma + 1.0) / (2.0 * (gamma - 1.0))
    )


@njit(cache=True, fastmath=True)
def _subsonic_throat_mdot_mag(P0, T0, Pb, A_throat, gamma, R_specific):
    """Magnitude of isentropic subsonic throat flow from reservoir to back pressure."""
    if (P0 <= 0.0 or T0 <= 0.0 or Pb <= 0.0 or Pb >= P0 or
            A_throat <= 0.0 or gamma <= 1.0 or R_specific <= 0.0):
        return 0.0
    pr = Pb / P0
    e1 = 2.0 / gamma
    e2 = (gamma + 1.0) / gamma
    f = pr ** e1 - pr ** e2
    if f <= 0.0:
        return 0.0
    coeff = 2.0 * gamma / ((gamma - 1.0) * R_specific * T0)
    return A_throat * P0 * np.sqrt(coeff * f)


@njit(cache=True, fastmath=True)
def _subsonic_outflow_dmdp(P_cell, T_cell, P_ambient, A_throat, gamma, R_specific):
    """Derivative d(mdot_out)/dP_cell for the subsonic outflow branch."""
    pr = P_ambient / P_cell
    e1 = 2.0 / gamma
    e2 = (gamma + 1.0) / gamma
    f = pr ** e1 - pr ** e2
    if f <= 1.0e-16:
        return 0.0
    coeff = 2.0 * gamma / ((gamma - 1.0) * R_specific * T_cell)
    root = np.sqrt(coeff * f)
    dfdpr = e1 * pr ** (e1 - 1.0) - e2 * pr ** (e2 - 1.0)
    dmdp = A_throat * (root - 0.5 * coeff * pr * dfdpr / root)
    if dmdp < 0.0:
        return 0.0
    return dmdp


@njit(cache=True, fastmath=True)
def _subsonic_inflow_dmdp(P_cell, T_ambient, P_ambient, A_throat, gamma, R_specific):
    """Derivative d(signed mdot)/dP_cell for the subsonic ambient-inflow branch."""
    pr = P_cell / P_ambient
    e1 = 2.0 / gamma
    e2 = (gamma + 1.0) / gamma
    f = pr ** e1 - pr ** e2
    if f <= 1.0e-16:
        return 0.0
    coeff = 2.0 * gamma / ((gamma - 1.0) * R_specific * T_ambient)
    root = np.sqrt(coeff * f)
    dfdpr = e1 * pr ** (e1 - 1.0) - e2 * pr ** (e2 - 1.0)
    dmdp = -0.5 * A_throat * coeff * dfdpr / root
    if dmdp < 0.0:
        return 0.0
    return dmdp


@njit(cache=True, fastmath=True)
def _nozzle_boundary_flow(
    P_cell, T_cell, A_throat, gamma, R_specific, P_ambient, T_ambient,
):
    """
    Signed quasi-steady isentropic throat boundary.

    Returns ``(mdot, dmdot_dP_cell, upstream_temperature, state_code)``.
    Positive ``mdot`` is chamber outflow. Negative ``mdot`` is ambient
    inflow through the throat.
    """
    if A_throat <= 0.0 or gamma <= 1.0 or R_specific <= 0.0:
        return 0.0, 0.0, T_cell, NOZZLE_STATE_BALANCED

    Pc = max(P_cell, 1.0)
    Tc = max(T_cell, 1.0)
    Pa = max(P_ambient, 1.0)
    Ta = max(T_ambient, 1.0)
    balance_tol = max(1.0e-6, 1.0e-10 * Pa)
    if abs(Pc - Pa) <= balance_tol:
        return 0.0, 0.0, Tc, NOZZLE_STATE_BALANCED

    gamma_fn = _critical_flow_function(gamma)
    crit_pr = _critical_pressure_ratio(gamma)

    if Pc > Pa:
        pr = Pa / Pc
        if pr <= crit_pr:
            coeff = A_throat * gamma_fn / np.sqrt(R_specific * Tc)
            return Pc * coeff, coeff, Tc, NOZZLE_STATE_CHOKED_OUT

        mdot = _subsonic_throat_mdot_mag(Pc, Tc, Pa, A_throat, gamma, R_specific)
        dmdp = _subsonic_outflow_dmdp(Pc, Tc, Pa, A_throat, gamma, R_specific)
        return mdot, dmdp, Tc, NOZZLE_STATE_SUBSONIC_OUT

    pr = Pc / Pa
    if pr <= crit_pr:
        mdot = Pa * A_throat * gamma_fn / np.sqrt(R_specific * Ta)
        return -mdot, 0.0, Ta, NOZZLE_STATE_CHOKED_IN

    mdot = _subsonic_throat_mdot_mag(Pa, Ta, Pc, A_throat, gamma, R_specific)
    dmdp = _subsonic_inflow_dmdp(Pc, Ta, Pa, A_throat, gamma, R_specific)
    return -mdot, dmdp, Ta, NOZZLE_STATE_SUBSONIC_IN


@njit(cache=True, fastmath=True)
def _piso_step_with_energy_diagnostics(
    rho, u, P, T, A_port, D_hyd,
    mass_source, thermal_source, momentum_source, f_darcy,
    dx, dt, gamma_arr, R_arr, Cp_arr, T_ceiling_arr,
    A_throat, P_ambient, T_ambient, N,
    port_mach_cap=0.0,
    limit_activation_counts=None,
    limit_duration_s=None,
    limit_abs_energy_j=None,
    limit_first_time_s=None,
    limit_last_time_s=None,
    pressure_floor_max_deficit_pa=None,
    mach_limit_activation_counts=None,
    mach_limit_duration_s=None,
    mach_limit_first_time_s=None,
    mach_limit_last_time_s=None,
    step_start_time_s=0.0,
    pressure_floor_pa=1.0e3,
):
    """
    One complete PISO time step on a staggered grid.

    This function advances the flow field by one time step. It takes the
    current state (ρ, u, P, T) and source terms (mass_source, f_darcy)
    and returns the updated state. It has no knowledge of where the mass
    source or friction come from — that is the simulation driver's job.

    v0.7.1 (Phase 3): the previously-scalar gas thermo (gamma, R, Cp,
    T_flame) is now per-cell. The nozzle boundary uses cell-N-1's
    mixture (γ, R); the EOS update uses cell-i's R; the pressure-correction
    transient term uses cell-i's R; and the energy equation advects
    sensible enthalpy ``Cp·T`` rather than scalar T so cells with
    different Cp conserve energy across faces. The T-clipping ceiling
    is also per-cell from ``T_ceiling_arr``.

    Parameters
    ----------
    rho : ndarray (N,)
        Cell-center density [kg/m³].
    u : ndarray (N+1,)
        Face velocities [m/s]. u[0] = head-end wall, u[N] = nozzle.
    P : ndarray (N,)
        Cell-center pressure [Pa].
    T : ndarray (N,)
        Cell-center temperature [K].
    A_port : ndarray (N,)
        Cell-center port cross-sectional area [m²].
    D_hyd : ndarray (N,)
        Cell-center hydraulic diameter [m].
    mass_source : ndarray (N,)
        Mass source per unit length [kg/(m·s)]. From propellant
        combustion and igniter.
    thermal_source : ndarray (N,)
        Enthalpy injection per unit length [W/m]. v0.7.1 unit change
        (Phase 3 step 1) from the previous mass-rate-times-T convention.
    momentum_source : ndarray (N+1,)
        Face-centered axial momentum source [N/m^3]. Positive values
        accelerate flow downstream in the momentum predictor.
    f_darcy : ndarray (N,)
        Darcy friction factor at cell centers [-].
    dx : float
        Cell width [m].
    dt : float
        Time step size [s].
    gamma_arr : ndarray (N,)
        Per-cell ratio of specific heats [-].
    R_arr : ndarray (N,)
        Per-cell specific gas constant [J/(kg·K)].
    Cp_arr : ndarray (N,)
        Per-cell specific heat at constant pressure [J/(kg·K)].
    T_ceiling_arr : ndarray (N,)
        Per-cell temperature ceiling [K] used for the T_raw clip.
        See ``simulation._compute_T_ceiling_arr``.
    A_throat : float
        Nozzle throat area [m²].
    P_ambient : float
        Ambient/back pressure [Pa].
    T_ambient : float
        Ambient reservoir temperature [K] used for reverse nozzle inflow.
    N : int
        Number of cells.

    Returns
    -------
    rho_new : ndarray (N,)
        Updated density.
    u_new : ndarray (N+1,)
        Updated face velocities.
    P_new : ndarray (N,)
        Updated pressure.
    T_new : ndarray (N,)
        Updated temperature.
    """
    # v0.7.1 Phase 3: nozzle boundary uses cell-N-1 mixture thermo
    gamma_bnd = gamma_arr[N - 1]
    R_bnd = R_arr[N - 1]
    # Face areas (interpolated from cell centers)
    A_face = np.zeros(N + 1)
    A_face[0] = A_port[0]
    for j in range(1, N):
        A_face[j] = 0.5 * (A_port[j - 1] + A_port[j])
    A_face[N] = A_port[N - 1]

    # -------------------------------------------------------
    # STEP 1: MOMENTUM PREDICTOR (at faces j = 1 .. N-1)
    # -------------------------------------------------------
    # For face j between cells j-1 and j:
    #   ρ_f · (u* - u^n)/dt = -(P[j] - P[j-1])/dx + convection + friction
    #
    # The pressure gradient (P[j] - P[j-1])/dx uses ADJACENT cells only.
    # This is the staggered grid's key property — no cell-skipping.

    u_star = np.zeros(N + 1)
    u_star[0] = 0.0  # Head-end wall BC

    for j in range(1, N):
        # Face j sits between cells j-1 and j
        rho_f = 0.5 * (rho[j - 1] + rho[j])
        A_f = A_face[j]

        if rho_f < 1e-10 or A_f < 1e-10:
            u_star[j] = 0.0
            continue

        # --- Pressure gradient (the clean staggered form) ---
        dPdx = (P[j] - P[j - 1]) / dx
        pres_force = -dPdx

        # --- Convective momentum flux ---
        # Momentum CV for face j spans from cell-center j-1 to cell-center j.
        # Fluxes cross at cell centers j-1 and j.
        # Velocity at cell center i ≈ 0.5*(u[i] + u[i+1])
        u_at_jm1 = 0.5 * (u[j - 1] + u[j])      # velocity at cell center j-1
        u_at_j = 0.5 * (u[j] + u[j + 1])          # velocity at cell center j

        # Mass flux at cell centers
        mdot_jm1 = rho[j - 1] * u_at_jm1 * A_port[j - 1]
        mdot_j = rho[j] * u_at_j * A_port[j]

        # Upwind velocity for momentum flux
        if mdot_jm1 >= 0:
            u_upwind_left = u[j - 1] if j > 1 else 0.0  # face j-1
        else:
            u_upwind_left = u[j]

        if mdot_j >= 0:
            u_upwind_right = u[j]
        else:
            u_upwind_right = u[j + 1] if j < N - 1 else u[j]

        conv = -(mdot_j * u_upwind_right - mdot_jm1 * u_upwind_left) / dx

        # --- Friction ---
        # Interpolate friction factor and hydraulic diameter to face
        f_f = 0.5 * (f_darcy[j - 1] + f_darcy[j])
        D_f = 0.5 * (D_hyd[j - 1] + D_hyd[j])
        friction = -f_f / (2.0 * D_f) * rho_f * abs(u[j]) * u[j]

        # --- Update ---
        RHS = pres_force + conv / A_f + friction + momentum_source[j]
        u_star[j] = u[j] + dt / rho_f * RHS

    # Nozzle face: extrapolate for now (corrected by pressure step)
    u_star[N] = u_star[N - 1]

    # -------------------------------------------------------
    # STEP 2: PRESSURE CORRECTION 1
    # -------------------------------------------------------
    # Continuity at cell i:
    #   (ρ_new - ρ_old)·A·dx/dt + (ṁ*[i+1] - ṁ*[i]) = S·dx
    #
    # With ρ_new = (P + P')/(R·T) and velocity correction
    #   u'[j] = -dt/ρ_f[j] · (P'[j] - P'[j-1]) / dx
    #
    # This gives the tridiagonal system for P'.

    # Momentum equation inverse diagonal coefficient at each face
    d_face = np.zeros(N + 1)
    for j in range(1, N):
        rho_f = 0.5 * (rho[j - 1] + rho[j])
        d_face[j] = dt / max(rho_f, 1e-6)

    a_sub = np.zeros(N)
    a_diag = np.zeros(N)
    a_sup = np.zeros(N)
    b_rhs = np.zeros(N)

    for i in range(N):
        RT_local = R_arr[i] * T[i]
        a_t = A_port[i] * dx / (RT_local * dt)  # Transient density term

        # West face coefficient (face i)
        if i > 0:
            coeff_w = A_face[i] * A_face[i] * d_face[i] / dx
        else:
            coeff_w = 0.0  # Wall: no flux correction

        # East face coefficient (face i+1)
        if i < N - 1:
            coeff_e = A_face[i + 1] * A_face[i + 1] * d_face[i + 1] / dx
        else:
            coeff_e = 0.0  # Nozzle: handled separately

        a_sub[i] = -coeff_w
        a_sup[i] = -coeff_e
        a_diag[i] = a_t + coeff_w + coeff_e

        # Nozzle BC: P' at last cell drives nozzle flow correction
        if i == N - 1:
            _mdot_bc, dmdp_bc, _T_bc, _state_bc = _nozzle_boundary_flow(
                P[i], T[i], A_throat, gamma_bnd, R_bnd, P_ambient, T_ambient
            )
            a_diag[i] += dmdp_bc

        # Continuity residual using u*
        if i > 0:
            rho_w = 0.5 * (rho[i - 1] + rho[i])
            mdot_star_w = rho_w * u_star[i] * A_face[i]
        else:
            mdot_star_w = 0.0

        if i < N - 1:
            rho_e = 0.5 * (rho[i] + rho[i + 1])
            mdot_star_e = rho_e * u_star[i + 1] * A_face[i + 1]
        else:
            mdot_star_e, _dmdp, _T_bc, _state = _nozzle_boundary_flow(
                P[i], T[i], A_throat, gamma_bnd, R_bnd, P_ambient, T_ambient
            )

        b_rhs[i] = mass_source[i] * dx - (mdot_star_e - mdot_star_w)

    P_prime = thomas_solve(a_sub, a_diag, a_sup, b_rhs, N)

    # Update pressure
    P_new = P + P_prime

    # Correct face velocities: u = u* - d_face · (P'[j] - P'[j-1])/dx
    u_new = np.copy(u_star)
    for j in range(1, N):
        u_new[j] = u_star[j] - d_face[j] * (P_prime[j] - P_prime[j - 1]) / dx
    u_new[0] = 0.0

    # -------------------------------------------------------
    # STEP 3: PRESSURE CORRECTION 2
    # -------------------------------------------------------
    rho_new_1 = np.zeros(N)
    for i in range(N):
        rho_new_1[i] = P_new[i] / (R_arr[i] * T[i])

    # Recompute d_face with updated density
    for j in range(1, N):
        rho_f = 0.5 * (rho_new_1[j - 1] + rho_new_1[j])
        d_face[j] = dt / max(rho_f, 1e-6)

    for i in range(N):
        RT_local = R_arr[i] * T[i]
        a_t = A_port[i] * dx / (RT_local * dt)

        if i > 0:
            coeff_w = A_face[i] * A_face[i] * d_face[i] / dx
        else:
            coeff_w = 0.0

        if i < N - 1:
            coeff_e = A_face[i + 1] * A_face[i + 1] * d_face[i + 1] / dx
        else:
            coeff_e = 0.0

        a_sub[i] = -coeff_w
        a_sup[i] = -coeff_e
        a_diag[i] = a_t + coeff_w + coeff_e
        if i == N - 1:
            _mdot_bc, dmdp_bc, _T_bc, _state_bc = _nozzle_boundary_flow(
                P_new[i], T[i], A_throat, gamma_bnd, R_bnd, P_ambient, T_ambient
            )
            a_diag[i] += dmdp_bc

        if i > 0:
            rho_w = 0.5 * (rho_new_1[i - 1] + rho_new_1[i])
            mdot_w2 = rho_w * u_new[i] * A_face[i]
        else:
            mdot_w2 = 0.0

        if i < N - 1:
            rho_e = 0.5 * (rho_new_1[i] + rho_new_1[i + 1])
            mdot_e2 = rho_e * u_new[i + 1] * A_face[i + 1]
        else:
            mdot_e2, _dmdp, _T_bc, _state = _nozzle_boundary_flow(
                P_new[i], T[i], A_throat, gamma_bnd, R_bnd, P_ambient, T_ambient
            )

        b_rhs[i] = mass_source[i] * dx - (mdot_e2 - mdot_w2)

    P_prime2 = thomas_solve(a_sub, a_diag, a_sup, b_rhs, N)

    P_new = P_new + P_prime2
    for j in range(1, N):
        u_new[j] = u_new[j] - d_face[j] * (P_prime2[j] - P_prime2[j - 1]) / dx
    u_new[0] = 0.0
    mdot_boundary, _dmdp, T_boundary, _state = _nozzle_boundary_flow(
        P_new[N - 1], T[N - 1], A_throat, gamma_bnd, R_bnd,
        P_ambient, T_ambient,
    )
    if mdot_boundary >= 0.0:
        rho_boundary = rho_new_1[N - 1]
    else:
        rho_boundary = P_ambient / (R_bnd * max(T_boundary, 1.0))
    u_new[N] = mdot_boundary / max(rho_boundary * A_face[N], 1.0e-12)

    # -------------------------------------------------------
    # STEP 3c: AERODYNAMIC-CHOKING VELOCITY LIMITER
    # -------------------------------------------------------
    # The pressure-based PISO has no choking limit on the interior port
    # (only the nozzle THROAT chokes via _nozzle_boundary_flow) and is not
    # shock-capturing, so during the violent ignition fill the cold/hot
    # contact velocity can blow up to grid-divergent supersonic values (max
    # fill Mach 3.4 / 4.8 / 12 for 50 / 100 / 200 cells; see
    # docs/v0_7_4/IGNITION_SPIKE_REOPENED.md §6). Physically a constant- or
    # diverging-area port cannot exceed ~Mach 1 under distributed mass
    # addition without a converging-diverging passage. When port_mach_cap > 0,
    # clamp the interior face velocities |u[j]| to port_mach_cap times the
    # local sound speed. No-op at the subsonic plateau (Mach ~0.36), so steady
    # operation is untouched; applied BEFORE the energy advection below so the
    # convective enthalpy flux is bounded too. The head wall (face 0) and the
    # nozzle face (N, set by the choked BC) are left as-is.
    if port_mach_cap > 0.0:
        for j in range(1, N):
            gamma_f = 0.5 * (gamma_arr[j - 1] + gamma_arr[j])
            R_f = 0.5 * (R_arr[j - 1] + R_arr[j])
            T_f = 0.5 * (T[j - 1] + T[j])
            a_f = (gamma_f * R_f * max(T_f, 1.0)) ** 0.5
            u_lim = port_mach_cap * a_f
            if u_new[j] > u_lim:
                u_new[j] = u_lim
                if mach_limit_activation_counts is not None:
                    face_index = j - 1
                    mach_limit_activation_counts[face_index] += 1
                    mach_limit_duration_s[face_index] += dt
                    if mach_limit_first_time_s[face_index] < 0.0:
                        mach_limit_first_time_s[face_index] = step_start_time_s
                    mach_limit_last_time_s[face_index] = step_start_time_s + dt
            elif u_new[j] < -u_lim:
                u_new[j] = -u_lim
                if mach_limit_activation_counts is not None:
                    face_index = j - 1
                    mach_limit_activation_counts[face_index] += 1
                    mach_limit_duration_s[face_index] += dt
                    if mach_limit_first_time_s[face_index] < 0.0:
                        mach_limit_first_time_s[face_index] = step_start_time_s
                    mach_limit_last_time_s[face_index] = step_start_time_s + dt

    # -------------------------------------------------------
    # STEP 3b: ENERGY EQUATION  (v0.7.1 Phase 3: sensible-enthalpy form)
    # -------------------------------------------------------
    # We advect h_i = Cp_arr[i] * T[i] (specific sensible enthalpy)
    # rather than T directly. Across a face between cells with different
    # Cp, the conserved scalar that matches the mass flux is h, not T —
    # advecting T directly loses energy proportional to ΔCp / Cp. See
    # docs/v0_7_1/DESIGN.md §7.
    #
    # Cell update (per mass-flow at faces, mass-conservative form):
    #   m_new * h_new = m_old * h_old
    #                 + dt * (flux_h_in - flux_h_out + thermal_source[i] * dx)
    # T_new = h_new / Cp_arr[i] (Cp stays at its start-of-step value within
    # the PISO step — species composition is updated by _advect_species
    # after PISO returns, then Cp_arr is refreshed for the next step.)
    gas_energy_before = 0.0
    gas_energy_after = 0.0
    convective_scalar_flux_power = 0.0
    nozzle_scalar_flux_power = 0.0
    thermal_source_power = 0.0
    clipping_correction_power = 0.0
    T_new = np.zeros(N)
    for i in range(N):
        Cp_i = Cp_arr[i]
        old_mass = rho[i] * A_port[i] * dx
        new_mass = max(rho_new_1[i] * A_port[i] * dx, 1e-16)
        h_old = Cp_i * T[i]
        gas_energy_before += old_mass * h_old

        # West face enthalpy flux (upwind using staggered face velocity)
        if i == 0:
            flux_h_w = 0.0
        else:
            rho_w = 0.5 * (rho_new_1[i - 1] + rho_new_1[i])
            mdot_w = rho_w * u_new[i] * A_face[i]
            if mdot_w >= 0:
                flux_h_w = mdot_w * Cp_arr[i - 1] * T[i - 1]
            else:
                flux_h_w = mdot_w * Cp_i * T[i]

        # East face enthalpy flux
        if i < N - 1:
            rho_e = 0.5 * (rho_new_1[i] + rho_new_1[i + 1])
            mdot_e = rho_e * u_new[i + 1] * A_face[i + 1]
            if mdot_e >= 0:
                flux_h_e = mdot_e * Cp_i * T[i]
            else:
                flux_h_e = mdot_e * Cp_arr[i + 1] * T[i + 1]
        else:
            mdot_e, _dmdp, T_upstream, _state = _nozzle_boundary_flow(
                P_new[i], T[i], A_throat, gamma_bnd, R_bnd,
                P_ambient, T_ambient,
            )
            # Outflow carries cell-N-1's mixture Cp; inflow uses cell-N-1's
            # Cp as a stand-in for ambient Cp (mdot_e is small in that
            # regime so the choice is negligible).
            flux_h_e = mdot_e * Cp_i * T_upstream
            nozzle_scalar_flux_power = flux_h_e

        convective_h_flux = flux_h_w - flux_h_e
        convective_scalar_flux_power += convective_h_flux
        # thermal_source is enthalpy injection [W/m] (Phase 3 step 1).
        thermal_source_power += thermal_source[i] * dx

        h_scalar = old_mass * h_old + dt * (
            convective_h_flux + thermal_source[i] * dx
        )
        h_new = h_scalar / new_mass
        T_raw = h_new / Cp_i if Cp_i > 0.0 else h_new
        T_floor = max(1.0, T_ambient)
        T_ceiling_i = T_ceiling_arr[i]
        T_clipped = T_raw
        limit_index = -1
        if T_clipped < T_floor:
            T_clipped = T_floor
            limit_index = LIMIT_TEMPERATURE_FLOOR
        if T_clipped > T_ceiling_i:
            T_clipped = T_ceiling_i
            limit_index = LIMIT_TEMPERATURE_CEILING
        correction_power_i = (T_clipped - T_raw) * new_mass * Cp_i / dt
        clipping_correction_power += correction_power_i
        if limit_index >= 0 and limit_activation_counts is not None:
            limit_activation_counts[limit_index, i] += 1
            limit_duration_s[limit_index, i] += dt
            limit_abs_energy_j[limit_index, i] += abs(correction_power_i) * dt
            if limit_first_time_s[limit_index, i] < 0.0:
                limit_first_time_s[limit_index, i] = step_start_time_s
            limit_last_time_s[limit_index, i] = step_start_time_s + dt
        T_new[i] = T_clipped
        gas_energy_after += new_mass * Cp_i * T_new[i]

    # -------------------------------------------------------
    # STEP 4: UPDATE DENSITY
    # -------------------------------------------------------
    for i in range(N):
        if P_new[i] < pressure_floor_pa:
            deficit = pressure_floor_pa - P_new[i]
            if (pressure_floor_max_deficit_pa is not None
                    and deficit > pressure_floor_max_deficit_pa[i]):
                pressure_floor_max_deficit_pa[i] = deficit
            if limit_activation_counts is not None:
                limit_activation_counts[LIMIT_PRESSURE_FLOOR, i] += 1
                limit_duration_s[LIMIT_PRESSURE_FLOOR, i] += dt
                if limit_first_time_s[LIMIT_PRESSURE_FLOOR, i] < 0.0:
                    limit_first_time_s[LIMIT_PRESSURE_FLOOR, i] = step_start_time_s
                limit_last_time_s[LIMIT_PRESSURE_FLOOR, i] = step_start_time_s + dt
            P_new[i] = pressure_floor_pa
    rho_new = np.zeros(N)
    for i in range(N):
        rho_new[i] = P_new[i] / (R_arr[i] * T_new[i])

    gas_sensible_dE_dt = (gas_energy_after - gas_energy_before) / dt
    energy_residual = (
        gas_sensible_dE_dt
        - convective_scalar_flux_power
        - thermal_source_power
        - clipping_correction_power
    )

    return (rho_new, u_new, P_new, T_new,
            gas_energy_before, gas_energy_after, gas_sensible_dE_dt,
            convective_scalar_flux_power, nozzle_scalar_flux_power,
            thermal_source_power, clipping_correction_power,
            energy_residual)


@njit(cache=True, fastmath=True)
def piso_step(
    rho, u, P, T, A_port, D_hyd,
    mass_source, thermal_source, momentum_source, f_darcy,
    dx, dt, gamma_arr, R_arr, Cp_arr, T_ceiling_arr,
    A_throat, P_ambient, T_ambient, N,
    port_mach_cap=0.0,
):
    """Public PISO step returning only the updated flow state.

    v0.7.1 (Phase 3): the scalar (gamma, R_specific, Cp_gas, T_flame)
    arguments are now per-cell arrays. See
    ``_piso_step_with_energy_diagnostics`` for the full signature.

    ``port_mach_cap`` (default 0.0 = disabled) bounds interior port face
    velocities to that Mach number; see the limiter block in
    ``_piso_step_with_energy_diagnostics``.
    """
    out = _piso_step_with_energy_diagnostics(
        rho, u, P, T, A_port, D_hyd,
        mass_source, thermal_source, momentum_source, f_darcy,
        dx, dt, gamma_arr, R_arr, Cp_arr, T_ceiling_arr,
        A_throat, P_ambient, T_ambient, N,
        port_mach_cap,
    )
    return out[0], out[1], out[2], out[3]


# ================================================================
# CFL Time Step
# ================================================================

@njit(cache=True, fastmath=True)
def compute_dt_cfl(u, a_sound, dx, N, cfl_target, dt_max):
    """
    Compute the adaptive time step from the CFL condition.

    For compressible flow, information propagates at the speed of sound
    relative to the flow. The CFL condition requires that the numerical
    domain of dependence contains the physical domain of dependence:

        CFL = (|u| + a) · dt / dx ≤ CFL_target

    Rearranging: dt = CFL_target · dx / max(|u| + a)

    Parameters
    ----------
    u : ndarray (N,)
        Velocity field (face or cell-center, depending on caller).
    a_sound : float
        Speed of sound [m/s] (max across all cells).
    dx : float
        Cell width [m].
    N : int
        Length of the velocity array.
    cfl_target : float
        Target CFL number. 0.5 is stable for explicit schemes;
        0.3 is conservative.
    dt_max : float
        Maximum allowed time step [s]. Prevents excessively large
        steps during low-velocity phases (e.g., before ignition).

    Returns
    -------
    dt : float
        Time step [s].
    """
    max_wave_speed = a_sound  # At minimum, the sound speed
    for i in range(N):
        ws = abs(u[i]) + a_sound
        if ws > max_wave_speed:
            max_wave_speed = ws
    dt = cfl_target * dx / max_wave_speed
    return min(dt, dt_max)


@njit(cache=True, fastmath=True)
def compute_dt_source_cap(
    rho, A_port, thermal_source, mass_source, N,
    Cp_arr, T_flame, T_ambient, source_cfl_factor, dt_floor,
):
    """Source-aware CFL: cap dt so per-step per-cell energy injection
    cannot change a cell's gas temperature by more than
    ``source_cfl_factor * (T_flame - T_ambient)``.

    Required to break a Hasegawa-A-class numerical resonance documented
    in srm_1d/docs/v0_7_0/audits/2026-05-20: at certain emissivities,
    radiation-driven ignition cascades fire 5-50 cells per millisecond,
    injecting energy faster than the wavespeed-CFL dt can resolve. The
    resulting Mach-front oscillation at the ignition front trips the
    numerical-collapse abort.

    v0.7.1 (Phase 3): ``thermal_source`` carries W/m (enthalpy injection
    per unit length) and ``Cp_arr`` is per-cell. The per-cell power is
    ``thermal_source[i] * dx`` (W); available thermal capacity is
    ``rho[i] * A_port[i] * dx * Cp_arr[i]`` (J/K):

        dt_cap[i] = source_cfl_factor * (T_flame - T_ambient)
                    * rho[i] * A_port[i] * Cp_arr[i] / |thermal_source[i]|

    Take the minimum across all cells with nonzero source. ``dt_floor``
    prevents an infinite spike from collapsing dt to zero.

    Parameters
    ----------
    rho : ndarray (N,)
        Cell-centered density [kg/m^3].
    A_port : ndarray (N,)
        Cell port (gas) cross-sectional area [m^2].
    thermal_source : ndarray (N,)
        Per-cell enthalpy source [W/m].
    mass_source : ndarray (N,)
        Per-cell mass source [kg/(s*m^3)]. Used only to skip cells where
        the source is purely advective (no propellant injection).
    Cp_arr : ndarray (N,)
        Per-cell gas mixture specific heat [J/(kg*K)].
    T_flame, T_ambient : float
        Bounds on the temperature range used for the per-step cap. Stay
        scalar for v0.7.1 — the cap range is a global control knob, not
        a per-cell mixture property.
    source_cfl_factor : float
        Fraction of the (T_flame - T_ambient) range allowed per step.
        0.10 (10%) is a standard CFL-like default.
    dt_floor : float
        Lower bound on the returned dt to keep the loop from stalling.
    """
    dT_max = source_cfl_factor * (T_flame - T_ambient)
    if dT_max <= 0.0:
        return 1.0e10
    dt_min = 1.0e10
    for i in range(N):
        thermal_abs = abs(thermal_source[i])
        if thermal_abs > 0.0 and mass_source[i] > 0.0:
            dt_i = dT_max * rho[i] * A_port[i] * Cp_arr[i] / thermal_abs
            if dt_i < dt_min:
                dt_min = dt_i
    if dt_min < dt_floor:
        dt_min = dt_floor
    return dt_min
