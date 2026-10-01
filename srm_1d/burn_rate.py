"""
burn_rate.py — Ma et al. (2020) Erosive Burning Model
=======================================================

PURPOSE:
    Computes the total burn rate (normal + erosive) at each cell in the
    motor grain port. This is the physics that distinguishes a 1D solver
    from a lumped-parameter approach: the burn rate varies along the
    grain because gas velocity (and thus convective heat transfer) varies.

THE MODEL (Ma et al. 2020, Eqs. 2-7):
    The total burn rate at any axial station is:

        r = r₀ + r_e

    where r₀ = a·Pⁿ is the normal (Saint-Robert) burn rate, and r_e is
    the erosive increment driven by convective heat transfer from the
    crossflow gas to the propellant surface.

    The erosive component comes from an energy balance at the surface:

        r_e = (T_flame - T_surface) / (ρ_p · C_ps · (T_surface - T_initial)) · h

    where h is the convective heat transfer coefficient, modified by
    transpiration (mass injection through the burning surface).

THE HEAT TRANSFER CHAIN:
    1. Haaland friction factor (Eq. 15 in Ma) — Darcy-Weisbach friction
       from the Haaland (1983) explicit correlation. Feeds into Gnielinski.

    2. Gnielinski Nusselt number (Eqs. 8-10 in Ma) — Convective heat
       transfer from the Gnielinski (2013) correlation, covering laminar,
       transition, and turbulent regimes. Includes entrance effects and
       a temperature-ratio correction.

    3. Transpiration correction (Eq. 16 in Ma) — The burning surface
       injects mass into the boundary layer, thickening it and reducing
       heat transfer. The correction is h/h₀ = β/(exp(β) - 1) where
       β = ρ_p · r · C_p / h₀.

    4. Bisection solver — The total burn rate r appears on both sides of
       the equation (through the transpiration correction), making it
       implicit. Fixed-point iteration diverges when the Jacobian > 1
       (which happens at moderate erosion). Bisection on the residual
       F(r) = r - r₀ - r_e(r) is monotonic and always converges.

ZERO ARBITRARY CONSTANTS:
    Every parameter in this model traces to either a named correlation
    (Haaland, Gnielinski) or a physical property (from CEA/RPA or
    propellant characterization). There are no tuning factors.

REFERENCES:
    Ma, Y. et al. (2020). "A New Erosive Burning Model of Solid
    Propellant Based on Heat Transfer." Int. J. Aerospace Eng., 2020.

    Haaland, S.E. (1983). "Simple and explicit formulas for the friction
    factor in turbulent pipe flow." J. Fluids Eng., 105(1), 89-90.

    Gnielinski, V. (2013). "On heat transfer in tubes." Int. J. Heat
    Mass Transfer, 63, 134-140.
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
# Step 1: Friction Factor
# ================================================================

@njit(cache=True, fastmath=True)
def haaland_friction(Re, roughness, D):
    """
    Darcy friction factor from the Haaland (1983) correlation.

    Computes the Darcy-Weisbach friction factor for pipe flow across
    the full Reynolds number range. This feeds into the Gnielinski
    Nusselt number correlation — the friction factor appears in the
    turbulent heat transfer formula.

    Parameters
    ----------
    Re : float
        Reynolds number based on hydraulic diameter [-].
    roughness : float
        Absolute surface roughness [m]. For propellant grain surfaces,
        typical values are 10-50 μm.
    D : float
        Hydraulic diameter [m].

    Returns
    -------
    f : float
        Darcy friction factor [-].

    Notes
    -----
    Three regimes:

    Laminar (Re < 2300):
        f = 64/Re  — exact Hagen-Poiseuille solution for fully
        developed laminar pipe flow.

    Turbulent (Re > 4000):
        Haaland's explicit approximation to the implicit Colebrook-White
        equation: 1/√f = -1.8·log₁₀((ε/D)^1.11/3.7^1.11 + 6.9/Re)
        Accurate to within 1.5% of Colebrook-White for all Re > 4000.

    Transition (2300 < Re < 4000):
        Linear blend between laminar and turbulent values. This is a
        simplification — the real transition involves intermittency —
        but it provides a smooth, continuous function that prevents
        numerical oscillations in the solver.

    References
    ----------
    Haaland, S.E. (1983). J. Fluids Eng., 105(1), 89-90.
    Ma et al. (2020), Eq. 15.
    """
    Re = max(Re, 1.0)
    eps_rel = roughness / max(D, 1e-10)

    f_lam = 64.0 / Re

    inv_sqrt_f = -1.8 * np.log10(eps_rel**1.11 / 3.7**1.11 + 6.9 / Re)
    inv_sqrt_f = max(inv_sqrt_f, 1.0)
    f_turb = 1.0 / (inv_sqrt_f * inv_sqrt_f)

    if Re < 2300.0:
        return f_lam
    elif Re > 4000.0:
        return f_turb
    else:
        gamma_t = (Re - 2300.0) / 1700.0
        return (1.0 - gamma_t) * f_lam + gamma_t * f_turb


# ================================================================
# Step 2: Nusselt Number (Heat Transfer)
# ================================================================

@njit(cache=True, fastmath=True)
def gnielinski_nusselt(Re, Pr, D, L, f, T_gas, T_surface, kappa):
    """
    Nusselt number from the Gnielinski (2013) correlation.

    Computes the convective heat transfer coefficient (as a Nusselt
    number) for internal pipe flow. The Nusselt number relates the
    heat transfer coefficient to the thermal conductivity:
        h₀ = Nu · k / D

    This is the heat transfer WITHOUT transpiration correction — the
    "bare" coefficient before accounting for mass injection through
    the burning surface.

    Parameters
    ----------
    Re : float
        Reynolds number based on hydraulic diameter [-].
    Pr : float
        Prandtl number [-]. For combustion gases, use the EFFECTIVE
        value (0.3-0.5), not the frozen value (0.4-0.7).
    D : float
        Hydraulic diameter [m].
    L : float
        Distance from the head end [m]. Appears in the entrance
        correction: heat transfer is enhanced in the developing
        region near the entrance.
    f : float
        Darcy friction factor from haaland_friction [-].
    T_gas : float
        Core gas temperature [K] (= flame temperature).
    T_surface : float
        Propellant surface temperature [K].
    kappa : float
        Temperature-ratio exponent for the gas-heating correction.
        Gnielinski recommends κ = 0.45 for gases being heated
        (the wall is cooler than the gas). This accounts for the
        variation of gas properties across the thermal boundary layer.

    Returns
    -------
    Nu : float
        Nusselt number [-]. Minimum value clamped to 3.66 (the
        exact solution for fully developed laminar flow with
        constant wall temperature).

    Notes
    -----
    Three regimes:

    Laminar (Re < 2300) — Ma Eq. 8:
        Composite of three asymptotic Nusselt contributions:
        Nu₁ = 3.66 (fully developed)
        Nu₂ = 1.615·(Re·Pr·D/L)^(1/3) (thermally developing)
        Nu₃ = f(Pr)·(Re·Pr·D/L)^(1/2) (simultaneously developing)
        Combined: Nu_lam = (Nu₁³ + 0.7³ + (Nu₂-0.7)³ + Nu₃³)^(1/3)

    Turbulent (Re > 4000) — Ma Eq. 9:
        Modified Petukhov-Kirillov with Gnielinski's correction:
        Nu_turb = (f/8)·(Re-1000)·Pr / (1 + 12.7·√(f/8)·(Pr^(2/3)-1))
                  · (1 + (D/L)^(2/3)) · (T_gas/T_surface)^κ

    Transition (2300 < Re < 4000) — Ma Eq. 10:
        Linear blend: Nu = (1-γ)·Nu_lam + γ·Nu_turb
        where γ = (Re - 2300) / 1700.

    References
    ----------
    Gnielinski, V. (2013). Int. J. Heat Mass Transfer, 63, 134-140.
    Ma et al. (2020), Eqs. 8-10.
    """
    Re = max(Re, 1.0)
    L = max(L, D)

    Re_Pr_D_L = Re * Pr * D / L

    # Laminar (Eq. 8)
    Nu1 = 3.66
    Nu2 = 1.615 * Re_Pr_D_L ** (1.0 / 3.0)
    Nu3 = (2.0 / (1.0 + 22.0 * Pr)) ** (1.0 / 6.0) * Re_Pr_D_L ** 0.5
    Nu2_term = max(Nu2 - 0.7, 0.0)
    Nu_lam = (Nu1**3 + 0.7**3 + Nu2_term**3 + Nu3**3) ** (1.0 / 3.0)

    # Turbulent (Eq. 9)
    if kappa > 0.0:
        K = (T_gas / max(T_surface, 1.0)) ** kappa
    else:
        K = 1.0
    entrance = 1.0 + (D / L) ** (2.0 / 3.0)
    f8 = f / 8.0
    num = f8 * (Re - 1000.0) * Pr
    den = 1.0 + 12.7 * np.sqrt(f8) * (Pr ** (2.0 / 3.0) - 1.0)
    Nu_turb = max(num / max(den, 1e-10), 0.0) * entrance * K

    # Transition blend (Eq. 10)
    gamma_t = min(max((Re - 2300.0) / 1700.0, 0.0), 1.0)
    Nu = (1.0 - gamma_t) * Nu_lam + gamma_t * Nu_turb
    return max(Nu, 3.66)


# ================================================================
# Step 3: Transpiration Correction
# ================================================================

@njit(cache=True, fastmath=True)
def transpiration_correction(beta):
    """
    Transpiration (blowing) correction for heat transfer.

    When propellant burns, mass is injected through the surface into
    the boundary layer. This thickens the boundary layer and reduces
    heat transfer compared to the non-blowing case. The correction is:

        h / h₀ = β / (exp(β) - 1)

    where β = ρ_p · r · C_p / h₀ is the blowing parameter.

    Parameters
    ----------
    beta : float
        Blowing parameter [-]. β > 0 means mass injection (burning
        surface), which reduces heat transfer. β = 0 gives h/h₀ = 1.
        Large β (intense burning) drives h/h₀ toward 0.

    Returns
    -------
    ratio : float
        h/h₀, the ratio of actual to non-blowing heat transfer
        coefficient. Clamped to [0, 1].

    Notes
    -----
    For small β (< 1e-6), the Taylor expansion β/(exp(β)-1) ≈ 1 - β/2
    is used to avoid 0/0.

    This correction creates the self-limiting nature of erosive burning:
    higher burn rate → more mass injection → lower heat transfer → less
    erosive burning. This negative feedback is why erosive burning
    reaches a finite enhancement rather than running away.

    References
    ----------
    Ma et al. (2020), Eq. 16.
    """
    if abs(beta) < 1e-6:
        return 1.0 - beta / 2.0
    beta_c = min(max(beta, -500.0), 500.0)
    result = beta_c / (np.exp(beta_c) - 1.0)
    return min(max(result, 0.0), 1.0)


# ================================================================
# Step 4: Single-Cell Burn Rate (Bisection Solver)
# ================================================================

@njit(cache=True, fastmath=True)
def select_tab_idx(P, tab_min_p, tab_max_p, n_tabs):
    """
    Hard-switchover tab lookup matching openMotor's
    Propellant.getCombustionProperties: strict containment first
    (min < P < max), else closest boundary fallback.
    """
    for k in range(n_tabs):
        if tab_min_p[k] < P < tab_max_p[k]:
            return k
    best = 0
    d_lo0 = abs(P - tab_min_p[0])
    d_hi0 = abs(P - tab_max_p[0])
    best_dist = d_lo0 if d_lo0 < d_hi0 else d_hi0
    for k in range(1, n_tabs):
        d_lo = abs(P - tab_min_p[k])
        d_hi = abs(P - tab_max_p[k])
        d = d_lo if d_lo < d_hi else d_hi
        if d < best_dist:
            best = k
            best_dist = d
    return best


@njit(cache=True, fastmath=True)
def saint_robert_from_tabs(P, tab_min_p, tab_max_p, tab_a, tab_n, n_tabs):
    """r₀ = a(P)·P^n(P) using tab-lookup."""
    P_pos = max(P, 0.0)
    if P_pos == 0.0:
        return 0.0
    k = select_tab_idx(P, tab_min_p, tab_max_p, n_tabs)
    return tab_a[k] * P_pos ** tab_n[k]


@njit(cache=True, fastmath=True)
def _burn_rate_closure_residual(r, r0, alpha, ebf, h0):
    """Residual of the implicit Ma burn-rate closure at candidate ``r``."""
    beta = alpha * r
    if beta < 1e-6:
        h = h0 * (1.0 - beta / 2.0)
    elif beta > 500.0:
        h = 0.0
    else:
        h = h0 * beta / (np.exp(beta) - 1.0)
    return r - r0 - ebf * h


@njit(cache=True, fastmath=True)
def burn_rate_cell(
    P, Re_local, D_hyd, x_from_head, roughness,
    Pr, k_thermal, Cp_gas, T_flame, T_surface,
    rho_p, Cps, T_initial,
    tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
    kappa,
):
    """
    Compute total burn rate at a single cell using bisection.

    This solves the implicit equation for the total burn rate r:

        r = r₀ + (T_f - T_s)/(ρ_p · C_ps · (T_s - T_i)) · h(r)

    where h(r) = h₀ · β/(exp(β)-1) depends on r through the blowing
    parameter β = ρ_p · r · C_p / h₀.

    WHY BISECTION (not fixed-point iteration):
        The fixed-point iteration r_{n+1} = f(r_n) has a Jacobian
        df/dr that depends on the transpiration correction slope. When
        erosion is moderate to strong (erosion ratio > ~1.5), the
        Jacobian exceeds 1 and fixed-point iteration diverges. This was
        discovered empirically during development and confirmed by
        computing df/dr analytically.

        Bisection on F(r) = r - r₀ - (erosive contribution at r) works
        because F is monotonically increasing in r: higher r means more
        transpiration, which reduces heat transfer, which reduces the
        erosive increment, so F grows. F(r₀) < 0 (pure normal rate has
        positive erosive deficit) and F(∞) > 0, guaranteeing a unique
        root in [r₀, ∞).

        30 iterations of bisection give precision of (r_hi - r_lo)/2³⁰
        ≈ 1e-9 relative to the initial bracket width, which is more
        than sufficient.

    Parameters
    ----------
    P : float
        Local pressure [Pa].
    Re_local : float
        Local Reynolds number [-].
    D_hyd : float
        Local hydraulic diameter [m].
    x_from_head : float
        Distance from the head end [m]. Used in entrance correction.
    roughness : float
        Surface roughness [m].
    Pr : float
        Prandtl number [-].
    k_thermal : float
        Gas thermal conductivity [W/(m·K)].
    Cp_gas : float
        Gas specific heat [J/(kg·K)].
    T_flame : float
        Gas flame temperature [K].
    T_surface : float
        Propellant surface temperature [K].
    rho_p : float
        Propellant density [kg/m³].
    Cps : float
        Propellant solid specific heat [J/(kg·K)].
    T_initial : float
        Initial propellant temperature [K].
    tab_min_p, tab_max_p, tab_a, tab_n : ndarray (n_tabs,)
        Per-tab Saint-Robert parameters and operating-pressure ranges.
    n_tabs : int
        Number of tabs.
    kappa : float
        Gnielinski temperature-ratio exponent [-].

    Returns
    -------
    r_total : float
        Total burn rate (normal + erosive) [m/s].
    r_erosive : float
        Erosive component only [m/s].

    References
    ----------
    Ma et al. (2020), Eqs. 2-7.
    """
    # Normal burn rate (Saint-Robert's law) with tab lookup
    r0 = saint_robert_from_tabs(P, tab_min_p, tab_max_p, tab_a, tab_n, n_tabs)

    if Re_local < 100.0:
        return r0, 0.0
    if D_hyd <= 0.0:
        return r0, 0.0

    # Friction factor
    f = haaland_friction(Re_local, roughness, D_hyd)

    # Nusselt number without transpiration
    Nu = gnielinski_nusselt(Re_local, Pr, D_hyd, x_from_head,
                             f, T_flame, T_surface, kappa)

    # Heat transfer coefficient without transpiration
    h0 = Nu * k_thermal / D_hyd

    if h0 < 1e-6:
        return r0, 0.0

    # Energy balance factor:
    #   ebf = (T_flame - T_surface) / (ρ_p · C_ps · (T_surface - T_initial))
    # This has units of [K / (kg/m³ · J/(kg·K) · K)] = [m³·s²/kg] ... no.
    # Actually: [K] / [kg/m³ · J/(kg·K) · K] = [K·m³·kg·K] / [kg·J·K]
    # = [m³/J] ... still wrong. Let me just check: ebf * h0 has units of
    # [m³/J?] * [W/m²] = ... The point is the algebra works out to [m/s].
    T_diff_gas = T_flame - T_surface
    T_diff_prop = T_surface - T_initial
    if T_diff_prop < 1.0:
        return r0, 0.0
    ebf = T_diff_gas / (rho_p * Cps * T_diff_prop)

    # Precompute α = ρ_p · C_p_gas / h0 (appears in β = α·r)
    alpha = rho_p * Cp_gas / h0

    # -----------------------------------------------------------
    # Bisection on F(r) = r - r₀ - ebf · h₀ · β/(exp(β) - 1)
    # F is monotonically increasing in r.
    # -----------------------------------------------------------

    # Preserve the historical 20x bracket when it contains the root. Expand
    # only when its upper residual is still negative.
    r_lo = r0
    r_hi = r0 * 20.0
    if r_hi <= r_lo:
        r_hi = r_lo + max(ebf * h0, 1e-12)

    F_lo = _burn_rate_closure_residual(r_lo, r0, alpha, ebf, h0)

    # If F_lo >= 0, no erosive burning possible (transpiration too strong)
    if F_lo >= 0.0:
        return r0, 0.0

    F_hi = _burn_rate_closure_residual(r_hi, r0, alpha, ebf, h0)
    for _ in range(60):
        if F_hi >= 0.0:
            break
        r_hi = r_hi * 2.0
        F_hi = _burn_rate_closure_residual(r_hi, r0, alpha, ebf, h0)

    # h/h0 is bounded by one, so this is a guaranteed physical upper bound.
    if F_hi < 0.0:
        r_hi = r0 + max(ebf * h0, 0.0)

    # Bisection: 30 iterations gives precision of r_range / 2^30 ≈ 1e-9
    for _ in range(30):
        r_mid = 0.5 * (r_lo + r_hi)
        F_mid = _burn_rate_closure_residual(r_mid, r0, alpha, ebf, h0)

        if F_mid < 0.0:
            r_lo = r_mid
        else:
            r_hi = r_mid

        if (r_hi - r_lo) < 1e-8:
            break

    r_total = 0.5 * (r_lo + r_hi)
    r_erosive = r_total - r0
    return r_total, max(r_erosive, 0.0)


# ================================================================
# Vectorized Wrapper
# ================================================================

@njit(cache=True, fastmath=True)
def compute_burn_rates(
    P, Re, D_hyd, x_centers, is_burning, roughness,
    Pr, k_thermal, Cp_gas, T_flame, T_surface,
    rho_p, Cps, T_initial,
    tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
    kappa, N,
):
    """
    Compute burn rates for all cells in the domain.

    Calls burn_rate_cell for each burning cell. Non-burning cells
    (burned out or gap cells) get zero burn rate.

    Parameters
    ----------
    P : ndarray (N,)
        Cell-center pressure [Pa].
    Re : ndarray (N,)
        Cell-center Reynolds number [-].
    D_hyd : ndarray (N,)
        Cell-center hydraulic diameter [m].
    x_centers : ndarray (N,)
        Cell-center axial positions [m].
    is_burning : ndarray (N,), bool
        True for cells with active combustion.
    roughness : float
        Surface roughness [m].
    Pr : float
        Prandtl number [-].
    k_thermal : float
        Gas thermal conductivity [W/(m·K)].
    Cp_gas : float
        Gas specific heat [J/(kg·K)].
    T_flame : float
        Flame temperature [K].
    T_surface : float
        Surface temperature [K].
    rho_p : float
        Propellant density [kg/m³].
    Cps : float
        Propellant specific heat [J/(kg·K)].
    T_initial : float
        Initial propellant temperature [K].
    tab_min_p, tab_max_p, tab_a, tab_n : ndarray (n_tabs,)
        Per-tab Saint-Robert parameters and operating-pressure ranges.
    n_tabs : int
        Number of tabs.
    kappa : float
        Gnielinski temperature-ratio exponent [-].
    N : int
        Number of cells.

    Returns
    -------
    r_total : ndarray (N,)
        Total burn rate at each cell [m/s].
    r_erosive : ndarray (N,)
        Erosive component at each cell [m/s].
    """
    r_total = np.zeros(N)
    r_erosive = np.zeros(N)
    for i in range(N):
        if is_burning[i]:
            r_total[i], r_erosive[i] = burn_rate_cell(
                P[i], Re[i], D_hyd[i], x_centers[i], roughness,
                Pr, k_thermal, Cp_gas, T_flame, T_surface,
                rho_p, Cps, T_initial,
                tab_min_p, tab_max_p, tab_a, tab_n, n_tabs,
                kappa,
            )
    return r_total, r_erosive
