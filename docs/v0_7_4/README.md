# v0.7.4 — Ignition-Transient Spike: Research Synthesis

> **NEXT SESSION START HERE → [`SPIKE_REOPEN_BRIEF.md`](SPIKE_REOPEN_BRIEF.md)**
> (curated kickoff for the deliberate re-open). Current state of the
> investigation is `IGNITION_SPIKE_REOPENED.md` §8–§9; the empirical probes are
> in [`probes/`](probes/). This README is the original 2026-05-30 synthesis
> (historical).

**Status** (2026-05-30, historical): research complete; since then Phase F/Z
shipped opt-in, the spike was closed-out then **re-opened** — see the brief.
**Evidence status update (2026-11-20):** the Chunc workbook/CSV contain only
time and pressure, start at 0.010 s / 8.466 MPa, and have no independent event
time or calibration/baseline metadata. Its historical “zero spike” reading is
not eligible as a validation conclusion under protocol v0.1; see the current
case registry in `docs/v0_8_0/data/validation_case_registry.json`.
**Goal**: assess the historical ~2× ignition pressure-spike over-prediction on
high-L:D motors (Chunc/machbusterNew is the primary candidate, pending source
qualification) without
numerical tuning — using physically-grounded, literature-backed
mechanisms.

This document is the durable record of three parallel literature
agents + a code-level pipeline analysis. The companion implementation
plan is [`TASKS.md`](TASKS.md). Primary PDFs live in
[`../references/`](../references/).

---

## 1. The problem, localized

Chunc (`machbusterNew.ric`) — high L:D fore-tapering conical, low-erosivity
high-performance propellant, MTV-pellet fore-end ignition. The historical
interpretation treated the trace as a zero-spike rise to a flat ~8.8 MPa
plateau, while simulations showed an early peak. The current evidence audit
finds no independent event-time marker, baseline, or pressure calibration in
the source workbook; the startup-spike classification and measured comparison
therefore remain unresolved under protocol v0.1. Preserve the earlier model
results as repository-reported analysis, not a validated error estimate.

### The pipeline, and where the spike is born

Per-step loop in `_run_time_loop` ([simulation.py:1150](../../simulation.py)):

```
dt (CFL + source cap)
 └─ STEP 1 geometry      advance_bore_regression(r_total[prev])      L1362
 └─ STEP 2 burn rates    compute_burn_rates → burn_rate_cell         L1389 / burn_rate.py:343
                         r = a·Pⁿ + r_erosive(Re)   ← MEMORYLESS, instantaneous
 └─ STEP 3 ignition+src  _goodman_ignition_sources_and_mass          L1532 / L910
                         per unignited cell: h_c·(T_gas−T_surf) → Goodman ODE → T_surf
                         gate: T_surf > T_ignition  ← INSTANT STEP    solid_thermal.py:118
                         per burning cell: mass_source = ρ_p·r·C_burn  L1114
 └─ STEP 4 PISO          mass_source → P                              L1659
 └─ STEP 5 post-PISO     → Re (feeds STEP 2 next step)                L1669
```

The loop `r → mass → P → Re → r` is **positive feedback with no state /
memory anywhere in the gas-coupling path** (only geometry `regress` and
the pre-ignition Goodman `delta` are stateful). Two independent
structural roots compound multiplicatively:

- **Root A — ignition simultaneity.** The gate `T_surf > T_ignition`
  ([solid_thermal.py:118](../../solid_thermal.py)) is a hard step, and
  each cell's Goodman solver is **axially uncoupled in the solid**. Once
  PISO fills the bore with hot gas (1–2 ms, sonic traverse), every cell
  sees similar `T_gas` and `h_c` and crosses threshold within a tight
  window → near-simultaneous lighting. Real motors light as a **front**
  propagating over tens of ms.
- **Root B — memoryless burn-rate response.** `burn_rate_cell`
  ([burn_rate.py:343](../../burn_rate.py)) is a pure algebraic function
  of *instantaneous* `(P, Re)`. The erosive term is a steady convective
  energy balance with **no thermal-relaxation lag**. The instant a cell
  ignites and `G = ṁ/A_port` jumps, full augmentation fires that step.

Compounded: simultaneous ignition (A) → simultaneous mass injection →
`G` rises bore-wide at once → instantaneous erosive augmentation
everywhere (B) → mass-production surge outruns the throat → spike.

`tau_establishment` ([simulation.py:1103](../../simulation.py)) is the
only existing damper — a *linear* ramp on `r_total` with a *fit*
timescale. It "didn't fit / needed per-motor tuning" historically and is
the numerical-tuning trap to avoid (`feedback_no_unfounded_smoothing`).

---

## 2. Three research threads (agents, 2026-05-30)

### Thread 1 — Z-N dynamic burn rate (top-down)
- Memoryless `r(P, Re)` is valid only in the quasi-steady limit; it
  breaks during rapid transients. The missing physics is **condensed-phase
  thermal inertia**, relaxation time **τ ≈ α_s/r²** with an O(1) prefactor.
- Numeric: α_s ≈ 1.5×10⁻⁷ m²/s, r = 3–8 mm/s → **τ ≈ 1.5–11 ms** — the
  same order as the observed artifact, so the lag is first-order
  non-negligible. **τ ∝ 1/r²** means smoothing is strongest at low burn
  rate (early ignition) and self-attenuates at the high-r plateau —
  it will *not* smear steady-state.
- **Confirmed directly by Lengellé (1975)**, line 642: "for v_b = 0.5
  cm/s τ_p is about 4×10⁻³ sec" → τ ≈ 4 ms at r = 5 mm/s. ✔
- Validation corpus is **depressurization/extinction, L\*-instability,
  T-burner response** — *not* ignition-spike-specific. Honest framing:
  we transfer a mechanism validated for one transient sign to the
  ignition transient by symmetry of the underlying condensed-phase
  physics. (Strongest anchor: Greatrix 2008; experimental reality of
  the lag: Strand 1974; extinction agreement: De Luca vol. 143 / ADA143573.)

### Thread 2 — Ma / SPINBALL erosive-burning literature
- **Ma 2020's burn rate is steady-state/memoryless** (our
  [burn_rate.py](../../burn_rate.py) transcribes it faithfully). The
  paper's "transient scheme" refers to the surrounding 1-D CFD fill, not
  the burn-rate closure. Nothing to borrow.
- **SPINBALL** (Cavallini 2009, DiGiacinto 2008) models the *full*
  ignition phenomenology — igniter source terms, **flame spreading**,
  induction interval, throat-seal rupture, chamber filling — but its
  burn rate is **quasi-steady Lenoir-Robillard**, and it *explicitly
  dismisses* Z-N as negligible for VEGA-scale motors (long fill, low
  bore mass flux). Key insight: **τ ≈ α/r² is fixed by propellant
  properties, not motor scale** — Cavallini's dismissal is correct for
  VEGA (τ_rise ≫ τ_ZN) and *wrong for our high-mass-flux regime*
  (τ_rise ≈ τ_ZN). SPINBALL's lesson is its *bottom-up* ignition
  architecture, not its burn-rate model.
- **Ma 2020 remains our erosive closure — do NOT replace it.** The
  Lenoir-Robillard, King (1993), and Lengellé (1975) erosive
  correlations are *inferior* to Ma 2020 and must not be implemented
  over it. They are cited here only for the *concept* that an erosive /
  velocity-coupled burn rate also exhibits a dynamic relaxation lag
  (Lengellé 1975, "Model Describing the Erosive Combustion and Velocity
  Response of Composite Propellants," AIAA J 13(3):315–322, DOI
  10.2514/3.49697 = [`lengelle1975.pdf`](../references/lengelle1975.pdf)).
  Z-N lags Ma's output; it does not change the erosive model.

### Thread 3 — SRM ignition-transient simulation architecture (bottom-up)
- The artifact maps onto the canonical **three-interval starting
  transient** (induction → **flame spreading** → chamber filling). Our
  model collapses the flame-spreading interval to ~0.
- Established transient codes (SPP's IGT module, SPINBALL, the foundational
  **Peretz-Kuo-Caveny-Summerfield 1973**, [`19740005393.pdf`](../references/19740005393.pdf))
  do **not** let bulk gas temperature trip ignition everywhere. They
  compute a **flame-spread front** by **successive heating-to-ignition**
  (Peretz: "the flame spreading rate is calculated by local successive
  [ignition]"), driven by *local convective heat flux* + a surface
  autoignition temperature, with the front velocity (~1–10 m/s, up to
  ~200 m/s measured) **decoupled from the acoustic/fill speed** (~300–1000
  m/s). A model that ignites on hot-gas arrival over-speeds ignition by
  ~10–100× — fully consistent with a ~2× peak.
- Unnikrishnan (2001): "altered variation of the flame spread rate will
  alter the starting transient as well as the ignition peak" — flame-spread
  timing is *the* primary lever on the spike. Simplified all-surface
  ignition is documented to over-predict the peak.
- **Verdict: established codes fix simultaneity BOTTOM-UP (flame-spread
  front). Z-N is a complementary top-down refinement of per-cell rate,
  never the simultaneity cure.**

---

## 3. The tension, and its resolution

The two threads point at the two different roots:

| | Root it attacks | Mechanism | Lit anchor |
|---|---|---|---|
| **Bottom-up (flame-spread front)** | A — *spatial* simultaneity | gate ignition to a front advancing by successive heating-to-ignition | Peretz 1973; SPP IGT; SPINBALL; Unnikrishnan 2001 |
| **Top-down (Z-N relaxation)** | B — *temporal* rate lag | relax `r` toward `r_qs` over τ = κ·α_s/r² | Greatrix 2008; Lengellé 1975; Strand 1974 |

They are **complementary, not competing**. Fixing the rate (B) without
fixing the sequencing (A) leaves the dominant error — simultaneous
mass-addition — intact, which is consistent with the prior finding that
τ-ramps "didn't fit." Established-code consensus: **get the front right
first.**

**Decision (user, 2026-05-30): implement the bottom-up flame-front fix
first; then Z-N; evaluate each independently and then combined.** The
plan in [`TASKS.md`](TASKS.md) follows this ordering.

---

## 4. The Greatrix nuance (scope-critical)

Reading Greatrix 2008 verbatim ([`greatrix2008.pdf`](../references/greatrix2008.pdf))
corrects a conflation in the agent summaries. Greatrix's actual model is
**two layers**:

1. **Physics (per cell):** solve 1-D solid-phase heat conduction
   `k_s ∂²T/∂x² = ρ_s C_s ∂T/∂t` (Eq. 9/10, RK4, Δt ≈ 1×10⁻⁷ s,
   Fourier-limited Δx into the solid), integrate the temperature
   distribution, and back out the *unconstrained* instantaneous rate
   `r*_b = r_b,qs − [1/(T_s−T_i−ΔH_s/C_s)] · ∂/∂t ∫ΔT dx` (Eq. 6).
2. **Numerical limiter:** `dr_b/dt = K_b (r*_b − r_b)` (Eq. 14), where
   **K_b is a tuned damping coefficient** ("set below a maximum
   permissible value for a nondivergent solution, and adjusted further
   downward to match combustion response behavior"), calibrated to
   T-burner data. His examples use K_b = 6,700–170,000 s⁻¹ → 1/K_b ≈
   **6–150 µs** — a *numerical* timescale, NOT the physical relaxation.

The **physical** relaxation time is the thermal-wave penetration/diffusion
time τ ≈ α_s/r² (~ms; Lengellé Eq. analog, Greatrix Eq. 19
`T = T_i + (T_s−T_i)·exp(−r_b|x|/α_s)`). It emerges naturally from
layer 1's conduction solve; K_b is just a stabilizer on top.

**Implication for srm_1d:** the full layer-1 solid mesh (an inner
spatial dimension per cell at Δt ≈ 1e-7 s) is intractable in our Numba
time loop. We adopt the **lumped reduction** — one ODE per cell,
`dr/dt = (r_qs − r)/τ`, τ = κ·α_s/r², κ ≈ 1 — which collapses Greatrix's
two layers into the physical relaxation directly. This is parameter-free
at κ = 1 (uses existing `k_solid`, `ρ_p`, `Cps`), defensible as the
lumped limit of Eq. 14 with τ set by the penetration physics rather than
a fitted K_b. We cite Greatrix for the lag-ODE *form* and Lengellé / Eq. 19
for the *timescale*; we explicitly do **not** reproduce his full
conduction solve. (Greatrix 2011, [`greatrix2011.pdf`](../references/greatrix2011.pdf),
is the open-access sibling; MDPI Aerospace 2023,
[`aerospace-10-00767-v2.pdf`](../references/aerospace-10-00767-v2.pdf),
gives closed-form Z-N closures if a non-ODE form is later preferred.)

---

## 5. Reference index (local PDFs in `../references/`)

| File | Paper | Role |
|---|---|---|
| `greatrix2008.pdf` | Greatrix 2008, IJAE 826070 | Z-N lag-ODE form (Eq. 14), the two-layer architecture, K_b nuance |
| `greatrix2011.pdf` | Greatrix, Energies 4(1):90 | open-access sibling model |
| `lengelle1975.pdf` | Lengellé 1975, velocity response | τ ≈ α/r² ≈ 4 ms @ 5 mm/s confirmation; velocity-coupled (erosive) response |
| `strand1974.pdf` | Strand microwave-Doppler | experimental transient-r deviates from r_qs on ms scale |
| `aerospace-10-00767-v2.pdf` | MDPI Aerospace 2023 | closed-form Z-N closures |
| `ADA143573.pdf` | (De Luca-adjacent report) | extinction/transient validation surrogate |
| `19740005393.pdf` | Peretz-Kuo-Caveny-Summerfield 1973, AIAA J 11(12) | foundational bottom-up flame-spread (successive heating-to-ignition) |
| `cavallini2009.pdf`, `digiacinto2008.pdf` | SPINBALL | bottom-up ignition architecture; quasi-steady burn rate |
| `ma2020...Heat.pdf` | Ma 2020 erosive | our steady erosive closure (no transient term) |
| `hasegawa2006.pdf` | Hasegawa 2006 | motor geometry (Table 3) + propellant data (Table 1) |

**No missing references.** Earlier drafts listed a "King velocity
response" paper as wanted — that was an agent misattribution. The title
"...Erosive Combustion and Velocity Response of Composite Propellants"
and DOI **10.2514/3.49697** are **Lengellé (1975)**, already present as
`lengelle1975.pdf`. King's separate erosive work is the 1993 *J.
Propulsion & Power* 9(6):785–805 review (a survey, not a model we need).
Ma 2020 (`...Heat.pdf`) is the erosive model of record and is in the
references folder.

---

## 6. Open risks carried into implementation

- **Re-calibration**: existing knobs (roughness 37 µm, k_solid, kappa)
  were tuned against the *spiky* baseline. Both fixes touch the spike →
  expect a cross-motor re-LHS (frozen transport per the v0.7.3.3
  frozen-wins finding) on the fired set: Hasegawa A, Zerox, Chunc,
  BALLSstick.
- **Flame-front stall**: a hard front gate must not deadlock if the
  frontmost cell never reaches T_ignition (need an ignition-source floor).
- **Z-N plateau preservation**: τ ∝ 1/r² should self-attenuate at the
  plateau; verify it does not depress steady-state P.
- **Numba cache**: delete `srm_1d/__pycache__/` after every @njit edit
  (gotcha #1) or stale compiled code masks the change.
