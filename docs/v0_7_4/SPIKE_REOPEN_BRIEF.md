# Ignition-Spike — RE-OPEN BRIEF (next-session kickoff)

**Read this first.** It is the curated entry point for a fresh session
deliberately re-opening the high-L/D ignition over-spike. It exists to make you
efficient (don't re-derive dead ends) and accurate (don't hallucinate) — *not*
to fence you in.

> **Evidence status update (2026-11-19):** This brief predates scoring protocol
> v0.1. Treat its empirical ratios and “no spike” interpretations as historical
> hypotheses, not registered validation results. The current registry blocks
> Chunc absolute timing pending an independent event anchor and calibration
> review; BALLSStick is qualitative-only, with pressure ratios ineligible until
> channel provenance, calibration, linearity, and timing are resolved. See
> [`validation_case_registry.json`](../v0_8_0/data/validation_case_registry.json).

---

## 0. Rules of engagement (from the user, carry these)

- **You are deliberately re-opening a "closed" investigation.** Every prior
  "ruled out" below is **information, not a fence.** Re-open any of it *with
  cause*. Question whether to re-run old studies or modify current methods.
  Do **not** treat hard rule-outs dogmatically.
- **Keep the physics-grounded Ma-2020 erosive model as much as *reasonably*
  possible** — a strong preference, not an absolute. A move that introduces a
  tuned figure (e.g. blending Ma toward Mukunda-Paul) is *allowed* but must be a
  deliberate, user-visible tradeoff against the no-tuning shipping property —
  never slipped in. See [[feedback_keep_ma_erosive_model]].
- **Treat the WHOLE solver core as suspect and possibly coupled.** The ignition
  model, sonic propagation / shock modeling, the **energy balance wired into
  PISO**, and the Goodman gate may interact. **Assume no historical
  code/architecture is infallible** — audit it, don't inherit it.
- **Anti-hallucination discipline:** Read the local PDFs in
  [`../references/`](../references/) (Ma 2020, Mukunda 1997, Beddini, Keller
  1966, Liñán-Williams, Lengellé, Kuo/Summerfield, …) — **do not WebFetch**
  physics you can read locally. **Verify before asserting** (grep/read the code
  or the PDF first — this session's habit caught a wrong docstring and a
  `KeyError` in the README). **Reproduce + measure**, don't speculate. Known DOI
  trap: `10.2514/3.49697` is **Lengellé 1975**, not King.

---

## 1. The goal (quantitative)

Find a **physically-justified closure** that reduces the high-L/D ignition
pressure spike so the sim matches **empirical static-fire data across the full
short → long L/D range** — not just one motor.

The original working hypothesis was that over-prediction increased with motor
L/D (Chunc and BALLSstick are both high-L/D; Hasegawa A's peak is the *late
progressive* peak, not the ignition transient). The current case audit does not
establish that relationship: Chunc's event time and pressure eligibility remain
open, and BALLSstick's pressure ratios are not eligible. Treat L/D and geometry
as hypotheses to test after case qualification, not as established signatures.

---

## 2. Historical solver diagnosis (empirical inference still pending)

Chunc (`motors/machbusterNew.ric`) was used as a diagnostic case. Its trace
has a small peak-to-plateau ratio within the recorded samples, but starts at
0.010 s and lacks an independently registered ignition zero; pressure
calibration and startup-window review also remain pending. It therefore does
not yet establish that the real motor had zero ignition spike or quantify the
sim-to-experiment spike difference. The following solver findings are
repository-reported analyses from the 2026-06-18 re-open (REOPENED §9), not
registered validation results:

- The spike is a **genuine Ma quasi-steady erosive feedback loop**
  (G → r_erosive → P → ρ → G) responding to a **modestly-elevated** transient
  mass flux — **sim aft-G is only 1.29× experimental**, with **no gas-dynamic
  inflation** (sim ρ|u| = 0.87× its own quasi-steady ṁ/A; **aft Mach 0.23** at
  the spike — the supersonic fill is *over* by then). The historical analysis
  compared this with a 346-point static-fire trace; that case's event-time and
  pressure eligibility remain pending under protocol v0.1.
- Of the hump, **~2/3 is real QS erosive physics** (Mukunda-Paul universal law
  gives η ≈ 2.08 at the spike's g) and **~1/3 is Ma's high-g over-prediction**
  (Ma η ≈ 2.76 vs MP 2.08; they agree within ~6% at the plateau).
- **It grid-CONVERGES to ~14 MPa** — the canonical 100-cell run *under-resolves*
  it (11.76 / 12.64 / 13.20 MPa for 50/100/200 cells), so the "true" model spike
  is ~1.7×, *worse* vs experiment. Finer grid → sharper contact front → tighter
  ignition window → bigger synchronized aft-G surge (**simultaneity
  sharpening**).

### 2.1 Historical case comparison—not a registered empirical target

The following values were exploratory calculations from June 2026. Current
protocol v0.1 eligibility is shown separately:

| Case | Historical trace observation | Current evidence status |
|---|---|---|
| **Chunc** (high L/D, **fore-tapering conical**) | 1.015× peak/plateau ratio in the 346-point trace | Starts at 0.010 s; independent event zero, calibration, and scoring window remain pending. Do not classify it as a verified “no ignition spike.” |
| **BALLSStick 2″ subscale** | Approx. 1.35–1.50× ratios calculated from the logged PSI column | Channel provenance/calibration/linearity and event time are unresolved. Qualitative observation only; ratios are not scoreable. |

Do not infer a high-L/D trend, a geometry effect, or a model success target
from this pair yet. Chunc's fore-tapering geometry and a possible relationship
between port geometry and aft-G distribution remain hypotheses. First complete
the case-specific timing, calibration, geometry, and window reviews in the
registry; do not tune a model to these exploratory ratios.

---

## 3. Ruled out **in isolation** (informed — reopen with cause)

Do not silently redo these; *do* reconsider them under coupling / whole-core
re-analysis. Full detail in REOPENED §7, §9.4.

| Direction | Status | Note |
|---|---|---|
| Igniter mass / topology / IC | ruled out | + gas over-production **FIXED** (`bfc2f3f`: ProPep impetus + `gas_mass_fraction`) |
| Ignition kernel (Goodman surface-T gate) | validated | Keller-Baer-Ryan 1966 validates surface-T + the 756 K value |
| Ignition sequencing (flame-front, Phase F) | bottoms at floor | (a real DeMar-bypass bug was found + fixed; front now works, still floors) |
| Burn-rate magnitude / Beddini derivable lag | ruled out | no parameter-free τ from Beddini |
| Numerical resolution | ruled out | every-step burn+geom, cfl 0.15 identical |
| Velocity / gas-dynamic G inflation | ruled out | **thread D, direct-measured** (§2 above); `port_mach_cap` decoupled from P_peak |
| Transient h_c relaxation | probe-falsified | only **postpones** the spike (t_peak tracks the window) |
| **Ma → Mukunda-Paul high-g recal** | **the one floor-breaking lever** | cuts spike ~1.33× without touching plateau — but **rejected on no-tuning product grounds.** *The user now softens this: it's back on the table as a deliberate tradeoff, not forbidden.* |

**Explicitly do NOT retry as a *lone* spike fix** (they postpone or are
decoupled): transient h_c establishment relaxation; `port_mach_cap`;
integrated-energy/endothermic ignition criterion (Keller validates the existing
gate; it's a ~35% fidelity item, not a spike fix). *These may still matter **in
combination** — see thread B.*

---

## 4. Genuinely-open threads (NOT closed — the real leads)

From REOPENED §8 (deliberately left open) + this session's whole-core framing:

- **(B) Lever COUPLING is largely untested.** Everything above was eliminated
  *in isolation* or in pairs. Build a small Cartesian sweep over {igniter
  gas-fraction, flame_front_velocity, tau_establishment, ignition-criterion,
  port_mach_cap} and look for **non-additive** interactions. The user's meta-
  caution: *something is still missing; treat "eliminated" as "eliminated
  alone."*
- **(C) Grid-climb / simultaneity sharpening.** The spike grid-converges to
  ~14 MPa. Likely a sharper contact front tightening the whole-grain ignition
  window → bigger synchronized aft-G surge. Is there a **physical** simultaneity
  limiter (finite flame-spread + finite per-cell ignition energy) that the
  grid-sharpening is missing? (Ties to the erosive feedback at the source.)
- **(E) Sonic propagation / shock-aware transient — the user's explicit
  suspect.** The pressure-based PISO is **not shock-capturing**; the ignition
  fill is a shock-tube-like blowdown; `port_mach_cap` is a crude clamp. Consider
  a proper aerodynamic-choking source term, artificial viscosity / flux limiting
  at contacts, or a **density-based (Riemann) sub-step** for the transient —
  making the fill physical *without* the clamp (which may surface coupling the
  clamp currently hides). This is also the right long-term flow-field fix. See
  `docs/core_loop_opt/` (the acoustic-CFL "Lever B" work touches the same seam).
- **(F, user) The energy balance wired into PISO — audit it as a system.** How
  combustion + igniter enthalpy is injected and advected (sensible-enthalpy
  `Cp·T`, upwind face fluxes, per-cell `T_ceiling` clip, the convective
  wall-loss + pyrogen-radiation debits) is a plausible coupled contributor to
  the fast snap-on and the G surge that has **not** been examined as a whole.
  Does the fill's energy bookkeeping over-heat the near-front gas → over-drive
  the ignition/erosive response?
- **The Beddini turbulent-BL-development erosive factor (the Ma-preserving
  lever).** `EROSIVE_CLOSURE_RESEARCH.md` option A: ramp the *erosive*
  enhancement by `f_dev(t_since_flow_arrival)` over the BL-development time
  (~L_e/u ~ ms) — the erosive lift needs a *developed* turbulent BL, which the
  fill hasn't established. **Keeps Ma as the steady limit** (satisfies "keep Ma
  reasonably"). The h_c probe (§9.2) scaled the *shared* Nu·k/D and only
  postponed — but **decoupling the erosive-development factor from the
  ignition-gate h_c was never cleanly isolated.** Legitimate reopen under the
  coupling lens; expect *partial* per the lit dive.

---

## 5. Read-first order + entry points

1. **This brief.**
2. `IGNITION_SPIKE_REOPENED.md` — **§8 (open threads)** + **§9 (latest state:
   G reconstruction, Keller validation, h_c falsification, igniter fix)**. The
   current source of truth.
3. `EROSIVE_CLOSURE_RESEARCH.md` — transient erosive-closure lit dive + the
   Beddini lever.
4. `IGNITION_SPIKE_CLOSEOUT.md` — earlier snapshot (already banners "re-opened").
5. `README.md` (this dir) + `TASKS.md` — the original v0.7.4 synthesis.
6. Memories: [[project_ignition_model_audit]], [[project_ignition_fidelity_candidates]],
   [[project_chunc_erosive_spike_diagnostic]], [[feedback_keep_ma_erosive_model]],
   [[feedback_no_unfounded_smoothing]].

**Reproduce + measure:**
- Repro: `examples/chunc_ignition_2x2.py`, `examples/ignition_spike_diagnostic.py`,
  `examples/run_chunc.py`.
- Diagnostics tool: `srm_1d/tools/ignition_diagnostics.py`.
- **The investigation probes: [`probes/`](probes/)** (index in
  `probes/README.md`) — the exact scripts that produced the findings above
  (G reconstruction, relaminarization, Mukunda compare, mach convergence, h_c
  decomp, igniter, cfl). Machine-absolute paths; re-run, question, modify.
- **Validation data — inventory + per-file load recipes in
  [`../../static_fire_data/README.md`](../../static_fire_data/README.md)**
  (the files do NOT share a column layout — `zerox_data.csv` has **reversed**
  columns; read the recipes or you will silently swap time/pressure).
  - ⚠️ **Chunc**: `static_fire_data/thomas_chunc_firing.csv` has 346 points,
    but both it and its two-column workbook start at 0.010 s / 8.466 MPa and
    contain no event marker, baseline, or pressure-calibration metadata. The
    1.015 peak/previously reported plateau ratio is descriptive only; do not
    call it a verified no-spike result. See the current case registry.
  - ⚠️ **BALLSStick**: no 3″ firing exists. The 2″ subscale raw DAQ is in
    `static_fire_data/raw/` but is **NOT usable as-is**: the record holds only
    ~0.21 s of pre-ignition data and **PSI is already ramping monotonically from
    the first sample while the logged volts stay flat** — so there is no
    quiescent baseline and the PSI column is *not* a linear map of the volts
    (provenance unknown). **Do not apply a constant offset.** Load-cell impulse
    is also 0.85× expected and time isn't zeroed; a justified 2″→3″ scaling is
    still needed. See `raw/ballsstick_subscale_raw.notes.md`.
  - ❌ **Low/short-L/D coverage is still absent** (no clean data yet). This is the
    remaining coverage gap — a fix must be shown **not to break** the
    near-zero-spike low-L/D case, so treat conclusions drawn only from high-L/D
    motors as provisional.

**Hard constraints:** keep Ma reasonably (§0); no unfounded smoothing/dispersion
fudges ([[feedback_no_unfounded_smoothing]]); respect roughness/kappa/k_solid
physical bounds.

---

## 6. One-line kickoff (paste to start the session)

> Read `docs/v0_7_4/SPIKE_REOPEN_BRIEF.md` and follow its rules of engagement.
> We are re-opening the high-L/D ignition over-spike to find a physically-
> justified closure that matches static-fire data across the short→long L/D
> range, keeping Ma as much as reasonably possible, treating the whole solver
> core (ignition, sonic/shock, PISO energy balance, Goodman) as suspect and
> possibly coupled, and treating prior rule-outs as information rather than
> fences. Start by reviewing current data eligibility; do not treat Chunc's
> historical 1.015 ratio as a no-spike validation result until event time,
> baseline/calibration, and windows are established. Then reproduce the solver
> result and pick a thread from §4.
