# Ignition-Spike Investigation — RE-OPENED (2026-06-16/17)

> **Next session: start with [`SPIKE_REOPEN_BRIEF.md`](SPIKE_REOPEN_BRIEF.md)**
> (curated kickoff + rules of engagement). This doc is the full record; §8 (open
> threads) and §9 (latest state) are the live parts. Empirical probes:
> [`probes/`](probes/).

**This document updates `IGNITION_SPIKE_CLOSEOUT.md`.** The user re-opened the
question with the explicit instruction to *suppress the closeout's prior
assumptions* and test whether the spike is a structural solver artifact rather
than fundamental burn-rate physics. The verdict after a full multi-angle
re-investigation:

> **Evidence status update (2026-11-19):** The empirical comparisons below
> predate scoring protocol v0.1. The registry does not currently qualify
> Chunc's absolute timing or pressure amplitude, and it permits only qualitative
> inspection of BALLSStick; the latter's pressure ratios are not eligible.
> Historical “no spike” and “same residual” phrases below are not current
> validation conclusions. See
> [`validation_case_registry.json`](../v0_8_0/data/validation_case_registry.json).

> **STATUS: OPEN — deliberately left open for a new session (2026-06-17, see §8).**
> The re-opening found and fixed two genuine solver issues — the Phase F flame-front
> DeMar-cascade bug, and a grid-divergent PISO velocity artifact in the ignition fill
> — but **neither drives the spike.** The spike grid-CONVERGES to a real value
> (~14 MPa) and reads as genuine quasi-steady erosive physics (Mukunda-Paul gives
> η≈2.1 at the spike condition). Every non-tuning lever tested IN ISOLATION fails to
> close it. BUT lever COUPLING is largely untested and the ms-scale ignition snap-on
> is still physically suspect (§8), so the closeout's erosive verdict holds only for
> what was tested — it is not the last word.
>
> **Update 2026-06-18 (§9):** thread D is RESOLVED — the spike is the erosive feedback,
> not a G-inflation artifact (sim G only 1.29× experiment; no gas-dynamic inflation;
> direct mass-flux measurement). The ignition criterion is Keller-1966-validated and
> the transient-h_c lever is probe-FALSIFIED (it only postpones the spike). A separate
> igniter gas-generation over-production was found and FIXED (commit bfc2f3f, ProPep
> impetus + condensed-phase split). The one surviving lever for tamping the erosive
> spike (Ma-vs-Mukunda-Paul high-g model fidelity, §9.5) was **REJECTED by the user
> (2026-06-18): Ma's no-tuning property is a SHIPPING requirement** (it lets srm_1d ship
> without users quantifying their own erosive burning), and any high-g correction adds a
> fitted figure. **→ The high-L/D ignition-transient erosive over-spike is now an
> ACCEPTED, DOCUMENTED LIMITATION — there is no remaining non-dogma lever, and the spike
> thread is effectively CLOSED on a product decision, not a physics gap.**
> Simulations of BALLSstick (bpnv, high-L/D) still produce a large early peak
> post-fix. Because its experimental pressure channel is uncalibrated, this is a
> predicted model behavior, not an empirically measured residual.

Two genuine solver issues were found and fixed (Phase F DeMar-cascade flame-front
bypass; PISO fill velocity artifact, `port_mach_cap`). Two physically-grounded
fidelity upgrades were recorded for future implementation. The diagnostic motor
throughout is Chunc / `machbusterNew` (high-L/D, canonical knobs: bpnv pyrogen,
Sutton sizing, roughness 32 µm, kappa 0.44, T_ignition 756 K).

---

## 1. Baseline & decomposition

| quantity | value |
|---|---|
| P_peak (baseline) | 12.65 MPa @ 11.6 ms |
| plateau (0.4–0.9 s) | 8.23 MPa |
| spike ratio | **1.54×** |
| experimental trace (Chunc) | starts at 0.010 s and 8.466 MPa; independent event zero, pressure calibration, and startup window pending—no ignition-spike classification yet |
| erosive-disabled | early spike **gone** (peak at end) → the spike is 100% erosive |

The plateau and decay already match experiment well; the entire model error is
the 0–25 ms transient.

---

## 2. The sharp early spike is structural (the supersonic-fill ignition gate)

**Ignition is near-simultaneous: all 100 grain cells light in 0.61–1.47 ms**
(0.86 ms spread). A fill-resolved probe shows why — a hot (~3100 K) gas-dynamic
contact front sweeps the bore at **~800 m/s** and the convective Goodman gate
ignites each cell the instant the hot gas arrives.

Tracing one mid-grain cell (x = 435 mm) through its ignition:

```
 t_ms  Mach   u(m/s)     Re   Tgas  Tsurf   q''(MW/m2)
 0.92  0.06     21      5290   293   293       0.0    (cold, quiescent)
 0.97  0.66    226    136610   296   294       0.0    (fast gas, still COLD)
 1.02  1.55   1703    132780  3119   820      39.0    IGNITES, one step
```

The cell ignites at **Mach 1.55, Re ≈ 133,000**. Over the 0.6–1.5 ms ignition
window the bore is transonic-to-supersonic (**max Mach 1.3–3.3, Re up to
330,000**), and the velocities are continuity-consistent (`ρuA` = 0.5–0.8× the
cumulative source — honest gas dynamics, not a PISO blow-up).

**Root cause:** the ignition gate sets surface heat flux by feeding these
supersonic fill velocities into the **Gnielinski/Haaland steady, low-Mach,
developed-pipe-flow correlation** (`_bare_heat_transfer_coeff` → `gnielinski_nusselt`).
At Re ≈ 10⁵ the `Nu ∝ Re^0.8` scaling returns q″ ≈ 40–78 MW/m², overwhelming the
Goodman thermal-inertia lag and lighting the whole grain in ~1 ms. The
correlation is being used 1–3 orders outside its validity. This is structural,
not burn-rate physics — and it is the dominant cause of the *sharp* spike.

### Gate-Mach-clamp experiment (negative for the principled bound)
Bounding the ignition-gate Re to the Gnielinski validity Mach:

| clamp | P_peak | ratio | full ignition |
|---|---|---|---|
| off (baseline) | 12.65 | 1.54× | 1.5 ms |
| Mach 0.30 (validity bound) | 12.56 | 1.53× | 1.6 ms |
| Mach 0.05 | 11.66 | 1.42× | 4.4 ms |
| Mach 0.02 | 10.61 | 1.29× | 11.2 ms |

At the *defensible* Mach 0.3 the clamp is a no-op. Lesson: **hot-gas presence,
not flux magnitude, is the bottleneck** — once 3100 K gas occupies a cell (which
real gas dynamics deliver everywhere in ~1.5 ms), even a validity-clamped flux
ignites the cold surface within ~1 ms. The clamp was reverted (not kept).

---

## 3. The Phase F flame-front bug (found and FIXED)

The principled way to spread ignition is the v0.7.4 Phase F flame-front (withhold
surface heating from cells the front hasn't reached). But `flame_front_velocity`
from 15→100 m/s gave **identical** results (full ignition fixed at ~4.6 ms).

**Bug:** for `forward_plenum`, the DeMar surface-flux target is `head_grain_cell`
= the *first unignited grain cell*, and the gate exempted any cell carrying DeMar
flux. As cells ignite, that target **marches down the grain** — a self-propagating
DeMar cascade that bypassed the front entirely. (This is why the v0.7.4 closeout
found flame-front "ineffective on Chunc" — it was testing a bypassed front.)

**Fix (landed in the working tree, `simulation.py`,
`_goodman_ignition_sources_and_mass`; not yet committed):** removed the DeMar-flux
exemption from the `heat_cell` gate. The seed grain cell is
always in `ignitable`, so it still lights from the head-end plume; DeMar flux is
applied only at ignitable cells, so it cannot ignite ahead of the front.

Verification (default-off, no regression):

| case | P_peak | ratio | full ignition |
|---|---|---|---|
| default (ff off) | 12.65 | 1.54× | 1.5 ms (unchanged) |
| ff v=100 | 11.18 | 1.36× | 8.9 ms (≈ 0.86 m / 100) |
| ff v=50 | 10.47 | 1.27× | 17.0 ms |
| ff v=15 | 10.32 | 1.25× | 54.5 ms |

The front is now authoritative and `flame_front_velocity` controls ignition
spread as documented — **without** the `diagnostic_disable_pyrogen_surface_heating`
crutch. 406/406 non-working-tree tests pass.

---

## 4. The residual ~1.25× erosive hump is mostly physical

Flame-spread (any speed) plus a per-cell burn-establishment ramp both bottom out
at a hard **~1.25× floor** (10.3 MPa); the establishment ramp adds essentially
nothing on top of flame spread. `diagnostic_disable_erosive` removes the hump →
**it is erosive.** Two papers (now in `docs/references/`) settle its nature:

**Mukunda & Paul 1997 (`mukunda1997.pdf`) — the hump is largely real.** Their
universal, **Hasegawa-validated** law `η = 1 + 0.023(g^0.8 − 35^0.8)·H(g−35)`,
`g = (G/ρ_p r₀)(Re₀/1000)^−0.125`, computed per-cell on Chunc:

| | η (sim, Ma) | η (Mukunda-Paul) | Ma/MP | cells g>35 |
|---|---|---|---|---|
| spike (12 ms) | 2.76 | **2.08** | 1.33 | 89/100 (g≈181) |
| plateau (600 ms) | 1.27 | 1.19 | 1.06 | 53/100 |

The threshold (g_th=35) is **not tripped-suppressed** at the transient (g≈181).
The validated universal law itself predicts η≈2.1 at the spike condition — i.e.,
substantial erosive enhancement is genuine QS physics (their Fig. 10 shows the law
reproducing real-motor erosive humps). **Ma over-predicts only ~33% vs the
universal law at high g** (agreeing within 6% at the plateau).

**Liñán & Williams 1971 (`linan1971.pdf`) — the Goodman gate is validated.**
Ignition is a thermal runaway, but *to leading order it occurs when the inert
surface temperature reaches a chemistry-set threshold* — exactly the Goodman
`T_surf > T_ignition` structure. The reactive correction is one order-unity factor
(e^b = 0.65, ~35%). So the Goodman kernel is leading-order correct; an
integrated-energy/AP-melt upgrade is a *fidelity* refinement, not a spike fix.

---

## 5. The last floor-breaking lever (relaminarization) — TESTED AND DEAD

The only remaining non-tuning candidate to suppress the erosive term transiently
was favorable-pressure-gradient relaminarization (turbulence cannot be sustained
under strong flow acceleration; erosive burning is a turbulent effect). Criterion:
`K = (ν/u²)(du/dx) > ~3×10⁻⁶` (Kays-Crawford / Narasimha-Sreenivasan).

Instrumented on Chunc at the erosive-active cells:

| | peak-erosive cell K | median K | % above 3×10⁻⁶ |
|---|---|---|---|
| spike (12 ms) | **−4.6×10⁻⁷** | 7.6×10⁻⁸ | **0%** |
| plateau (600 ms) | −1.6×10⁻⁷ | 1.4×10⁻⁸ | 0% |

K ≈ 10⁻⁷ — two orders below threshold — at both spike and plateau, and **negative**
(adverse/decelerating, du/dx ≈ −6800 s⁻¹) at the peak-erosive aft cell (x = 825 mm).
The spike flow is fully turbulent. **Relaminarization is not active; the lever is
dead.** (The sibling Beddini head-end-laminar spatial gate is dead too: the spike
is at the aft, far past any head-end laminar zone, and it is a steady effect that
would alter the plateau.)

---

## 6. The supersonic fill is a grid-divergent NUMERICAL ARTIFACT (2026-06-17)

Re-examining §2: is the supersonic ignition fill (which over-drives the convective
ignition gate) even physical? **No — it is a numerical artifact of the pressure-based
PISO solver at the cold/hot contact, and it diverges under grid refinement.**

| grid (cfl 0.3) | max Mach (fill, t<3 ms) | P_peak |
|---|---|---|
| 50 cells | 3.36 | 11.77 |
| 100 cells | 4.81 | 12.65 |
| 200 cells | **12.04** | 13.22 |

| CFL (100 cells) | max Mach (fill) | P_peak |
|---|---|---|
| 0.30 | 4.81 | 12.65 |
| 0.15 | 4.31 | 12.66 |
| 0.05 | 4.27 | 12.67 |

The fill Mach **grows without bound as dx→0** (3.4→4.8→12) but is **converged in
CFL** → a spatial artifact, not a time-stepping one. The spike tracks it (P_peak
11.8→12.7→13.2). The v0.7.4 closeout's "numerical resolution ruled out" tested CFL
and burn-cadence only — never grid refinement — and missed this.

**Mechanism (mass-balance probe).** At the Mach peak (~1.0 ms): the bore is nearly
empty (0.37 g of gas), the nozzle is NOT yet flowing (the aft is still cold, so the
nozzle BC sees `P[N-1]`≈ambient → unchoked, outflow = 0), and the igniter is forcing
~130 g/s into the head. Hot igniter gas (~10:1 density ratio vs the cold bore) fills
a near-vacuum duct with no outlet. The PISO momentum predictor has **no aerodynamic-
choking limit** on the interior port (only the nozzle THROAT chokes), and a
pressure-based scheme is **not shock-capturing**, so the velocity that conserves mass
flux across the contact discontinuity blows up — and sharpens as the grid refines.
Once the bore fills and the nozzle passes flow (~1.5 ms), Mach collapses to 0.33–0.37
and the steady state is physical and converged: **the artifact is confined to the
violent fill.** Corroborating: pressurization is ~2.4× too fast (8 MPa at 4.2 ms vs
experimental ~10 ms), and the head port (127 mm²) is smaller than the nozzle throat
(228 mm²) — the worst case for the contact blow-up.

**But the artifact is DECOUPLED from the spike (cap implemented + tested).** A
physically-grounded aerodynamic-choking limiter — `port_mach_cap`, which clamps
interior face |u| to that Mach × local sound speed in
`_piso_step_with_energy_diagnostics` — was added (default 0.0 = off) and tested at
`port_mach_cap=1.0`. The fill Mach is then bounded and grid-converged
(3.4/4.8/12 → 1.03/1.07/1.16 for 50/100/200 cells): **the velocity artifact is
fixed.** But **the spike is UNCHANGED** — P_peak 12.64 vs 12.65 at 100 cells, and
still grid-climbing 11.76/12.64/13.20. So the supersonic velocity, though a genuine
grid-divergent artifact, does NOT drive the spike. The spike is set by hot-gas
PRESENCE igniting the grain (velocity-magnitude-independent — the §2 / Mach-0.3
ignition-gate clamp was also a no-op), and the hot gas fills the short 0.8 m bore in
<1 ms even at Mach 1. The spike P_peak **grid-CONVERGES** (increments 0.88, 0.56 →
~14 MPa), so it is a **real erosive feature, not a divergent artifact.** Flame spread
(§3) reduces the spike by slowing ignition TIMING, not by touching the velocity.

**Literature.** Peretz–Kuo–Caveny–Summerfield 1973 ("Starting Transient … with High
Internal Gas Velocities") confirms high-L/D startups run genuinely transonic — but
**bounded by choking (~Mach 1)**, not Mach 12. Salita lists "compression of chamber
gases during pressurization" as a distinct, real spike mechanism — choke-limited.
The sim exceeds that physical bound, confirming the artifact.

**Disposition:** `port_mach_cap` is kept as a default-off (0.0) solver-hygiene fix —
the Mach-12 velocity field is genuinely unphysical, the cap makes the velocity / Mach
diagnostics physical, and it is a no-op at the subsonic plateau — but it is **NOT a
spike fix.** The §4 conclusion (a real erosive hump) is unchanged. Lesson: the
supersonic fill and the ignition spike are two independent phenomena; the velocity
artifact was a red herring for the spike, though a genuine artifact worth bounding.

---

## 7. Conclusion & disposition

- **No parameter-free closure removes the residual erosive hump — among levers tested
  IN ISOLATION.** Each physically-grounded lever individually — ignition timing (flame
  spread), burn buildup (establishment ramp), ignition criterion (Liñán-Williams
  validates Goodman), erosive threshold (Mukunda-Paul: not tripped), relaminarization
  (K dead), velocity (cap decoupled) — failed to close it, and the hump reads as ~⅔
  genuine QS erosive physics + ~⅓ Ma high-g over-prediction. The closeout's erosive
  verdict holds **for what was tested.** It is NOT the last word: see §8 — lever
  coupling is largely untested, the ms snap-on is still suspect, and the spike
  grid-converges to ~14 MPa (worse than the 100-cell 1.54×).
- **A separate grid-divergent VELOCITY artifact exists but is DECOUPLED from the
  spike (§6).** The pressure-based PISO blows the cold/hot contact velocity up to a
  grid-divergent Mach 12 during the fill (no interior choking limit, not
  shock-capturing). A `port_mach_cap` limiter fixes it (Mach 12 → ~1, grid-converged)
  but leaves P_peak unchanged (12.64 vs 12.65) — the spike is the erosive response to
  hot-gas-presence ignition, grid-CONVERGING to ~14 MPa (real, not an artifact). Kept
  default-off as solver hygiene; not a spike fix.
- **The v0.7.4 closeout's spike conclusion stands.** "Faithful Ma erosive / numerical
  ruled out" holds — the spike grid-converges (real). The genuine corrections from the
  re-opening are the now-fixed DeMar-cascade flame-front bug and the (decoupled,
  separately-fixed) PISO velocity artifact — neither changes the erosive verdict.
- **Landed (working tree, uncommitted):** the DeMar-cascade flame-front fix
  (correctness; opt-in Phase F now functional; default-off, no regression).
- **Recorded for future (fidelity, NOT spike fixes):** (1) de-tune Phase Z to the
  parameter-free `τ_cond = ln(10²)·α_s/r² ≈ 4.6 α_s/r²` (Lengellé thermal-wave
  time); (2) integrated-energy / AP-melt ignition criterion with the endothermic
  AP-decomposition debit (Liñán-Williams §5 plateau mechanism). See the
  `project-ignition-fidelity-candidates` memory.

---

## 8. OPEN THREADS — investigation deliberately left OPEN (next session)

This investigation is **not closed.** The residual hump survived every lever tested
*in isolation*, and Mukunda-Paul makes a chunk of it look physical — but the
ms-scale ignition snap-on remains physically suspicious, lever coupling is largely
untested, and there is an unresolved tension with experiment. Leads for a fresh
session, roughly in priority order:

**(A) The ~ms surface-heating snap-on is still implausibly fast (the prime suspect).**
Each cell heats 293→756 K in <1 ms once hot gas arrives, at q″ ~ 40 MW/m² from the
Gnielinski h_c. We proved this is hot-gas-PRESENCE-driven, not velocity-driven (Mach
cap was a no-op). But the *heating closure itself* is untested:
- Is the steady, fully-developed-pipe-flow Gnielinski Nu valid for the transient,
  developing, contact-front fill even at subsonic Mach? It likely over-predicts q″.
- The Goodman gate has NO endothermic AP-decomposition sink and NO finite ignition
  energy (recorded fidelity candidate 2). Liñán-Williams validates the surface-T
  threshold to LEADING order, but its §5 endothermic *plateau* is a real,
  flux-insensitive delay we never implemented. If the real per-cell snap-on is ~10×
  slower, the grain never reaches whole-grain-min-port simultaneously → lower transient
  G → smaller hump. **TEST:** implement the integrated-energy / AP-melt criterion +
  endothermic debit; does the snap-on slow AND the spike drop (alone, and with flame
  spread)?

**(B) Lever COUPLING is largely untested — confounding is plausible.** We eliminated
levers mostly in isolation (or in pairs). A lever that is a no-op alone may matter in
combination, or one lever's effect may be masked by another. Untested:
- velocity cap + flame spread + establishment ramp together;
- endothermic / integrated-energy criterion + flame spread;
- whether `port_mach_cap` shifts the flame-spread result (it shouldn't — verify).
Build a small Cartesian sweep over {port_mach_cap, flame_front_velocity,
tau_establishment, ignition-criterion} and look for non-additive interactions.

**(C) The spike GRID-CONVERGES to ~14 MPa — the canonical 100-cell run UNDER-resolves
it.** P_peak 11.76 / 12.64 / 13.20 for 50 / 100 / 200 cells (increments 0.88, 0.56 →
~14). So the "true" model spike is ~1.7×, WORSE vs experiment than the 100-cell 1.54×.
Leading suspect for the grid-climb: the Gnielinski ENTRANCE correction `1+(D/L)^(2/3)`
in the *erosive* burn rate uses L = x_from_head, which is grid-sensitive near the head
(first cell L ~ dx/2 → enhancement grows as dx→0). **TEST:** is the grid-climb the
entrance term? Is L=x_from_head physical or a discretization artifact for erosive
enhancement (a possible *third* artifact, distinct from the bulk erosive physics)?

**(D) The Mukunda-Paul tension is unresolved.** The validated universal law predicts
η≈2.1 (a real hump) at the sim's spike *condition*, yet experimental Chunc shows ~no
spike. Two possibilities, not yet discriminated: (i) the sim's transient G is itself
too high (the ignition/fill produces a G the real motor never reaches — ties to A), or
(ii) the real erosive response is transiently suppressed by physics no validated model
captures. **Discriminator:** reconstruct the EXPERIMENTAL transient G from the measured
P(t) + geometry and compare to the sim's; if the real G is much lower, the problem is
upstream (ignition/fill), not the erosive closure.

**(E) Shock-aware / better compressible transient treatment.** `port_mach_cap` is a
crude clamp. The pressure-based PISO is not shock-capturing and the ignition fill is a
genuine shock-tube-like blowdown. Consider a proper aerodynamic-choking source term,
artificial viscosity / flux limiting at contacts, or a density-based (Riemann)
sub-step for the transient — making the fill velocity physical *without* the clamp,
which might also surface coupling the clamp currently hides. Lower priority than A–D
for the spike, but the right long-term fix for the flow field.

**Meta-caution (user, 2026-06-17):** treat every "eliminated" lever as eliminated
IN ISOLATION only. Something is still missing; the per-cell ms snap-on is the most
likely seat of it. Next session should **re-derive from the snap-on**, not from the
erosive term — and watch for coupling between levers and unknown confounders.

---

## 9. Session 2026-06-18 — G reconstruction (thread D), Keller validation + h_c probe (thread A), igniter fix

Threads A and D were worked and one structural fix (igniter gas generation) landed.
**The investigation stays OPEN, but the snap-on / h_c seat suspected in §8 is now
largely eliminated, and the steady erosive-model fidelity (Mukunda-Paul, §8(D)
option ii) is the one surviving lever for tamping the spike.**

### 9.1 Thread D RESOLVED — the spike is the erosive feedback, NOT a G-inflation artifact
A full-resolution Chunc static-fire trace (`ThomasMach5_edited.xlsx`, 346 pts @ 7 ms —
now the canonical Chunc experimental, superseding the 59-pt digitization) made the
discriminator executable. Reconstructing the aft-port mass flux G at the spike
(t ≈ 11 ms, P_sim = 12.64 MPa):

| quantity | value | meaning |
|---|---|---|
| sim actual ρ\|u\| (erosion-driving) | 3763 kg/m²s | what Ma's erosive rate sees |
| experimental aft-G ceiling (from measured P) | 2922 kg/m²s | real motor never exceeds |
| **sim G / exp G** | **1.29×** | modest excess, not gross |
| **sim ρ\|u\| / sim quasi-steady ṁ/A** | **0.87×** | **no transient gas-dynamic inflation** |
| aft Mach at spike | **0.23** | subsonic — the supersonic fill is over |

So §8(D) option (i) "sim's transient G is itself too high" is **false in the
gas-dynamic sense**: the spike-time G is quasi-steady and only 1.29× experiment, and
the supersonic fill is decoupled from the spike — now confirmed by *direct mass-flux
measurement*, independent of the §6 `port_mach_cap` test. The 1.54× pressure spike is
a genuine **erosive feedback loop** (G→r_erosive→P→ρ→G) amplifying a modest 1.29× G
excess through Ma's super-linear high-g response. (Scripts: `c:/tmp/
chunc_G_reconstruction.py`, `chunc_G_analyze.py`; plot `artifacts/chunc_G_reconstruction.png`.)

### 9.2 Thread A — ignition criterion VALIDATED, h_c over-predicted, but the h_c lever FALSIFIED
Literature pulled to `docs/references/` (keller1966, cain2006, kulkarni1982,
jacobs1969, bircumshaw1954, + Kuo/Summerfield `Fundamentals_of_Solid-Propellant_Combustion.pdf`;
Liñán-Williams re-read firsthand via pypdf):
- **Keller-Baer-Ryan 1966** — convective AP/HTPB shock-tube ignition, 20–160 cal/cm²·s,
  to Mach 1.0 (the motor-startup analog) — **validates the surface-T criterion AND the
  756 K value**: Eq. 17 `T_ign = 300 + 286.1·F_s^0.08 ≈ 664–730 K`, ~flux-independent.
  **Hermance 1984** (Kuo/Summerfield ch. 5) + **Kulkarni-Kumar-Kuo 1982** confirm
  surface-T is the standard convective-ignition criterion and that there is *no
  universal* criterion. **The Goodman gate is not the bug** (reinforces the closeout
  from a new, on-point paper). The earlier "energy-gap" idea (vs Cain's radiant
  E_ign ~30 J/cm²) is **RETRACTED** — radiant ignition is in-depth (volumetric)
  absorbed (Cain: "surface absorption not appropriate"), not comparable to convective
  surface flux.
- **The Gnielinski h_c IS ~7× Keller's measured convective h** at matched mass flux
  (14,254 vs ~2,000 W/m²K at G = 72 g/cm²s); ignition q″ median 31 MW/m² (peak 64),
  ~10× Keller's 6.7 MW/m² ceiling → the <1 ms snap-on. So §8(A)'s "Gnielinski
  over-predicts q″" is *confirmed*.
- **BUT the h_c lever is DEAD.** A transient h_c knock-down probe — scaling the shared
  `Nu·k/D` that feeds BOTH the ignition gate AND the Ma erosive rate, gated to a time
  window (default-off scaffold; since reverted) — only **POSTPONES** the spike:
  `t_peak` tracks the window (suppress 80 ms → spike at 99 ms, 1.14×). Reaching ratio
  ~1.0 needs ~120 ms suppression (unphysical; physical flow-establishment is ~1–2 ms
  via `L_e/u ≈ 10D/u`). A parameter-free h_c establishment law restores full h_c ~8 ms
  *before* the 11 ms spike forms → no effect. **§8(A)'s integrated-energy/endothermic
  TEST is therefore de-prioritized as a spike fix**: it would slow the snap-on, but
  flame-spread + establishment-ramp already bottom at the same ~1.25× erosive floor
  (§4), Keller validates the existing threshold, and the L-W endothermic correction is
  only ~35% (e^b = 0.65). (Scripts: `chunc_hc_probe_sweep.py`, `chunc_flux_vs_keller.py`,
  `chunc_hc_decomp.py`.)

### 9.3 Igniter gas-generation over-production — FOUND + FIXED (commit bfc2f3f)
Prompted by low-mass spikes in simulation and the historical comparison with
higher real-firing charges; that empirical spike classification is not
currently qualified under protocol v0.1.
The pyrogen gas physics ignored both the cited impetus and the condensed phase:
BPNV `M = 0.030` / 100%-gas gave `R·T/M = 6869 psi·in³/g` — **2.1× the ProPep
per-charge value (3245)** and 1.37× BPNV's own cited DeMar 5000. ProPep (BPNV 25:60:15
@1000 psi): gas MW 48.2, T 2719 K, γ 1.17, **gas_mass_fraction 0.781** (condensed =
B(liq) 0.746 + BN(s) 0.557 = 21.9 g/100 g). Fix (committed): bpnv gas properties + new
`Pyrogen.gas_mass_fraction` wired into both igniter paths (solid depletes fully, only
the gas fraction pressurizes). Effect: Chunc Sutton 0.9 g spike **1.50→1.36**, as-fired
6 g **5.08→3.22**, plenum P_ig **17.6→11.5 MPa**; **Hasegawa A unchanged (6.14 MPa** —
its peak is the late progressive peak, not the ignition transient). Resolves the
low-mass over-drive but **leaves the simulated transient** (1.36×). **BALLSstick**
(bpnv, high-L/D) still produces a large simulated early peak after the fix.
Because its raw experimental pressure channel is uncalibrated, this is not an
empirically measured residual or a validation result.

> **MARKED FOR LATER IMPLEMENTATION (igniter CEA extension).** MTV (impetus 1.47×) and
> Cu/Al thermite (6.54×) carry the *same* impetus + all-gas error that BPNV did, but were
> NOT fixed (the user supplied ProPep only for BPNV). To extend commit bfc2f3f to them:
> run ProPep/CEA on each at the 1000 psi chamber standard, read off the **gas-phase mean
> MW**, **T_flame**, **γ**, and **gas mass fraction** (gas_mass / charge_mass = gas mols ×
> gas MW / 100 g), then set those four fields in `mtv.yaml` / `mtv_fast.yaml` /
> `thermite.yaml`. This is a pure correctness fix (no dogma issue — it just honours
> measured thermochemistry, like BPNV). Until then those pyrogens over-produce igniter gas
> by their listed factors. Thermite is also architecturally mismatched with the 0D-plenum
> mass-injection model (low-gas/high-flux) — see [[reference_copper_thermite_igniter]].

### 9.4 Net state — every transient / ignition lever is now exhausted
Eliminated: igniter mass/topology/IC/**gas-over-production**, ignition kernel
(**Keller-validated**), sequencing (flame-front), Beddini, numerical resolution,
velocity/gas-dynamic G inflation (**thread D, direct-measured**), transient h_c
relaxation (**probe-falsified**). The spike is the **genuine Ma quasi-steady erosive
feedback to a modestly-elevated (1.29×) transient G**; ~2/3 of the hump is real physics
(Mukunda-Paul η ≈ 2.08 at the spike g), ~1/3 is Ma's high-g over-prediction.

### 9.5 DIRECTIONS FOR TAMPING THE EROSIVE SPIKE
The remaining lever is **not** transient/ignition — it is the **steady erosive model's
high-g fidelity**. Priority order:

1. **Ma-vs-Mukunda-Paul high-g recalibration — REJECTED (user decision 2026-06-18).**
   Ma over-predicts ~1.33× vs the Hasegawa-validated Mukunda-Paul universal law
   `η = 1 + 0.023(g^0.8 − 35^0.8)·H(g−35)` at the spike (η 2.76 vs 2.08), while
   agreeing within 6% at the plateau — so a high-g correction toward MP would cut the
   spike ~1.33× without altering the plateau. **It is the only lever shown able to break
   the erosive floor.** But the user rejects it: **Ma's lack of empirical tuning figures
   is a SHIPPING requirement** — the no-tuning guarantee is what lets srm_1d ship as an
   end product that does NOT require users to quantify their own erosive burning. Any
   high-g blend/correction toward MP introduces a fitted figure (even if MP is
   Hasegawa-validated) and erodes that guarantee. **Therefore there is no remaining
   non-dogma lever, and the high-L/D ignition-transient erosive over-spike is an
   ACCEPTED, DOCUMENTED LIMITATION** — the cost of the no-tuning property. (Read-only
   Ma-vs-MP *diagnostic* comparison is fine for understanding; do not wire any
   correction.) See [[feedback_keep_ma_erosive_model]].

2. **Grid-climb forensics (thread C, still open; secondary — understanding only).** The spike
   grid-CONVERGES to ~14 MPa (worse than the canonical 100-cell). Quick check first:
   the erosive entrance term `1 + (D/L)^(2/3)` uses `L = max(x_from_head, D_hyd)` — the
   `max(·, D)` clamp already bounds it ≤ 2, so it is probably **NOT** the grid driver
   (contra §8(C)). More likely the grid-climb is **simultaneity sharpening** (finer grid
   → sharper contact front → tighter ignition window → bigger synchronized aft-G surge),
   which ties back to thread D — in which case the MP high-g correction (lever 1)
   addresses it at the source.

3. **Post-igniter-fix coupling sweep (thread B; very low — likely moot).** With the MP
   high-g lever rejected (#1), the only remaining axes are the igniter
   (gas_mass_fraction/mass) and flame_front_velocity, both of which bottom at the
   erosive floor alone. A coupling sweep could check for non-additive interactions, but
   with no floor-breaking lever left, it can at best re-confirm the documented limitation.

**Explicitly NOT to retry** (recorded so they are not re-derived): transient h_c
establishment relaxation (probe-falsified — only postpones); velocity / `port_mach_cap`
as a spike fix (decoupled, thread D direct-measured); integrated-energy / endothermic
ignition criterion *as a spike fix* (Keller validates the existing surface-T gate +
756 K; the snap-on is not the controlling lever; the L-W endothermic refinement is a
~35% fidelity item, not a spike fix).

---

Diagnostic scripts used live under `c:/tmp/` (chunc_spike_probe, chunc_fill_probe,
chunc_gate_audit, chunc_mach_audit, chunc_cell_history, run_gate_clamp_experiment,
run_flamespread_experiment, run_flamespread_combo, chunc_mukunda_compare,
chunc_K_relaminarization, chunc_fill_physics, chunc_mach_convergence); comparison
plots under `artifacts/ignition_mach_clamp/`, `artifacts/ignition_flamespread*/`.
2026-06-18 session scripts (§9): `chunc_G_reconstruction.py`, `chunc_G_analyze.py`,
`chunc_G_plot.py`, `chunc_flux_vs_keller.py`, `chunc_hc_decomp.py`,
`chunc_hc_probe_sweep.py`, `igniter_survey.py`, `igniter_impetus_test.py`; plot
`artifacts/chunc_G_reconstruction.png`. Experimental trace: `srm_1d/plotting.py`
`CHUNC_EXPERIMENTAL` is the 59-pt digitization; the full 346-pt trace is
`ThomasMach5_edited.xlsx` (user-supplied).
