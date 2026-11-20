# v0.7.4 — Ignition-Spike Investigation Close-Out (2026-06-01)

> **UPDATE (2026-06-16/17): re-opened — see `IGNITION_SPIKE_REOPENED.md`.**
> The core conclusion below (no parameter-free closure for the residual erosive
> hump) was re-tested from multiple new angles and **stands, reinforced**
> (Mukunda-Paul: the hump is largely real QS physics; relaminarization lever
> tested and dead). Two refinements to this document's claims: (1) the *sharp
> early* spike is structural (supersonic-fill convective ignition gate), distinct
> from the residual hump; (2) "flame-front can't suppress at physical speeds" was
> partly an artifact of a now-fixed bug — the Phase F front was silently bypassed
> for `forward_plenum` by a marching DeMar-flux target.

> **Evidence status update (2026-11-20):** The Chunc source workbook and CSV
> contain only time and pressure, start at 0.010 s / 8.466 MPa, and provide no
> independent event marker, pre-event baseline, or pressure-calibration
> metadata. Under protocol v0.1, the historical “no spike” classification and
> model-to-experiment spike difference below are not registered validation
> conclusions. See
> [`validation_case_registry.json`](../v0_8_0/data/validation_case_registry.json).

**Historical verdict (pre-protocol v0.1):** the residual Chunc ignition
over-prediction was interpreted as the genuine, faithful Ma-2020 quasi-steady
erosive response to the real transient mass-flux during fast ignition of a
high-L/D motor — a regime Ma's paper explicitly excluded from every error
figure and never benchmarked. Under the project's no-tuning constraint, the
analysis found no literature-supported mechanism (Ma / King / Beddini) that
removed it. The simulator's agreement with Ma's own validation motors remains
repository-reported analysis, not a qualification of the Chunc measurement.

This document is the durable record of the differential diagnosis that
reached that verdict, so the question is not re-opened from scratch.

## The question

The historical comparison for Chunc (`machbusterNew`, high-L/D) described the
trace as rising with no spike to a flat ~8.8 MPa plateau while the simulation
produced a transient peak. The source now reviewed has no independent event
time, pre-event baseline, or pressure calibration and begins at 8.466 MPa at
0.010 s. Whether this trace resolves the startup event—and therefore the
reported peak/plateau error—remains open.

## What was FIXED and kept (committed, 69a6c06)

- **MTV burn-rate recalibration** (`mtv.yaml` a=3.0e-5→4.4e-5, n=0.50→0.35,
  Kubota & Serizawa 1987). The old seed over-predicted MTV burn rate 6–8×.
  Cut the Chunc head_basket spike 2.02×→1.55× at the Sutton 0.9 g charge.
  Old seed preserved as `mtv_fast.yaml` for A/B. **Necessary but not
  sufficient** — at the as-fired 6 g charge the spike is 3.17×.
- **Realistic basket cartridge geometry** (`basket_fill_fraction=0.5`,
  `pellet_packing_fraction=0.60`) replacing the solid-puck `L_cart`.
  Physical-correctness fix; **proven NEUTRAL on the spike** (total pyrogen
  mdot is conserved regardless of cartridge cell-distribution).

## The elimination chain (every candidate, with the test that killed it)

| candidate | test | verdict |
|---|---|---|
| igniter mass / concentration | mass×cartridge sweep | spike scales with mass, but distribution NEUTRAL |
| igniter **topology** | forward_plenum vs head_basket on Chunc | NEUTRAL (1.50× vs 1.56× @0.9g; 2.46× vs 3.17× @6g) |
| cold-bore **B.0 IC** | hot-bore (3105 K) vs cold-bore (293 K), Chunc 6 g | IDENTICAL (20.55 vs 20.62 MPa) |
| our **ignition kernel / implementation** | Hasegawa A/B/C (Ma's exact validation motors) | **FAITHFUL: A 1.15×, B 0.68× (no spike), C no startup spike** |
| **MTV burn-rate magnitude** | Kubota recal | fixed; necessary, not sufficient |
| ignition **sequencing** (flame front) | v_front sweep 3–200 m/s | only suppresses at ≥80 ms ignition (v≤10 m/s) = the rejected over-slow/slanted regime; convective speeds (30–200 m/s, 5–28 ms snap-on) give baseline 3.0–3.2× |
| transient **erosive lag (b)**, *derived* | Beddini 1986+1978 read | **no derivable τ exists** — only a spatial x/R gate, which (a) is steady → alters plateau, (b) doesn't suppress the *aft* spike, (c) like any flow-*state* gate reads "developed" at the high-flow spike. A working temporal lag needs a TUNED τ → violates no-tuning. |
| **numerical resolution** | every-step burn+geom update + CFL 0.15 | RULED OUT (3.17× → 3.17×, 852k steps; baseline dt≈1.9 µs already far finer than the transient) |

## The two load-bearing literature findings

1. **Ma 2020 (DOI 10.1155/2020/8889333)** — the erosive model of record —
   is **quasi-steady with no transient/lag/establishment term**, validated
   against Hasegawa A/B/C *full-motor transient traces using the model
   alone*. Its transient solver DOES produce a startup peak, but **every
   quantitative accuracy figure (≤10%) explicitly EXCLUDES the ignition
   process; the startup-peak magnitude is never benchmarked.** Ma cites
   King extensively; Beddini is never cited.
2. **Beddini 1986 + 1978** — the injection-driven-BL literature — gives
   **no parameter-free establishment *time***. The only derivable result
   is a spatial `x/R ≈ 5–10` turbulence-buildup window from the closed head
   end, which is steady (would alter the plateau) and does not gate the aft
   spike. The mean-flow transition criterion `Re_c,tr ∝ Re_s` is explicitly
   non-universal (tuned pseudoturbulence σ_v = 0.035–0.078). So a transient
   erosive closure cannot be derived without a tuned constant.

## Why high-L/D is the differentiator

Same faithful implementation, same `forward_plenum`, same knobs, same IC:
Hasegawa A is 1.15×, Chunc is 1.50–2.46×. The only difference is the motor.
Chunc's **high L/D** means that when the bore lights ~simultaneously during
the fast fill (~3 ms; real motors "snap on" in tens of ms — also fast), the
long port accumulates a large aft mass-flux G, and Ma's QS erosive term
fires hard. Hasegawa A's shorter port never builds that transient G, so the
*same* erosive term stays quiet. Ignition sequencing can only desynchronize
this by slowing ignition to ≥80 ms (unphysical, contradicts the snap-on),
because the aft G is the *sum* of all upstream mass addition. Hence the
spike is the erosive *response*, not ignition timing — and the response is
Ma-faithful.

## Decision

**Documented as a known QS-erosive limitation (user decision 2026-06-01).**
No tuned closure added; the no-tuning dogma is intact. The MTV burn-rate
and basket-geometry fixes stand (committed). Re-open only if a *cross-motor,
non-tuned* transient erosive closure ever emerges in the literature.

If the dogma is ever relaxed, the single remaining lever is an opt-in
(default-OFF) erosive-establishment relaxation with a τ_BL tuned within the
Beddini `L_e/u` bound and latched so steady-state is provably untouched —
explicitly acknowledged as one tuned constant. Not implemented.

## References pulled during this investigation (in `docs/references/`)

- `kubota1987.pdf` — MTV burn rate (DOI 10.2514/3.22990)
- `king1992.pdf` — erosive threshold = turbulent transition (NTRS 19920001907)
- `beddini1986.pdf` — injection-induced porous-duct flow development (DOI 10.2514/3.9522)
- `beddini1978.pdf` — reacting-TBL erosive burning (AIAA J 16(9):898)
- Greatrix 2010 (DOI 10.3390/en3111790, open access) — not pulled (MDPI block)
