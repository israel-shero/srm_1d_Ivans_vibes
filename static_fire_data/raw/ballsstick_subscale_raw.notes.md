# ballsstick_subscale_raw.csv — provenance & cleaning notes

**Status: RAW.** The logged pressure channel is not currently eligible for
absolute-amplitude or pressure-ratio scoring. Its provenance, calibration,
linearity, and ignition-time anchor remain unresolved; do **not** apply a
constant offset or treat the apparent excursion below as a validated result.

## Preliminary observation only: apparent excursion in the logged PSI column

A ratio would cancel a pure gain error and may be insensitive to a constant
offset, but those assumptions require a stable, linear channel. This file does
not establish that: PSI rises before ignition while the logged voltage is
nearly flat, and the column's derivation is unknown. The numerical ratios below
are therefore exploratory calculations from the PSI column, not validated
measurements or a scoring target.

Using the prior exploratory time origin (not independently anchored):

| feature | value |
|---|---|
| logged PSI peak | **1330 psi @ +0.070 s** |
| logged PSI post-peak trough | **885 psi @ +0.300 s** |
| logged PSI mid-burn value | **985 psi @ +1.03 s** |
| **logged peak / trough** | **1.50×** |
| **logged peak / mid-burn** | **1.35×** |

These calculations do **not** establish that the real 2″ motor had a 1.35–1.50×
pressure ratio, and they cannot yet be contrasted quantitatively with Chunc or
used to select a model success criterion. They are a reason to resolve the
channel and obtain an independent event-time anchor.

**Caveats:** it's the **2″ subscale**, not the 3″ motor; the absolute pressure is
unusable; and whether the apparent startup feature is **igniter-driven or
erosive cannot be determined from the pressure trace alone**. Until channel
linearity and timing are established, retain it as a qualitative observation,
not a validated pressure-ratio constraint.

## Provenance
Raw DAQ export from a **2″ subscale** firing of the BALLSStick design.
**No 3″-scale static fire exists** — so validating `motors/BALLSstick.ric` (3″)
against this requires both **cleaning** and a **justified subscale→3″ scaling**.

Columns: `Relative Time, Date, Time Stamp UTC, Volt, Volt, PSI, Pounds`
(2 voltage channels; PSI = pressure transducer; Pounds = load cell).

## Measured state (verified 2026-06-19)
- **417 rows**; raw time **[1388.16, 1392.32] s** — the **DAQ power-on clock**,
  not zeroed to ignition.
- **Ignition at ≈ 1388.37 s**; firing window ≈ **[1388.38, 1390.69] s** →
  **duration ≈ 2.31 s**. The record therefore contains only **~0.21 s of
  pre-ignition data.**
- **Peak = 1329.7 psi = 9.168 MPa** *as logged* (uncorrected — see below).
- **Load cell:** peak ≈ **2800 N**; integrated impulse over the firing window
  (naive baseline subtraction) ≈ **3971 N·s** vs the expected **≈ 4670 N·s**
  → **0.85×** (≈15 % low).

## Pressure channel — the real problem (read before correcting anything)
The user flagged that ambient reads far above the expected ~14.7 psi. On
inspection it is **worse than a simple offset**:

- **There is no quiescent baseline in the record.** PSI rises **monotonically
  from the very first sample**: `34.9 → 36.5 → 39.7 → 44.3 → 49.1 → … → 75.7`
  over the ~0.19 s before ignition. So there is nothing flat to tare against,
  and any "baseline" you compute is *already mid-ramp* (e.g. a median of the
  first 20 samples gives ~60.7 psi, which is meaningless).
- **PSI does not track the logged voltage.** Over that same pre-ignition span the
  Volt channels stay ~flat and merely oscillate (**0.092–0.105 V**, ~0.013 V
  span) while **PSI climbs ~41 psi monotonically.** A linear transducer map
  cannot produce that. So the PSI column is **not** a straightforward function
  of the logged Volt columns — its derivation (different channel? filter/lag?
  mis-scaled range? drift?) is **unknown**.
- Even the earliest sample (**34.9 psi**) is well above ambient.

**⇒ Do not subtract a constant offset and do not trust the logged PSI scale**
until the column's provenance and the transducer calibration are established.
The most useful next input is **more pre-trigger DAQ data** (enough quiescent
record to see a real zero), plus the transducer's range/calibration and how the
PSI column was derived from the raw volts.

## Known issues (user-reported + confirmed above)
1. **Time not zeroed** — starts at the DAQ power-on timestamp (~1388 s).
2. **Pressure channel provenance/calibration unresolved** — see above. This is
   the blocking issue.
3. **Load cell suspect.** Motor mass is **not** subtracted (and the static weight
   *drifts* as propellant burns, so a single tare is not exact); integrated
   impulse comes out 0.85× the expected ~4670 N·s.
4. **Subscale.** 2″ → 3″ scaling is required and must be **physically justified**
   (Kn / c\* / L\* reasoning), not eyeballed — see `feedback_no_unfounded_smoothing`.

## Cleaning TODO (in order)
1. **Establish the pressure channel's provenance + calibration** (blocking).
   Ideally obtain a longer pre-trigger record.
2. Zero the time axis to ignition; trim to the firing window.
3. Resolve the load cell (tare / motor-mass handling; reconcile the impulse gap).
4. Decide + justify the subscale→3″ scaling.
5. Emit a cleaned `../ballsstick_subscale.csv`, add its load recipe to
   `../README.md`, and record every step + justification here.

**Until 1–4 are resolved this trace cannot support any claim about the 3″
BALLSStick.** It is committed for provenance and so the cleaning is reviewable.
