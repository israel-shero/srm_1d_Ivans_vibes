# static_fire_data/ — experimental validation traces

Measured static-fire traces used to validate the simulator's pressure/thrust
history. **Primary consumer: the ignition-spike investigation**
(`docs/v0_7_4/SPIKE_REOPEN_BRIEF.md`), whose goal is to match empirical data
**across the full short→long L/D range** — so the value here is *coverage across
motors*, not just one trace. Also feeds the example plot overlays.

---

## How to load these (there is no single schema — read this)

`srm_1d.plotting.load_experimental_csv` is **generic**: you tell it which
columns hold time and pressure and what unit the pressure is in, and it returns
`{'time' (s), 'pressure' (MPa), 'label'}`.

```python
load_experimental_csv(path, time_col=0, pressure_col=1,
                      pressure_unit='MPa',   # 'Pa' | 'kPa' | 'MPa' | 'psi'
                      skip_header=1, label=None)
```

**The files here do NOT share a column layout.** Use these exact recipes:

| File | Recipe |
|---|---|
| `thomas_chunc_firing.csv` | `time_col=0, pressure_col=1, pressure_unit='MPa'` |
| `zerox_data.csv` | **`time_col=1, pressure_col=0`**, `pressure_unit='MPa'` — ⚠️ **columns are REVERSED** (pressure first). Loading with the defaults silently swaps time and pressure. |
| `pathfin54.csv` | `time_col=0, pressure_col=2, pressure_unit='Pa'` (col 1 is force in N) |
| `raw/ballsstick_subscale_raw.csv` | **do not load directly** — raw DAQ, needs cleaning (see below) |

A third form exists in code: hand-digitized dicts in `plotting.py`
(`HASEGAWA_MOTOR_A_EXPERIMENTAL`, `ZEROX_EXPERIMENTAL`, `CHUNC_EXPERIMENTAL`) —
`{'label', 'time' (s), 'pressure' (MPa), 'time_offset'?}`, where `time_offset`
shifts the experimental ignition to sim `t=0`. These are coarse (30–60 pt);
**prefer the CSVs above where one exists.**

---

## Inventory (last broad review 2026-06-19; Zerox workbook reviewed 2026-11-19)

| Motor | L/D | File / dict | Points | State | Notes |
|---|---|---|---|---|---|
| **Chunc / machbusterNew** | high | **`thomas_chunc_firing.csv`** (+ `.xlsx` source) | **346** | ✅ clean | t=[0.010, 2.353] s; plateau ≈ **8.77 MPa**; P_peak 8.894 @ 0.397 s → **spike ratio 1.015 = essentially NO ignition spike.** The high-res trace behind the G-reconstruction (REOPENED §9.1). **The key spike-diagnostic asset.** Supersedes the 59-pt `CHUNC_EXPERIMENTAL` dict. |
| **Zerox** | mid | `zerox_data.csv` + [`zerox_high_resolution.csv`](zerox_high_resolution.csv) + `Zerox Data.xlsx` source | 58 coarse / 630 native event samples | ⚠️ pressure calibration review pending | The workbook's `FL88` sheet has a `Py1 Fire` trigger and a denser pressure trace. The first observed trigger-high is at flight time 0.009 s; the preceding low is at -0.061 s. The old 58-point CSV is a workbook-selected downsample: using the workbook's own baseline formula, it reports a 3.9938 MPa peak 81 ms later and 5.7% lower than the first high-resolution peak. The extracted CSV preserves source `pTap 1(psi)` values without smoothing or resampling; the processed time-zero average is excluded. Calibration and the full tail remain under review. |
| **Hasegawa A** | ~42 (high) | `HASEGAWA_MOTOR_A_EXPERIMENTAL` (plotting.py) | 36 | ✅ clean (digitized) | Ma 2020 Fig. 10. Its peak is the *late progressive* peak, not an ignition spike. |
| **Pathfinder 54** | — | `pathfin54.csv` | 409 | ✅ clean | Has force (N) + pressure (Pa). |
| **BALLSStick (2″ subscale)** | high | **`raw/ballsstick_subscale_raw.csv`** | 417 | ⚠️ RAW: qualitative observation only; pressure ratios are not scoreable | The PSI column has an apparent startup excursion, but channel provenance/calibration, linearity, and event time remain unresolved. Do not treat the reported 1.35–1.50× ratios as validated or scoreable. Load-cell impulse is also unresolved, and the DAQ clock is not anchored. **No 3″ firing exists**; scaling to the 3″ motor is unsupported. See `raw/ballsstick_subscale_raw.notes.md`. |
| **54-2800 ST2.0** | — | *(CSV missing)* | — | ❌ absent | `examples/run_pathfin54st2_0.py` + `motors/2025.12.25 54-2800 ST2.0.ric` exist; the CSV does not. |
| **short / low-L/D** | low | *(none)* | — | ❌ absent | Held for now (no clean data). **Still the key coverage gap** — a fix must be shown not to break the near-zero-spike low-L/D case. |

**L/D coverage:** candidate high-L/D evidence includes Chunc's hi-res trace and
Hasegawa A's digitized trace; BALLSStick remains raw and qualitative-only.
**Low-L/D remains uncovered.**

### Zerox source review (2026-11-19)

The 630-row [`zerox_high_resolution.csv`](zerox_high_resolution.csv) contains
native pressure samples from the workbook's processed event interval
(flight time 0.009–8.153 s), with time shifted only to the first observed
`Py1 Fire` high. The previous low sample at -0.061 s brackets that transition
to 70 ms. This is a sampled trigger anchor, not an exact ignition instant.

The workbook's processed gauge-pressure column subtracts the mean of 1,564
pre-trigger `pTap 1(psi)` samples (`AVERAGE(K2:K1565)`, 9.8441379 psi) and
multiplies by 0.00689476 MPa/psi. The extracted CSV keeps the labeled raw psi
channel so the pressure-scale/calibration question is not hidden by that
processing. The workbook does not include independent transducer calibration
metadata. Its processed time-zero row is an averaged synthetic baseline, not a
native sample, and was not copied into the extracted trace.

The original 58-point CSV matches the workbook's selected every-tenth series
plus that synthetic baseline. Under the workbook's own baseline formula, its
peak is 3.9938 MPa at 0.421 s after the trigger, versus 4.23575 MPa with the
first maximum observed 0.340–0.360 s after the trigger in the denser trace.
That sampling difference is not a model score. The extracted series stops at
8.153 s, while the workbook retains later raw DAQ samples; do not score
tail-off until that remaining interval is reviewed.

---

## Adding a trace

1. **Cleaned →** drop `<motor>.csv` at the top level; add its load recipe to the
   table above (columns/units are per-file, so this is mandatory).
2. **Raw →** put it in [`raw/`](raw/) and write a sibling
   `<name>.notes.md` recording provenance + every processing step (calibration,
   offsets, filtering, scaling). Justify scalings physically — do not eyeball
   (`feedback_no_unfounded_smoothing`).

## Known issues (pre-existing, flagged not fixed)
- `examples/run_pathfin54.py` points at `…/pathfin54.csv/` (**trailing slash** —
  won't load).
- `examples/run_pathfin54st2_0.py` points at a nested `…/srm_1d/srm_1d/
  static_fire_data/…` path (wrong depth + missing file).

*(`.ric` motor files are openMotor's format — never hand-edit. Data files and
this README are fair game.)*
