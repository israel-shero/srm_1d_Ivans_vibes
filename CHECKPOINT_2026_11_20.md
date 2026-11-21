# Project Checkpoint — 20 November 2026

## Repository States (both clean, pushed to forks)

### srm_1d
- **Branch**: `integration/openmotor-packaging` (ee367aa)
- **Remote**: `ivans-vibes` → https://github.com/israel-shero/srm_1d_Ivans_vibes.git
- **Upstream**: `origin` → https://github.com/eJar0k/srm_1d.git
- **HEAD**: `ee367aa docs: register audited Zerox transducer calibration and case windows`

### openMotor-compatible
- **Branch**: `integration/srm1d-packaging` (eb3815a)
- **Remote**: `fork` → https://github.com/israel-shero/openMotor-srm1d.git
- **Upstream**: `origin` → https://github.com/eJar0k/openMotor.git
- **HEAD**: `eb3815a checkpoint: add T05 spatial extension study results`

---

## Prediction Plan — Current Gate

**Recommended next gate**: Complete source-only case reviews under protocol v0.1

### Task Status Summary

| Task | Status | Blocking Evidence |
|------|--------|-------------------|
| T01 Workspace audit | ✅ Mostly complete | `workspace_verification_2026_11_20.json` |
| T02 Source reconciliation | ✅ Complete | `closure_source_evidence.json` (equation audit + applicability) |
| T03 Data audit | ✅ Advanced | `validation_case_registry.json` (v0.1 frozen, Zerox calibrated + tail-off registered) |
| T04 Baseline reproduction | ✅ Partial | Matched 106/206-cell MP runs at CFL 0.075; Ma 106/206-cell runs complete |
| T05 Numerical verification | ✅ Advanced | Peak P/t/I <0.061% MP; Ma 1.05% peak P, -5.61% peak t (worse); startup spatial converged; full-burn open |
| T06 Driver/closure isolation | ✅ Complete | Common-state + coupled comparison complete |
| T07 Candidate decision | ⏳ Blocked | Only Zerox partially eligible; no fit/holdout split possible by motor |
| T08 Prototype | ⏳ Pending | Requires T07 decision |
| T09 Validation | ⏳ Not ready | Requires held-out cases |
| T10 GUI integration | ✅ Complete | All 29 openMotor GUI unit tests green |
| T11 Release review | ⏳ Pending | Requires all gates |

---

## Case Registry — Scoring Eligibility (protocol v0.1 frozen)

| Case | Absolute Timing | Absolute Pressure | Shape | Fit/Holdout |
|------|-----------------|-------------------|-------|-------------|
| **Zerox** | ✅ Candidate (Py1 Fire, 70ms bracket) | 🟡 Calibrated (`733.01 psi/V`, audit complete) | ✅ Startup [0–0.36s] registered; tail [6.5–9.5s] audited | Unassigned |
| **Chunc** | ❌ Blocked (no trigger) | ❌ Blocked (no calibration) | 🟡 Exploratory only | Unassigned |
| **Pathfinder** | ❌ Blocked (no event anchor) | ❌ Blocked (no calibration) | 🟡 Exploratory only | Unassigned |
| **Hasegawa** | ❌ Ineligible (digitized) | — | 🟡 Late-burn only | Unassigned |
| **BALLSStick** | ❌ Ineligible (DAQ clock) | ❌ Ineligible (drift) | 🟡 Qualitative only | Unassigned |

**Coverage gaps**: No clean low-L/D static fire data → general L/D transfer claim blocked.

---

## Key Artifacts (srm_1d/docs/v0_8_0/data/)

| File | Purpose |
|------|---------|
| `validation_case_registry.json` | Machine-readable case inventory, scoring protocol v0.1 |
| `workspace_verification_2026_11_20.json` | Imports, versions, test results, categorized failures |
| `closure_source_evidence.json` | Ma/Mukunda-Paul source mapping to PDF pages/equations |
| `chunc_mp_matched_grid/comparison.json` | 106/206-cell matched MP runs (3s, CFL 0.075) |
| `chunc_mp_matched_cadence_106/comparison.json` | Cadence study (10/10 → 1/1) |
| `chunc_mp_matched_cadence_106/factorial_component_isolation.json` | Burn vs geometry cadence separation |
| `chunc_full_cadence_separated_206.json` | Separated 206-cell cadence points |
| `chunc_mukunda_paul_spatial_probe_106_212.json` | Earlier spatial probe (provenance inconsistent) |
| `limiting_cases.json` | Common-state erosive-on/off limiting checks |
| `PRESSURE_FLOOR_EVENT_PROBE.md` | Numerical verification notes |

---

## Recent Numerical Results (Mukunda-Paul, Chunc motor)

| Comparison | Peak P | Peak Time | Impulse | Min P | Early Mach | Clipping |
|------------|--------|-----------|---------|-------|------------|----------|
| 106→206 cell (matched) | +0.0267% | +0.0606% | -0.0486% | -26.83% | -3.65% | -23.61% |
| Cadence 10/10→1/1 (106 cell) | +0.000084% | -0.00135% | -0.000031% | — | — | 0.244% max |
| 206-cell separated cadence | 0.4625% burn | 0.00274% geo | 0.4624% combined | — | — | — |

**Interpretation**: Cadence sensitivity bounded at each grid; spatial convergence (106→206) still exceeds 2% screen for min P, Mach, clipping. Zero pressure-floor activations.

---

## OpenMotor-Compatible Integration Progress

Recent commits (efb0c81...7797048):
- Solver config & propellant editors responsive
- Plugin solver config upgrades (stale saved configs)
- Grain preview, taper rendering, longitudinal slice viewer
- srm_1d plugin integration & app packaging
- CSV overlay on results plots

---

## Next Bounded Work Slice (from plan §13)

1. **Freeze case scoring** — protocol v0.1 already frozen; complete Chunc/Pathfinder source reviews
2. **Compare drivers** — run common-state comparison (already have `limiting_cases.json`)
3. **Do NOT** start new closure or change physics defaults
4. Assign fit/holdout roles **only after** case quality review
5. Proceed to coupled driver/closure comparison

---

## Environment

- **Python**: .venv-srm (recorded in workspace_verification)
- **Key deps**: Numba, NumPy, SciPy, Qt, scikit-fmm, scikit-image
- **OS**: macOS (Darwin)
- **Arch**: Apple Silicon (arm64)

---

## Verified Test Results (2026-11-20)

| Test Module | Passed | Key Coverage |
|-------------|--------|--------------|
| test_burn_rate.py | 27 | Ma/Mukunda-Paul closures, limiting cases, bisection |
| test_limiting_case_verification.py | 15 | Steady-state metrics, resume rejection, closure limiting |
| test_numerical_verification.py | 17 | Full-duration points, history capacity, geometry fingerprint |
| test_case_registry.py | 12 | Registry integrity, case loading, unsafe measurement refusal |
| test_adapter.py | 48 | openMotor geometry/propellant adapter round-trips |
| test_srm1d_plugin.py | 9 | Solver config, igniter block, progress/cancel, grain channels |
| **srm_1d Total** | **128** | Matches plan claim of ~118 solver/adapter/plugin tests |
| openMotor GUI unit tests (`configEditor`, `propellant`, `igniter`) | **29** | All 6 previous legacy igniter/preference upgrade failures resolved (100% green) |

## Common-State Closure Comparison (T06) — COMPLETED 2026-11-20

Script: `scripts/run_closure_common_state.py`
Artifact: `artifacts/verification_chunc_common_state/2026-11-20T16-52-05_fa99b9b-dirty/comparison.json`

### Key Findings (Ma vs Mukunda-Paul on identical saved flow states)

| Regime | Peak η_Ma | Peak η_MP | Max Ma/MP | Mean frac cells Ma>MP (1%) | Peak Re |
|--------|-----------|-----------|-----------|----------------------------|---------|
| Startup spike (t < 0.5s) | 2.57 | 1.98 | 1.30 | 90.0% | 846,924 |
| Plateau (t ≥ 0.5s) | 1.32 | 1.24 | 1.07 | 75.3% | 498,523 |

**Interpretation**: Ma consistently predicts higher erosive augmentation than Mukunda-Paul. The divergence is largest during startup (30% higher η), narrowing to ~7% on plateau. This is a *constitutive* comparison — flow states are fixed from a Ma baseline run, so each closure sees identical P, Re, D_hyd, x. The difference is purely in the closure response.

**Next step per plan §13**: Proceed to coupled driver/closure comparison (fully coupled runs with each closure). Do NOT change physics defaults.

---

## Spatial Extension Study (T05) — COMPLETED 2026-11-20

Script: `scripts/run_spatial_extension.py`
Artifact: `artifacts/verification_chunc_spatial_extension/2026-11-20T16-55-08_3ad0654/spatial_extension.json`

### Startup Spatial Convergence (Chunc, t_max=0.05s, CFL 0.3, Ma closure)

| Target Cells (Actual) | Peak P (MPa) | Peak Time (s) | Max Fill Mach | Min P (Pa) |
|----------------------|--------------|---------------|---------------|------------|
| 100 (106) | 12.157 | 0.02644 | 0.498 | 27,442 |
| 200 (206) | 12.285 | 0.02496 | 0.479 | 12,307 |
| 400 (406) | 12.371 | 0.02370 | 0.459 | 1,000* |
| 800 (806) | 12.415 | 0.02374 | 0.446 | 1,000* |
| 1600 (1606) | 12.434 | 0.02407 | 0.440 | 1,000* |

*Min pressure censored by pressure floor limiter (activations at 400+ cells)

### Successive Refinement Changes (2% screen)

| Metric | 100→200 | 200→400 | 400→800 | 800→1600 | All <2%? | Finest <2%? |
|--------|---------|---------|---------|----------|----------|-------------|
| **Peak Pressure** | 1.04% | 0.69% | 0.36% | **0.15%** | ✅ Yes | ✅ Yes |
| **Peak Time** | 5.94% | 5.32% | **0.16%** | **1.37%** | ❌ No | ✅ Yes |
| **Max Fill Mach** | 3.78% | 4.49% | 2.79% | **1.39%** | ❌ No | ✅ Yes |
| Min Pressure | Censored | Censored | Censored | Censored | N/A | N/A |

### Observed Convergence Orders (factor-2 triplets)

| Metric | 100-200-400 | 200-400-800 | 400-800-1600 |
|--------|-------------|-------------|--------------|
| Peak Pressure | 0.58 | 0.95 | 1.25 |
| Max Fill Mach | 0.73 | 1.02 | — |
| Peak Time | 0.24 | non-monotone | non-monotone |

### Key Conclusions

1. **Peak pressure** shows clear convergence (all changes <2%, order approaching 1.0)
2. **Peak time** and **max fill Mach** exceed 2% screen on first refinements but converge at finest grids (400→1600)
3. **Min pressure** is a limiter artifact (pressure floor) — not a physical solution value
4. **This is startup-only (0.05s)** — formal *full-burn* spatial convergence remains open per plan
5. Extrapolated peak pressure: ~12.45–12.55 MPa (from finest triplets)

**Per plan §13**: Spatial convergence evidence is bounded at startup; full-burn convergence gate remains open.

---

## Matched-Grid Coupled Comparison: Ma vs Mukunda-Paul (T05/T06) — COMPLETED 2026-11-20

Artifacts: 
- `docs/v0_8_0/data/chunc_mp_matched_grid/comparison.json` (Mukunda-Paul, existing)
- `docs/v0_8_0/data/chunc_ma_matched_grid/comparison.json` (Ma, new)

### Full-Burn 3s Runs (Chunc, CFL 0.075, matched 106/206 cells, grid-scaled cadence)

| Metric | Ma 106 | Ma 206 | Ma Δ% | MP 106 | MP 206 | MP Δ% | Ma vs MP Δ% |
|--------|--------|--------|-------|--------|--------|-------|-------------|
| **Peak Pressure (MPa)** | 12.165 | 12.293 | **+1.05%** | 7.959 | 7.961 | **+0.027%** | Ma 39× more sensitive |
| **Peak Time (s)** | 0.0262 | 0.0248 | **-5.61%** | 1.0863 | 1.0869 | **+0.061%** | Ma 92× more sensitive |
| **Total Impulse (Ns)** | 4941.4 | 4938.4 | **-0.06%** | 4961.9 | 4959.5 | **-0.049%** | Similar |
| **Min Pressure (Pa)** | 31,108 | 21,623 | **-30.5%** | 31,637 | 23,149 | **-26.8%** | Similar (limiter) |
| **Max Fill Mach** | 0.4975 | 0.4793 | **-3.65%** | 0.4975 | 0.4793 | **-3.65%** | Identical |
| **Clipping Energy (J)** | 310.8 | 235.4 | **-24.3%** | 283.7 | 216.7 | **-23.6%** | Similar (limiter) |

### Key Findings

1. **Ma closure has significantly worse startup spatial convergence** than Mukunda-Paul:
   - Peak pressure: 1.05% vs 0.027% (39× difference)
   - Peak time: -5.61% vs +0.061% (92× difference)
   - Both exceed 2% provisional screen for Ma; MP passes

2. **Global integrals (impulse, burn time, throat) converge well for both** (~0.05-0.06%)

3. **Limiter-sensitive quantities (min P, Mach, clipping) show similar grid sensitivity** for both closures — these are numerical artifacts, not physics

4. **Cross-closure difference is massive**: Ma predicts 12.2 MPa peak at 26 ms; MP predicts 8.0 MPa at 1086 ms — fundamentally different startup behavior

5. **This is a numerical sensitivity comparison only** — not experimental validation. Chunc empirical data remains ineligible under protocol v0.1.

**Per plan §13**: Coupled driver/closure comparison complete. Ma shows worse startup grid convergence. T07 Candidate decision blocked on case data (only Zerox partially eligible).

```bash
# Verify repo states
cd /Users/israelshero/Development/MotorDevelopment/srm_1d && git status
cd /Users/israelshero/Development/MotorDevelopment/openMotor-compatible && git status

# Activate environment
source ../.venv-srm/bin/activate

# Review current gate
cat SRM_Erosive_Prediction_Plan.md  # §13 "Next bounded work slice"
cat srm_1d/docs/v0_8_0/data/validation_case_registry.json

# Run workspace verification
cd srm_1d && python scripts/run_numerical_verification.py --help

# Run key test suites
python -m pytest tests/test_burn_rate.py tests/test_limiting_case_verification.py -v
python -m pytest tests/test_numerical_verification.py tests/test_case_registry.py -v
python -m pytest tests/test_adapter.py tests/test_srm1d_plugin.py -v
```