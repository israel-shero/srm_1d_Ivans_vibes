# Project Checkpoint — 20 November 2026

## Repository States (both clean, pushed to forks)

### srm_1d
- **Branch**: `integration/openmotor-packaging` (7a8908f)
- **Remote**: `ivans-vibes` → https://github.com/israel-shero/srm_1d_Ivans_vibes.git
- **Upstream**: `origin` → https://github.com/eJar0k/srm_1d.git
- **HEAD**: `7a8908f docs: freeze scoring protocol v0.1 and register Zerox case windows`

### openMotor-compatible
- **Branch**: `integration/srm1d-packaging` (efb0c81)
- **Remote**: `fork` → https://github.com/israel-shero/openMotor-srm1d.git
- **Upstream**: `origin` → https://github.com/eJar0k/openMotor.git
- **HEAD**: `efb0c81 Make solver config editors responsive`

---

## Prediction Plan — Current Gate

**Recommended next gate**: Complete source-only case reviews under protocol v0.1

### Task Status Summary

| Task | Status | Blocking Evidence |
|------|--------|-------------------|
| T01 Workspace audit | ✅ Mostly complete | `workspace_verification_2026_11_20.json` |
| T02 Source reconciliation | 🟡 Partial | `closure_source_evidence.json` |
| T03 Data audit | 🟡 Partial | `validation_case_registry.json` (protocol v0.1 frozen) |
| T04 Baseline reproduction | 🟡 Partial | Matched 106/206-cell MP runs at CFL 0.075 |
| T05 Numerical verification | 🟡 Partial | Peak P/t/I <0.061%; min P, Mach, clipping >2% |
| T06 Driver/closure isolation | 🟡 Partial | Common-state limiting checks at 5 MPa |
| T07 Candidate decision | ⏳ Pending | Requires T06 completion |
| T08 Prototype | ⏳ Pending | Requires T07 decision |
| T09 Validation | ⏳ Not ready | Requires held-out cases |
| T10 GUI integration | 🟡 Partial | Headless/GUI parity unverified |
| T11 Release review | ⏳ Pending | Requires all gates |

---

## Case Registry — Scoring Eligibility (protocol v0.1 frozen)

| Case | Absolute Timing | Absolute Pressure | Shape | Fit/Holdout |
|------|-----------------|-------------------|-------|-------------|
| **Zerox** | ✅ Candidate (Py1 Fire, 70ms bracket) | ❌ Blocked (calibration) | ✅ Startup window registered (0–0.36s) | Unassigned |
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
| **Total** | **128** | Matches plan claim of ~118 solver/adapter/plugin tests |

## Resume Instructions for Next Agent

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