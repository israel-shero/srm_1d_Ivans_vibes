# Zerox Pressure Transducer Calibration & Full-Tail Audit

**Date**: 2026-11-20  
**Data Source**: `srm_1d/static_fire_data/Zerox Data.xlsx` (Sheet `FL88`, SHA-256: `0c2b3810a81dc76f6bac60a8156fab7c537ed772d98cc1b629b22a1204ad9662`)  
**Scope**: Resolve the open data audit gaps recorded in `validation_case_registry.json` (protocol v0.1):
1. `zerox_pressure_calibration` (missing / unresolved channel provenance)
2. `zerox_complete_tail` (unregistered tail-off window, trace previously truncated at 8.153 s)

---

## 1. Primary Channel Inventory & Transducer Transfer Function

The workbook contains 10,790 raw logged DAQ samples spanning Flight Time `-10.000 s` to `+120.000 s`.

Key raw columns:
- Column A: `Flight Time(s)`
- Column J: `Load Cell 1(N)` (thrust)
- Column K: `pTap 1(psi)` (engineering unit pressure)
- Column L: `pTap 1(volts)` (raw voltage from transducer signal pin)
- Column Q: `Py1 Fire(0/1)` (firing command / trigger)

### Transducer Transfer Function Fit

A linear regression across all 10,790 samples confirms an exact linear relationship between raw voltage and engineering pressure:

$$\text{Pressure [psi]} = 733.011321 \times \text{Voltage [V]} - 721.016233$$

- **Maximum absolute residual**: $3.89 \times 10^{-2}\text{ psi}$ (negligible roundoff).
- **Zero-pressure intercept**: $V_0 = 0.9836\text{ V}$.
- **Sensitivity slope**: $733.01\text{ psi/V}$ ($5.0539\text{ MPa/V}$ or $50.539\text{ bar/V}$).
- **Sensor Type**: Corresponds to an industrial 1.0–5.0 V (or 0.5–4.5 V) ratiometric pressure transducer with a ~3000 psi full scale. At 1.000 V, output is 12.00 psi (~atmospheric pressure at test elevation).

---

## 2. Pre-Trigger Baseline Calibration

Pre-trigger DAQ samples extend from Flight Time `-10.000 s` to `+0.009 s`:
- **Sample count**: Exactly 1,564 samples (matching workbook formula `AVERAGE(K2:K1565)` in column Z).
- **Mean baseline**: $9.8441\text{ psi}$ ($0.06787\text{ MPa}$).
- **Standard deviation**: $1.0919\text{ psi}$ ($0.00753\text{ MPa}$).
- **Trigger event**:
  - Sample 1563: Flight Time $-0.061\text{ s}$, `Py1 Fire` = $0.0$.
  - Sample 1564: Flight Time $+0.009\text{ s}$, `Py1 Fire` = $1.0$.
  - Transition bracket: $70\text{ ms}$. Extracted $t=0$ anchor is set to the first observed high ($+0.009\text{ s}$).

---

## 3. Burn Profile and Full Tail-Off Audit

The previously extracted CSV (`zerox_high_resolution.csv`) ended at Flight Time 8.153 s ($t_\text{rel} = 8.144\text{ s}$), where pressure was still $4.01\text{ psi}$ gauge ($13.86\text{ psi}$ absolute).

Inspection of the full 10,790 samples confirms the complete tail-off trajectory:

| Flight Time [s] | $t_\text{rel}$ [s] | $p_\text{abs}$ [psi] | $p_\text{gauge}$ [psi] | $p_\text{gauge}$ [MPa] | Thrust [N] | State Description |
|:---|:---|:---|:---|:---|:---|:---|
| 0.349 | 0.340 | 624.19 | 614.35 | 4.2358 | — | Peak initial pressure |
| 4.002 | 3.993 | 425.10 | 415.26 | 2.8631 | — | Mid-burn plateau |
| 6.500 | 6.491 | 185.30 | 175.46 | 1.2097 | — | Web burnout begins |
| 7.001 | 6.992 | 62.78 | 52.94 | 0.3650 | 235.1 | Tail-off steep descent |
| 7.507 | 7.498 | 30.07 | 20.23 | 0.1395 | 53.4 | Tail-off deceleration |
| 8.003 | 7.994 | 17.25 | 7.40 | 0.0510 | -23.3 | Chamber depressurizing |
| 8.153 | 8.144 | 13.86 | 4.01 | 0.0277 | -17.4 | Truncation point of earlier CSV |
| 8.521 | 8.512 | 13.40 | 3.55 | 0.0245 | -32.2 | Residual combustion gas venting |
| 9.008 | 8.999 | 11.56 | 1.72 | 0.0119 | -48.0 | Tail-off complete ($\approx 1.5\sigma$) |
| 10.007 | 9.998 | 10.74 | 0.90 | 0.0062 | -61.1 | Returned to baseline noise floor |
| 15.009 | 15.000 | 10.65 | 0.80 | 0.0055 | -60.2 | Stable post-fire atmospheric level |
| 30.010 | 30.001 | 9.00 | -0.85 | -0.0058 | -54.2 | Long-term thermal drift |

### Tail-Off Window Definition

- **Burnout end / depressurization completion**: Flight Time $9.50\text{ s}$ ($t_\text{rel} = 9.49\text{ s}$).
- Gauge pressure drops below $1.5\text{ psi}$ ($0.010\text{ MPa}$), well within the pre-fire ambient noise envelope ($1\sigma = 1.09\text{ psi}$).
- Post-fire steady baseline ($t_\text{rel} = 10\text{ to }30\text{ s}$) shows zero persistent pressure offset.

---

## 4. Recommendations for Validation Case Registry

1. **Unblock Pressure Linearity**: The logged `pTap 1(psi)` channel is verified as a direct linear conversion of `pTap 1(volts)` with slope $733.01\text{ psi/V}$.
2. **Register Tail-Off Window**: Register the Zerox tail-off window from $t_\text{rel} = [6.50, 9.50]\text{ s}$.
3. **Register Plateau Window**: Register mid-burn plateau window from $t_\text{rel} = [1.50, 5.50]\text{ s}$.
