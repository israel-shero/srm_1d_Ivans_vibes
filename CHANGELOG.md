# Changelog

## Unreleased

### Fixed

- Preserve the pyrogen gas mass fraction when saving and loading igniter blocks. Older blocks that omit the field retain the established default of 1.0.
- Reject limiting-case resumes when source content, dependency versions, effective settings, or completed-run identities differ. Legacy manifests without these identities require a fresh run.
- Require single-point numerical verification to reach its configured duration; retain diagnostics when a run truncates or fails.
- Update the workspace lock to the compatible openMotor checkout used by the integration branch.

### Added

- Implement solver config upgrade (`T_ignition`, `kappa`, `roughness`) and legacy uncustomized BPNV pyrogen block upgrade hooks for openMotor plugin integration.
- Explicit `--single-erosive-closure` selection for numerical verification, with source, runtime, and generated-geometry checksums in each point manifest. Ma remains the default.
- Complete-duration Mukunda–Paul grid-sensitivity records and bounded primary-source evidence. These records do not establish formal convergence or experimental validation.
- Ma matched-grid coupled comparison (106/206 cells, CFL 0.075, Chunc motor, full 3 s): all selected global metrics within provisional 2% screen. Ma shows worse startup grid sensitivity than Mukunda-Paul (peak pressure 1.05% vs 0.027%, peak time −5.61% vs +0.06%). This is a numerical sensitivity comparison, not a validation claim.
- Chunc startup spatial extension (100→200→400→800→1600 cells, CFL 0.075, t_max 0.05 s): peak pressure monotone-converging (all successive changes below 2% screen); peak time and Mach finest-interval below screen but not all; minimum pressure censored (floor activated from 400 cells). Formal convergence not established.
- Zerox pressure transducer calibration audit: linear transfer function confirmed (slope 733.011 psi/V, intercept −721.016 psi), plateau window [1.50, 5.50] s and tail-off window [6.50, 9.50] s formally registered in validation case registry under protocol v0.1. Audit document at `docs/v0_8_0/ZEROX_TRANSDUCER_CALIBRATION_AUDIT.md`.
- openMotor GUI unit tests: all 29 tests now green (configEditor, propellant, igniter) after igniter upgrade hook and string-property blanking fix.
