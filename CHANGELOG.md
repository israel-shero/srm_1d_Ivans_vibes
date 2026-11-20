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
