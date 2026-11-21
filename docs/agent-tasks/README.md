# SRM 1D — agent task queue

## Mandatory before starting any task
Read [CLAIM_PROTOCOL.md](CLAIM_PROTOCOL.md). **Claim the task atomically on GitHub before any task-specific work.** If a lock already exists or create-only behavior cannot be guaranteed, do not start; choose another eligible task. Issue labels/comments alone are not sufficient for concurrency safety.

**Implementation repository:** https://github.com/israel-shero/srm_1d_Ivans_vibes
**Separate GUI implementation repository:** https://github.com/israel-shero/openMotor_Ivan (do not write there from these tasks).

## Instructions for Codex/agents
1. Fetch current branches, issues, PRs, `AGENTS.md` and relevant tests. Do not assume this index is current.
2. Select one dependency-ready issue; check that its work has not already been implemented, and announce issue number, scope, proposed branch and PR base.
3. Follow the issue's detailed specification and Codex handoff. Make bounded commits with regression tests, numerical/physical evidence, changelog and exact command/results.
4. Prefer `main` for solver-native changes and `integration/openmotor-packaging` for adapter/plugin-specific changes, only after verifying branch state. Never change defaults or claim physical validation without evidence.
5. Coordinate interface contracts with openMotor_Ivan through references and separate PRs. Never silently edit both repos as part of one issue.

## Existing scoped tasks
| Issue | Area | Dependencies / readiness |
|---|---|---|
| [#2](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/2) | Adaptive timestep controller, phased baseline/policy/rollback | Start with baseline PR A; preserve legacy defaults |
| [#3](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/3) | Continuous axial profiles, conservative mesh | Start with geometry audit/contract |
| [#4](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/4) | Remaining propellant mass, CG, inertia | Initial audit can run; coordinate geometry #3 |
| [#7](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/7) | T07 candidate decision | Blocked until independent fit/holdout evidence; #8 feeds it |
| [#8](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/8) | Case source and calibration audit | Ready; parallel with #2/#3 |
| [#9](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/9) | Full-burn spatial convergence | Audit ready; coordinate baseline work #2 |
| [#10](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/10) | Adapter/serialization regression contract | Audit ready; coordinate #3 and openMotor_Ivan M0-0B/0C |

## Parallel work
- Parallel low-conflict starts: #2 baseline audit, #3 geometry contract, #8 source audit, #10 adapter contract.
- #4 mass-properties contract can be researched concurrently, but geometry-sensitive implementation must reconcile #3.
- #9 can review verification design concurrently, but expensive runs should use frozen #2 baselines.
- #7 is a blocked decision gate, not permission to tune a closure or assign unverified fit/holdout roles.

**Historical checkpoint:** The 2026-11-20 project checkpoint records T02 source reconciliation, T06 common-state/coupled comparison and T10 GUI tests as completed, while full-burn convergence and empirical validation remain incomplete. Verify all statuses against current source and evidence.
