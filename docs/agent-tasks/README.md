# SRM 1D — agent task queue

## Mandatory before starting any task
Read [CLAIM_PROTOCOL.md](CLAIM_PROTOCOL.md). **Claim the task atomically on GitHub before any task-specific work.** If a lock already exists or create-only behavior cannot be guaranteed, do not start; choose another eligible task. Issue labels/comments alone are not sufficient for concurrency safety.

**Implementation repository:** https://github.com/israel-shero/srm_1d_Ivans_vibes
**Separate GUI implementation repository:** https://github.com/israel-shero/openMotor_Ivan (do not write there from these tasks).

## Instructions for Codex/agents
1. Fetch current branches, issues, PRs, `AGENTS.md` and relevant tests. Do not assume this index is current.
2. Select one dependency-ready issue; check that its work has not already been implemented, then acquire a verified create-only exclusive claim per CLAIM_PROTOCOL.md before announcing issue number, scope, proposed branch and PR base.
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

## Reliability and coverage backlog (new)
| Issue | Work | Dependency |
|---|---|---|
| [#12](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/12) | INF-01 atomic claims and safe lease recovery | Highest priority for concurrent agents |
| [#13](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/13) | INF-02 registry and scheduler | #12 for live claims; schema can begin now |
| [#14](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/14) | INF-03 cross-repo contract CI | Coordinate #10 and openMotor INF-03 |
| [#15](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/15) | INF-04 independent verification | Can begin now |
| [#16](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/16) | INF-05 performance baselines | Can begin now |
| [#17](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/17) | QA-01 property/metamorphic tests | Coordinate #2/#3/#4 |

**Recommended parallel starts:** #12 design/adversarial tests, #14 contract inventory, #15 verifier report schema, #17 invariant test design. Do not dispatch implementation without a verified exclusive claim.

## Uncertainty quantification and identification workstream
| Issue | Task | Dependencies |
|---|---|---|
| [#18](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/18) | UQ-01 parameter schema, distributions, provenance | Start here |
| [#19](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/19) | UQ-02 batch runner, failure accounting and checkpointing | UQ-01 schema; runner scaffold can start separately |
| [#20](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/20) | UQ-03 Monte Carlo/LHS/quasi-Monte Carlo | UQ-01, UQ-02 |
| [#21](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/21) | UQ-04 local/global sensitivity | UQ-01, UQ-02 |
| [#22](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/22) | UQ-05 identification and identifiability | UQ-01, UQ-02, evidence #8 |
| [#23](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/23) | UQ-06 uncertainty and numerical error reporting | UQ-02, UQ-03, UQ-04 |

**Architecture:** SRM 1D owns all UQ computation; openMotor_Ivan owns UI-08. First validate on analytic benchmark functions, then pinned solver fixtures. Never silently discard failed samples or treat fitted parameters as independently validated.

## Full-conversation coverage expansion
| Issue | Scope |
|---|---|
| [#24](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/24) | PHY-01 conservation/boundaries |
| [#25](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/25) | PHY-02 thermal/ignition/erosive verification |
| [#26](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/26) | VAL-01 experimental evidence and discrepancy |
| [#27](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/27) | DATA-01 drawing evidence inventory |
| [#28](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/28) | DATA-02 geometry fixture completeness |
| [#29](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/29) | ARCH-01 API/release/diagnostics |
| [#30](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/30) | AGENT-01 merge queue/checkpoints |
| [#31](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/31) | RESEARCH-01 paper-to-code traceability |
| [#32](https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/32) | QA-03 numerical/release acceptance |

These are bounded follow-ons; #8/#9/#10 and existing UQ/INF tasks remain authoritative for their overlapping scopes. Check live issue and PR status before claiming.

## Implementation readiness and evidence gates
- [Implementation readiness and execution handoff](IMPLEMENTATION_READINESS.md)
- [Validation evidence matrix](VALIDATION_MATRIX.md)
- [Machine-readable initial requirements coverage](REQUIREMENTS_COVERAGE.json)

**Important:** These are planning artifacts, not proof of code coverage or passing tests. No multi-agent task execution until verified atomic claims are available; task announcement alone does not reserve work.

## Claim coordinator implementation and remaining manual deployment
- [Claim coordinator deployment, live tests and required manual steps](CLAIM_COORDINATOR_SETUP.md)
- Prototype acquisition CLI: [`tools/claim-coordinator/claim.py`](../../tools/claim-coordinator/claim.py)
- Offline tests: [`tools/claim-coordinator/test_claim.py`](../../tools/claim-coordinator/test_claim.py)

**Status:** code committed but not executed in this session; no live atomic race test, identity/fencing registry, branch ruleset or PR ownership enforcement yet. Do not run unattended multi-agent tasks until manual checklist passes.
