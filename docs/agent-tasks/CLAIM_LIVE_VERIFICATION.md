# Claim coordinator live-test evidence (2026-11-22 UTC)

**Repository:** israel-shero/srm_1d_Ivans_vibes  
**Deployed branch:** `main`  
**Merge PR:** https://github.com/israel-shero/srm_1d_Ivans_vibes/pull/38  
**Merge commit:** `af47c491fa945d7d3fffaa327af7dc99c8e2c7d4`  
**Active immutable-claim ruleset:** https://github.com/israel-shero/srm_1d_Ivans_vibes/rules/27825164

## Results
| Test | GitHub Actions evidence | Observed output | Result |
|---|---|---|---|
| Two simultaneous attempts (0 s stagger) | https://github.com/israel-shero/srm_1d_Ivans_vibes/actions/runs/44983142799 | Exactly one `ACQUIRED`, one `CONFLICT`; resulting test ref SHA matches merge commit | PASS |
| Two competing attempts (5 s stagger) | https://github.com/israel-shero/srm_1d_Ivans_vibes/actions/runs/44983279733 | Exactly one `ACQUIRED`, one `CONFLICT`; resulting test ref SHA matches merge commit | PASS |

The tests used two competing threads within the **same Actions process and credential**, calling GitHub's actual create-only `POST /git/refs` API. These tests **did not** exercise independent users, different tokens, permission boundaries, process failures, interrupted network responses or lease recovery. Offline mocked coordinator tests passed in GitHub Actions. Provenance CI passed.

## Durable test refs
- `refs/heads/claims/TEST-CLAIM-PILOT-20261121-A`
- `refs/heads/claims/TEST-CLAIM-PILOT-20261121-B`

Both were independently fetched via GitHub REST after the workflows completed, and both resolve to `af47c491fa945d7d3fffaa327af7dc99c8e2c7d4`. These are **testing-only** references. Do not reuse them as production claims, steal them, or delete them automatically. Any cleanup requires authorized administrator action and explicit approval.

## Remaining rollout gates
- [ ] Design verified task owner identity, unique generation/fencing token, heartbeat and recovery with safe stale-claim reconciliation.
- [ ] Implement mandatory PR check that confirms the actor and generation owns the exact task claim, not merely that the ref exists.
- [ ] Run independent-process / independently authenticated concurrent acquisition tests.
- [ ] Run deliberate network failure, expired-owner and recovery simulations.
- [ ] Enable appropriate default-branch required checks and limited bypass policies for PR acceptance.
- [ ] Validate cross-repository dependencies, task status reconciliation and abort/cancellation behavior.

**Current readiness:** Atomic create-ref acquisition is **verified for the two controlled race cases above**. Protected claim refs are configured. Full unattended parallel task ownership and safe recovery are **NOT VERIFIED**. Do not treat this as permission to run autonomous agents without supervision.
