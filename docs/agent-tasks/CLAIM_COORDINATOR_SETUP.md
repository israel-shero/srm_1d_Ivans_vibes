# Claim coordination — deployment, tests and manual completion

## What has been prepared
- `tools/claim-coordinator/claim.py`: stdlib-only **create-only** GitHub reference acquisition and inspection, fail-closed on errors.
- `tools/claim-coordinator/test_claim.py`: offline unit tests for success, conflict, ambiguity, invalid identifiers and malformed responses.
- This guide describes the operations that cannot be completed safely through the current connector.

**Not deployed:** No authenticated create-ref request, live race test, workflow execution, required status check, branch protection change, claim lease service or PR merge has been performed. This prototype is not a complete multi-agent scheduler. Do not launch unattended parallel agents until the remaining gates pass.

## Step 1 — review and merge the documentation PR
Review the task queue and coordinator code in the `docs/agent-task-queue` branch. The branch contains planning artifacts and is not yet necessarily the default branch. Merge only after human review. Do not copy this code into a privileged production workflow without reviewing permissions and untrusted inputs.

## Step 2 — choose the deployment identity
Use a GitHub App installation token with narrowly scoped **Contents: read/write** on only the intended repository. A short-lived fine-grained token is an alternative for local manual operation. Never commit credentials or print them to logs. For GitHub Actions, use the ephemeral `GITHUB_TOKEN` with `permissions: contents: write` only in a trusted, manually triggered workflow; fork PRs and untrusted issue text must not execute privileged coordinator code. Ensure branch/ref rules allow creating `claims/*` and deny unauthorized modification/deletion.

## Step 3 — run offline tests
From the repository root:
```bash
python3 -m unittest discover -s tools/claim-coordinator -p 'test_*.py' -v
```
Record Python version, commands and logs. These tests mock HTTP and do **not** prove GitHub concurrency behavior.

## Step 4 — run a controlled live acquisition
Obtain the intended base commit SHA from GitHub, not from a stale checkout:
```bash
git rev-parse origin/main
```
For the GUI repository, use the verified integration base instead of `main` if that remains the approved target.

With `GH_TOKEN` in the environment:
```bash
python3 tools/claim-coordinator/claim.py acquire --repo OWNER/REPO --task TEST-CLAIM-001 --base-sha FULL_40_CHARACTER_SHA
python3 tools/claim-coordinator/claim.py inspect --repo OWNER/REPO --task TEST-CLAIM-001
```
**Do not reuse a production task ID for testing.** Confirm a second acquisition of the same task fails with a conflict. Because claims are persistent, arrange explicit review and authorization for cleaning up the test ref. This helper intentionally contains no delete operation.

## Step 5 — concurrency test
Using two independently authenticated processes, synchronize simultaneous `acquire` attempts against the **same fresh** `TEST-CLAIM-002` ref. Exactly one must receive a verified HTTP 201 response; the other must fail. Repeat across different clients and a 5-second stagger. Test network timeout after a possible successful POST: both callers must stop and inspect rather than retry blindly. Log redacted results, never tokens.

## Step 6 — ownership publication and fencing (NOT YET IMPLEMENTED)
After successful acquisition, publish an issue comment or designated task record containing task ID, actor/run ID, claim ref, base SHA, work branch and UTC time. Verify that this publication succeeds before work. **The ref alone does not encode agent identity**; it only prevents duplicate acquisitions. Build a signed/authorized claim registry and generation/fencing checks before enabling automated lease expiry or reassignment. Existing claims must never be deleted merely because a heartbeat is old. Require explicit coordinator reconciliation of stale claims and open PRs.

## Step 7 — enforce at PR and branch protection (NOT YET IMPLEMENTED)
Create a trusted required status check that validates task ID, claim ref, authorized owner/generation, approved base and task-specific test evidence. Make it a required check through Settings → Branches / Rulesets, restrict bypass permissions, require review and passing CI. A check that only tests whether a claim ref exists is **insufficient** to prove PR ownership. Check API token permissions and rulesets with a repository administrator. Do not auto-merge.

## Step 8 — operational rollout
1. Start with one manual claim and one small PR.
2. Add structured registry and heartbeat with fencing, tested crash/recovery behavior and audit trail.
3. Integrate task dependency readiness and change-file conflict prediction.
4. Enable a small parallel pilot after live concurrency and PR enforcement pass.
5. Record evidence in `docs/agent-tasks/VALIDATION_MATRIX.md` and update `REQUIREMENTS_COVERAGE.json` with actual links, not inferred PASS statuses.

## Recovery playbook
- **HTTP 422/conflict:** another claim or invalid request; stop, inspect ref, never overwrite.
- **Timeout/5xx/unknown outcome:** stop, inspect ref; do not automatically retry POST.
- **Claim created but notice failed:** stop; report orphan claim for authorized reconciliation.
- **Agent vanished:** retain lock, inspect heartbeat, PR, work branch and ownership; coordinator-approved release only.
- **Unexpected ref moved/deleted:** stop all work for that task; audit permissions and event logs.
- **Cross-repo task:** claim each repository separately and pin compatible contract versions.

## Outstanding manual actions checklist
- [ ] Review and merge documentation PR
- [ ] Configure least-privilege GitHub App or trusted Actions token
- [ ] Run and archive offline tests
- [ ] Perform and document real GitHub atomic race tests
- [ ] Build ownership registry, fencing and stale-claim recovery
- [ ] Add required PR ownership check and repository ruleset
- [ ] Test pilot and document results
- [ ] Approve production multi-agent rollout

**Security boundary:** Never treat branch-name existence, issue labels, comments or heartbeat timestamps alone as verified ownership. This is a safe initial acquisition primitive, not complete coordination.
