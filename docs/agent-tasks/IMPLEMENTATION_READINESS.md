# Implementation readiness and execution gates

This document is a planning and review artifact. No source audit, CI run, physical validation or task claim is implied.

## Required gates
1. **Inventory:** read latest main/integration branch, all open PRs, live claims and code; record feature as present/partial/missing/blocked. Do not trust a task title as evidence of missing functionality.
2. **Ownership:** verified create-only atomic claim before task-specific implementation. The GitHub connector's create-or-update branch action is NOT sufficient. If unavailable, fail closed and use manual single-agent assignment.
3. **Bounded handoff:** record task ID, source SHA, target branch, dependency status, exact files likely to change, baseline tests, measurable acceptance, and rollback plan.
4. **Verification:** run unit, integration, negative and regression tests; capture exact commands, environment, versions, inputs, outputs and logs. Mark unrun checks as NOT RUN, never PASS.
5. **Scientific validity:** distinguish code correctness, numerical verification, calibration and independent physical validation. Keep numerical, measurement and model-form uncertainty distinct.
6. **Review:** independent reviewer checks changed behavior, security, reproducibility, dependency contracts, failure modes and evidence; no merge without explicit approval.
7. **Integration:** rebase/retest after dependent PRs land, validate both repositories against pinned versions, and update traceability matrix and changelog.

## Evidence record template
```yaml
task_id: null
repository: israel-shero/srm_1d_Ivans_vibes
base_sha: null
claim_ref: null
claim_verified_atomic: false
implementation_pr: null
owner: null
requirements: []
changed_files: []
unit_tests: {status: NOT_RUN, command: null, log_url: null}
integration_tests: {status: NOT_RUN, command: null, log_url: null}
negative_tests: {status: NOT_RUN, command: null, log_url: null}
regression_tests: {status: NOT_RUN, command: null, log_url: null}
scientific_validation: {status: NOT_RUN, evidence: null}
review: {status: PENDING, reviewer: null}
blockers: []
rollback_plan: null
```

## Initial launch sequence
- Before multi-agent coding: implement and independently race-test exclusive claims; verify cancellation, ambiguity and stale-claim recovery.
- In parallel without implementation claims: inventory source, map dependencies and design test fixtures.
- After claim gate: assign low-conflict bounded PRs, pin contract versions, run independent review and serial integration.

## Known inconsistencies and cautions
- README legacy 'announce' language is not a substitute for atomic claims.
- Markdown tasks in openMotor are not GitHub Issues; status must be reconciled from claims/PRs.
- No task is considered implemented merely because its issue, prompt or PR exists.
