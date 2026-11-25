# Agent development instructions — israel-shero/srm_1d_Ivans_vibes

## Select and claim your own task

Choose the highest-priority eligible task in **this repository**, respecting live dependencies, claims, existing PRs and actual implemented code. If the user specifically assigns a task (for example INF-01), verify and prioritize that task rather than choosing a different one. Review live GitHub issues first; the more extensive docs/agent-task-queue branch may include unmerged planning material.

Read [agent self-assignment](docs/agent-tasks/AGENT_SELF_ASSIGNMENT.md), [supervised pilot policy](docs/agent-tasks/SUPERVISED_PILOT.md) and [claim protocol](docs/agent-tasks/CLAIM_PROTOCOL.md).

**You are responsible for claiming your task.** Using your own authorized **non-admin** GitHub credential with Actions: write, dispatch `claim-request.yml` yourself on `main` with task ID, agent ID, work branch and scope. Do not ask the human operator to press the Run workflow button. A dispatch is not a claim: verify completed success, the exact `claims/TASK-ID` ref and the matching published owner record before touching task code. If missing permissions, stop and report the one-time access needed.

For an eligible unclaimed task, implement only on `work/TASK-ID/AGENT-ID` (or the exact implementation branch recorded in your successful claim); test it and submit a scoped PR to `dev`. Preserve numerical correctness and compatibility. No task hopping, claim release or reassignment without a fresh authorized process.

**Do not bypass branch protection, approve your own PR to satisfy review, use administrator bypass, release/delete claims, or change GitHub rulesets.** After opening the PR, report evidence. An authorized integrator may merge into `dev` under the permanent branch policy below; `main` promotions require owner approval.

## Permanent branch policy

- `dev` is the primary development and integration branch. Start task branches from `dev` and target ordinary development pull requests to `dev`. Tested, scoped changes may be merged into `dev` without a separate owner approval for each merge, subject to existing GitHub permissions, required checks, review rules and task-claim controls. Never bypass protections or approve your own PR to satisfy required review.
- `main` is the stable milestone branch. Promote coherent, tested integration milestones from `dev` to `main` by a deliberate pull request; require explicit owner approval. Do not make routine development commits to `main`.
- Preserve both branches and their history. Do not force-push, delete either branch, or alter permissions or rulesets without explicit authorization. Existing trusted workflow dispatch refs are separate from development PR targets; do not change privileged workflow execution refs without a security review.
