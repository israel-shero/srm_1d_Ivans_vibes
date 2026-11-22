# Agent-initiated task selection and exclusive claiming

**Default:** The agent chooses a task, claims it, and develops it without seeking approval to start. Human approval is required for **merging the PR**, modifying shared-branch protections, releasing/reassigning stale claims, or other privileged changes.

## One-time identity and permissions setup

Provide each agent a GitHub identity/credential with access to the **specific repository only**. To dispatch `claim-request.yml`, GitHub requires **Actions: write**. The workflow itself creates the claim and records it using its ephemeral `GITHUB_TOKEN`; the agent does **not** need unrestricted direct access to `POST /git/refs` merely to claim. For work-branch commits and PRs, the agent separately needs the permissions its development environment uses (typically Contents: write and Pull requests: write). It must have **no repository administrator identity, ruleset bypass, or permission to merge into protected `main`/`staging`**. Do not give agents a human owner's administrator PAT: an owner account may carry admin bypass even with narrowly scoped credentials.

A GitHub App installation or distinct non-admin bot identity is preferable; its installation token must be provisioned through the agent platform and available to `gh`. GitHub's Actions: write is a repository-level capability, **not limited to only this one workflow**. Only provision it to agents you trust to run appropriate actions; if that is unacceptable, use a reviewed comment-to-claim workflow with narrower trigger permissions instead. Credentials are configured once, not for every task.

## Agent workflow (NO manual claim button)

1. Read the current repository's `AGENTS.md`, this guide, `SUPERVISED_PILOT.md`, and `CLAIM_PROTOCOL.md`. Inspect live issues, open PRs, existing `claims/*`, dependencies and the task index (some older index copies reside on separate documentation branches).
2. Choose the **highest-priority actionable, dependency-ready, unclaimed** task in this repository. If the user gave an explicit task (such as INF-01), treat that as the intended task and verify eligibility. A task with existing code should be narrowed to documented remaining work rather than rebuilt.
3. Pick a distinct agent label, e.g. `codex-A`, and `work/TASK-ID/codex-A` implementation branch. Select a brief single-line scope (max 400 characters). Run, using **your own authorized GitHub CLI identity**:

```bash
gh workflow run claim-request.yml \
  --repo OWNER/REPOSITORY \
  --ref main-or-staging \
  --field task_id=TASK-ID \
  --field agent_id=codex-A \
  --field work_branch=work/TASK-ID/codex-A \
  --field scope='Bounded task scope'
```

4. **Do not treat a successful dispatch as a successful claim.** Find the corresponding `Request exclusive task claim` workflow run. Check that it completed successfully and that `refs/heads/claims/TASK-ID` exists, then independently find the ownership notice created by that same workflow run: SRM 1D in issue #12 comments; openMotor as a commit comment on the run's base SHA. Verify **all** of task, agent label, authenticated dispatch actor, work branch, scope, workflow URL and base SHA. Never infer ownership from a label or mere branch existence.
5. If the workflow returned a **confirmed conflict** and the claim is owned by someone else, choose the next eligible task and begin again (at most one successful active task per agent). If the dispatch/API status is **ambiguous**, the claim exists but no notice is published, or the notice does not match, **STOP** and request operator reconciliation. Never force-create, retry the same task blindly, or seize an expired claim.
6. Once success is independently verified, implement only the claimed scope on the recorded work branch. Use existing requirements and tests. Open a PR against the designated base branch with CI/test evidence and the claim-run link. **Stop after PR submission**; do not merge it, release the claim, or start another task until the supervised pilot owner approves a new work cycle.

### How to find the record

- SRM 1D: <https://github.com/israel-shero/srm_1d_Ivans_vibes/issues/12>; workflow branch `main`.
- openMotor: use GitHub's **commit comments** on the workflow's verified `github.sha`; workflow branch `staging`. The repository has Issues disabled, so do not depend on issue comments.
- Both repositories: `.github/workflows/claim-request.yml` is the only approved create-only task-claim entrypoint for this pilot.

## Stop conditions

Missing `gh` authentication or Actions: write, an administrator token, a claimed/blocked task, an unverified workflow, any attempted direct push to `main`/`staging`, missing CI, overlapping work, or any request to bypass review are reasons to stop and report. Do not ask the human to press the claim button as the standard fallback: instead, report exactly what **one-time credential capability** is missing.

**Security boundary:** Atomic create-ref prevents competing claims. Agent labels are not authenticated identities; `github.actor` identifies the authenticated dispatch principal. Self-selection and self-claiming are allowed, but **self-merging and automatic claim reassignment remain prohibited**.
