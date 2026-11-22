#!/usr/bin/env python3
"""Publish a durable, human-readable GitHub claim notice after atomic acquisition.

This record is a supervised-pilot audit trail, not cryptographic agent authentication.
On publication failure the existing claim remains locked and must be reconciled.
"""
import argparse
import os
import re
import sys
from claim import ClaimError, TASK, SHA, request

AGENT = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,79}$")
BRANCH = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._/-]{1,159}$")
ACTOR = re.compile(r"^[A-Za-z0-9][A-Za-z0-9-]{0,79}$")

def publish(repo, task, agent, work_branch, scope, actor, run_id, base_sha, ledger, token, api="https://api.github.com"):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+",repo):
        raise ClaimError("Invalid repository")
    if not TASK.fullmatch(task) or not AGENT.fullmatch(agent):
        raise ClaimError("Invalid task or agent identifier")
    if not BRANCH.fullmatch(work_branch) or ".." in work_branch or "//" in work_branch or work_branch.endswith("/"):
        raise ClaimError("Invalid work branch")
    if not ACTOR.fullmatch(actor) or not re.fullmatch(r"[0-9]{1,20}",str(run_id)):
        raise ClaimError("Invalid GitHub actor/run")
    if not SHA.fullmatch(base_sha):
        raise ClaimError("Invalid base SHA")
    if not isinstance(ledger,int) or ledger <= 0:
        raise ClaimError("Invalid ledger comment target")
    if not isinstance(scope,str) or len(scope)>400 or any(ord(c)<32 for c in scope):
        raise ClaimError("Scope must be one line (at most 400 characters)")
    scope = scope.replace("`","'")
    body = (
        "CLAIMED — supervised agent pilot\n\n"
        f"- Task: `{task}`\n"
        f"- Claimed ref: `refs/heads/claims/{task}`\n"
        f"- Agent label (operator supplied, not authenticated): `{agent}`\n"
        f"- GitHub dispatch actor: `{actor}`\n"
        f"- Workflow run (auditable): https://github.com/{repo}/actions/runs/{run_id}\n"
        f"- Intended work branch: `{work_branch}`\n"
        f"- Base commit: `{base_sha}`\n"
        f"- Scope: {scope or '(not provided)'}\n\n"
        "Status: RESERVED; do not code until this record is visible and ownership is reconfirmed. "
        "The claim must not be automatically released, transferred, or merged. "
        "Agent changes require a separate PR, passing CI, and explicit human approval before merge."
    )
    url=api.rstrip("/")+"/repos/"+repo+"/issues/"+str(ledger)+"/comments"
    status, result=request("POST",url,token,{"body":body})
    if status!=201 or not isinstance(result,dict) or not isinstance(result.get("id"),int):
        raise ClaimError("Claim notice not confirmed; STOP, do not code, inspect lock and ledger")
    return {"status":"PUBLISHED","comment_id":result["id"],"comment_url":result.get("html_url")}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("repo","task","agent","work-branch","scope","actor","run-id","base-sha","ledger"):
        p.add_argument("--"+key,required=True)
    args=p.parse_args()
    token=os.getenv("GH_TOKEN") or os.getenv("GITHUB_TOKEN")
    if not token: p.error("GH_TOKEN or GITHUB_TOKEN is required")
    try:
        result=publish(args.repo,args.task,args.agent,args.work_branch,args.scope,args.actor,args.run_id,args.base_sha,int(args.ledger),token)
        print("Claim notice published:",result)
    except (ValueError,ClaimError) as exc:
        print("CLAIM NOTICE FAILED: Lock remains held. STOP and reconcile manually:",exc,file=sys.stderr)
        return 1
    return 0
if __name__=="__main__":
    sys.exit(main())
