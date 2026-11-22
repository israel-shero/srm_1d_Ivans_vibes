#!/usr/bin/env python3
"""Fail-closed GitHub exclusive task claim helper (stdlib only).

Creates fixed refs/heads/claims/<TASK> with GitHub REST POST /git/refs.
This is create-only; never uses create-or-update branch operations.
Does not release, renew or steal locks. A claim is persistent until reviewed.
"""
import argparse
import json
import os
import re
import sys
import urllib.error
import urllib.parse
import urllib.request

TASK = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$")
SHA = re.compile(r"^[a-fA-F0-9]{40}$")

class ClaimError(Exception):
    pass

def request(method, url, token, data=None):
    payload = None if data is None else json.dumps(data).encode("utf-8")
    req = urllib.request.Request(url, data=payload, method=method, headers={
        "Authorization": "Bearer " + token,
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
        "User-Agent": "srm-agent-claim/0.1",
        "Content-Type": "application/json",
    })
    try:
        with urllib.request.urlopen(req, timeout=20) as response:
            return response.status, json.load(response)
    except urllib.error.HTTPError as exc:
        raise ClaimError("GitHub HTTP %d: %s" % (exc.code, exc.read(600).decode("utf-8", "replace"))) from exc
    except (OSError, ValueError) as exc:
        raise ClaimError("Ambiguous network/API outcome; STOP and inspect claim ref: %s" % exc) from exc

def claim_url(api, repo, task):
    owner, name = repo.split("/")
    prefix = api.rstrip("/") + "/repos/" + urllib.parse.quote(owner) + "/" + urllib.parse.quote(name)
    ref = "refs/heads/claims/" + task
    lookup = prefix + "/git/ref/" + urllib.parse.quote("heads/claims/" + task, safe="/")
    return prefix, ref, lookup

def run(repo, task, base, token, api, action):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repo):
        raise ClaimError("Invalid owner/repository")
    if not TASK.fullmatch(task):
        raise ClaimError("Invalid task ID")
    prefix, ref, lookup = claim_url(api, repo, task)
    if action == "inspect":
        _, obj = request("GET", lookup, token)
        return {"status":"EXISTS","ref":obj.get("ref"),"sha":obj.get("object",{}).get("sha")}
    if not base or not SHA.fullmatch(base):
        raise ClaimError("Acquire requires explicit verified 40-character base SHA")
    # No preliminary existence check can substitute for the atomic POST.
    status, obj = request("POST", prefix + "/git/refs", token, {"ref":ref,"sha":base})
    if status != 201 or obj.get("ref") != ref or obj.get("object",{}).get("sha") != base:
        raise ClaimError("Unexpected create response: STOP, inspect ref manually")
    return {"status":"ACQUIRED","ref":ref,"sha":base,"next":"Publish ownership notice and verify claim before work; no automatic release"}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("acquire","inspect"))
    p.add_argument("--repo", required=True)
    p.add_argument("--task", required=True)
    p.add_argument("--base-sha")
    p.add_argument("--api", default="https://api.github.com")
    args = p.parse_args()
    token = os.getenv("GH_TOKEN") or os.getenv("GITHUB_TOKEN")
    if not token:
        p.error("Set GH_TOKEN or GITHUB_TOKEN with Contents:write permission")
    try:
        print(json.dumps(run(args.repo,args.task,args.base_sha,token,args.api,args.action),indent=2))
    except ClaimError as exc:
        print("CLAIM NOT GRANTED: " + str(exc),file=sys.stderr)
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
