#!/usr/bin/env python3
"""Live GitHub atomic ref claim race test. Intentionally leaves a TEST ref."""
import argparse
import concurrent.futures
import os
import re
import threading
import time
import sys
from claim import ClaimError, run

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",required=True)
    parser.add_argument("--task",required=True,help="Unique TEST-CLAIM-* ID; never reuse")
    parser.add_argument("--base-sha",required=True)
    parser.add_argument("--stagger-seconds",type=float,default=0)
    args=parser.parse_args()
    if not re.fullmatch(r"TEST-CLAIM-[A-Za-z0-9._-]{1,55}",args.task):
        parser.error("Only TEST-CLAIM-* IDs are allowed")
    if not 0 <= args.stagger_seconds <= 10:
        parser.error("stagger must be between 0 and 10 seconds")
    token=os.getenv("GH_TOKEN") or os.getenv("GITHUB_TOKEN")
    if not token:
        parser.error("GH_TOKEN/GITHUB_TOKEN not set")
    gate=threading.Barrier(2)
    def attempt(i):
        gate.wait(timeout=10)
        if i==1:
            time.sleep(args.stagger_seconds)
        try:
            result=run(args.repo,args.task,args.base_sha,token,"https://api.github.com","acquire")
            return {"attempt":i,"result":result["status"]}
        except ClaimError as exc:
            msg=str(exc)
            if "GitHub HTTP 422" in msg or "GitHub HTTP 409" in msg:
                return {"attempt":i,"result":"CONFLICT"}
            return {"attempt":i,"result":"UNEXPECTED_ERROR","reason":msg[:240]}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures=[pool.submit(attempt,i) for i in range(2)]
        results=[f.result(timeout=60) for f in futures]
    print("Claim attempt results:",results)
    if sorted(x["result"] for x in results)!=["ACQUIRED","CONFLICT"]:
        print("FAIL: Expected exactly one ACQUIRED and one CONFLICT. STOP.",file=sys.stderr)
        return 1
    try:
        inspected=run(args.repo,args.task,None,token,"https://api.github.com","inspect")
    except ClaimError as exc:
        print("FAIL: Cannot inspect created ref:",exc,file=sys.stderr)
        return 1
    if inspected.get("sha")!=args.base_sha:
        print("FAIL: Ref SHA changed unexpectedly; stop.",file=sys.stderr)
        return 1
    print("PASS: exactly one winner; claim ref remains for authorized cleanup.",inspected)
    return 0
if __name__=="__main__":
    sys.exit(main())
