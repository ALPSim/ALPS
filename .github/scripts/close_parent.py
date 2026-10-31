#!/usr/bin/env python3
"""Close the parent of issue ISSUE once all of its sub-issues are closed.

Sub-issues close when their PRs merge, so the parent closes when the last of
its PRs is done. Grandparents follow, because closing the parent with this
workflow's token triggers the workflow again for it.
"""
import json
import os
import subprocess

OWNER, REPO = "ALPSim", "ALPS"


def gh(*args, body=None):
    out = subprocess.run(["gh", *args], input=None if body is None else json.dumps(body),
                         capture_output=True, text=True, check=True).stdout
    return json.loads(out) if out.strip() else None


query = """query($o: String!, $r: String!, $n: Int!) { repository(owner: $o, name: $r) { issue(number: $n) {
    parent { number state repository { nameWithOwner } subIssuesSummary { total completed } } } } }"""
parent = gh("api", "graphql", "--input", "-", body={"query": query, "variables": {
    "o": OWNER, "r": REPO, "n": int(os.environ["ISSUE"])}})["data"]["repository"]["issue"]["parent"]
if not parent or parent["state"] != "OPEN" or parent["repository"]["nameWithOwner"] != f"{OWNER}/{REPO}":
    print("no open parent in this repository")
elif parent["subIssuesSummary"]["completed"] < parent["subIssuesSummary"]["total"]:
    print(f"#{parent['number']}: {parent['subIssuesSummary']['completed']}"
          f"/{parent['subIssuesSummary']['total']} sub-issues closed")
else:
    gh("api", f"repos/{OWNER}/{REPO}/issues/{parent['number']}/comments", "--input", "-",
       body={"body": "All sub-issues are closed, so this issue is done."})
    gh("api", "-X", "PATCH", f"repos/{OWNER}/{REPO}/issues/{parent['number']}", "--input", "-",
       body={"state": "closed", "state_reason": "completed"})
    print(f"closed #{parent['number']}")
