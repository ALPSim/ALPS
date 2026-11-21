#!/usr/bin/env python3
"""Assign each issue on the board to the person whose name its card is under.

For every open ALPSim/ALPS issue on the project board
(github.com/orgs/ALPSim/projects/5) whose Owner is a person rather than
Bucket, that person is assigned to the issue, unless they were unassigned
after the card was put under their name: someone did that on purpose.

Moving a card from one name to another unassigns the previous person if a
workflow had assigned them. People assigned by hand stay assigned, and so does
everyone on an issue whose card goes back to Bucket.

BOARD_TOKEN reads the board; GH_TOKEN edits the issues. With DRY_RUN=1 it only
prints what it would do.
"""
import json
import os
import subprocess

ORG, PROJECT = "ALPSim", 5
REPO = "ALPSim/ALPS"
BUCKET = "Bucket"
DRY_RUN = os.environ.get("DRY_RUN") == "1"
BOT = "github-actions"  # actor of assignments made with GITHUB_TOKEN

ITEMS = """
query($org: String!, $number: Int!, $after: String) {
  organization(login: $org) {
    projectV2(number: $number) {
      owner: field(name: "Owner") { ... on ProjectV2SingleSelectField { options { name } } }
      items(first: 100, after: $after) {
        pageInfo { hasNextPage endCursor }
        nodes {
          owner: fieldValueByName(name: "Owner") { ... on ProjectV2ItemFieldSingleSelectValue { name updatedAt } }
          content {
            ... on Issue {
              number state
              repository { nameWithOwner }
              assignees(first: 20) { nodes { login } }
              timelineItems(itemTypes: [ASSIGNED_EVENT, UNASSIGNED_EVENT], last: 50) {
                nodes {
                  __typename
                  ... on AssignedEvent { createdAt actor { login } assignee { ... on User { login } } }
                  ... on UnassignedEvent { createdAt assignee { ... on User { login } } }
                }
              }
            }
          }
        }
      }
    }
  }
}
"""


def gh(*args, token, body=None):
    """gh's output; `gh issue edit` prints the issue's URL, not JSON."""
    return subprocess.run(["gh", *args], input=None if body is None else json.dumps(body),
                          capture_output=True, text=True, check=True,
                          env={**os.environ, "GH_TOKEN": token}).stdout


def load_board():
    names, items, after = None, [], None
    while True:
        page = json.loads(gh("api", "graphql", "--input", "-", token=os.environ["BOARD_TOKEN"],
                             body={"query": ITEMS, "variables": {"org": ORG, "number": PROJECT, "after": after}}))
        project = page["data"]["organization"]["projectV2"]
        names = names or [o["name"] for o in project["owner"]["options"]]
        items += project["items"]["nodes"]
        if not project["items"]["pageInfo"]["hasNextPage"]:
            return names, items
        after = project["items"]["pageInfo"]["endCursor"]


def last_events(issue):
    """Latest AssignedEvent and UnassignedEvent per login (ISO times sort as strings)."""
    assigned, unassigned = {}, {}
    for e in issue["timelineItems"]["nodes"]:
        login = ((e.get("assignee") or {}).get("login") or "").lower()
        if not login:
            continue
        if e["__typename"] == "AssignedEvent":
            assigned[login] = e
        else:
            unassigned[login] = e
    return assigned, unassigned


def main():
    names, items = load_board()
    people = {n.lower(): n for n in names if n != BUCKET}
    failed = False
    for item in items:
        issue = item["content"] or {}
        owner = (item["owner"] or {}).get("name")
        moved = (item["owner"] or {}).get("updatedAt") or ""
        if (issue.get("repository") or {}).get("nameWithOwner") != REPO or issue.get("state") != "OPEN":
            continue
        if not owner or owner.lower() not in people:
            continue
        assigned = [a["login"] for a in issue["assignees"]["nodes"]]
        by, off = last_events(issue)
        add = []
        if owner.lower() not in (a.lower() for a in assigned):
            if off.get(owner.lower(), {}).get("createdAt", "") < moved:
                add = [owner]
        # Only undo assignments made by a workflow, and only after the card moved.
        remove = [a for a in assigned
                  if a.lower() in people and a.lower() != owner.lower()
                  and (by.get(a.lower(), {}).get("actor") or {}).get("login") == BOT
                  and by[a.lower()]["createdAt"] < moved]
        if not add and not remove:
            continue
        print(f"#{issue['number']}: add {add or '-'}, remove {remove or '-'}")
        if DRY_RUN:
            continue
        args = ["issue", "edit", str(issue["number"]), "--repo", REPO]
        for login in add:
            args += ["--add-assignee", login]
        for login in remove:
            args += ["--remove-assignee", login]
        try:
            gh(*args, token=os.environ["GH_TOKEN"])
        except subprocess.CalledProcessError as e:
            # Someone without access to the repository can't be assigned; carry on with the rest.
            print(f"#{issue['number']}: {e.stderr.strip()}")
            failed = True
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
