#!/usr/bin/env python3
"""Open a draft PR for every sub-issue of an issue in "Generate PRs".

For each ALPSim/ALPS issue on the project board (github.com/orgs/ALPSim/projects/5)
whose Status is "Generate PRs", or only the issue numbered ISSUE if that is set,
every open sub-issue without a linked PR gets
- a branch issue-<number>, linked to the sub-issue, holding one empty commit,
- a draft PR from that branch whose body closes the sub-issue,
- a card on the board in "Generate PRs" with the parent's Owner.
The parent gets one comment listing each sub-issue's PRs, kept up to date, so
clicking its card shows every PR it branched into.

Issues in "Make Subissues" get a "Generate PRs for the sub-issues" checkbox in
their description; ticking it runs this for that issue at once, and the box is
unticked again afterwards.

Running it again does nothing for sub-issues that already have a PR or branch.
At most MAX_PRS PRs are opened per run. With DRY_RUN=1 it only prints what it
would do.
"""
import json
import os
import subprocess

ORG, PROJECT = "ALPSim", 5
OWNER, REPO = "ALPSim", "ALPS"
# Status options, by ID so renaming the columns doesn't break this.
STATUS = "1d217230"  # Generated PRs
MAKE_SUBISSUES = "47fc9ee4"
BOX = "Generate PRs for the sub-issues"
ISSUE = int(os.environ.get("ISSUE") or 0)
MAX_PRS = 10
DRY_RUN = os.environ.get("DRY_RUN") == "1"
MARKER = "<!-- subissue-prs -->"

ITEMS = """
query($org: String!, $number: Int!, $after: String) {
  organization(login: $org) {
    projectV2(number: $number) {
      id
      status: field(name: "Status") { ... on ProjectV2SingleSelectField { id options { id name } } }
      owner: field(name: "Owner") { ... on ProjectV2SingleSelectField { id options { id name } } }
      items(first: 100, after: $after) {
        pageInfo { hasNextPage endCursor }
        nodes {
          status: fieldValueByName(name: "Status") { ... on ProjectV2ItemFieldSingleSelectValue { optionId } }
          owner: fieldValueByName(name: "Owner") { ... on ProjectV2ItemFieldSingleSelectValue { name } }
          content {
            ... on Issue {
              number state body
              repository { nameWithOwner }
              subIssues(first: 50) {
                nodes {
                  id number title state
                  repository { nameWithOwner }
                  assignees(first: 10) { nodes { login } }
                  closedByPullRequestsReferences(first: 10, includeClosedPrs: true) { nodes { number } }
                  linkedBranches(first: 1) { totalCount }
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


def gh(*args, body=None):
    out = subprocess.run(["gh", *args], input=None if body is None else json.dumps(body),
                         capture_output=True, text=True, check=True).stdout
    return json.loads(out) if out.strip() else None


def graphql(query, **variables):
    result = gh("api", "graphql", "--input", "-", body={"query": query, "variables": variables})
    return result["data"]


def project_items():
    after = None
    while True:
        project = graphql(ITEMS, org=ORG, number=PROJECT, after=after)["organization"]["projectV2"]
        yield project, project["items"]["nodes"]
        page = project["items"]["pageInfo"]
        if not page["hasNextPage"]:
            return
        after = page["endCursor"]


def option_id(field, name):
    return next((o["id"] for o in field["options"] if o["name"] == name), None) if name else None


def open_pr(project, sub, parent, owner):
    repo = f"{OWNER}/{REPO}"
    branch = f"issue-{sub['number']}"
    if DRY_RUN:
        print(f"  would open draft PR {branch} for #{sub['number']} {sub['title']!r}")
        return None

    # A PR needs at least one commit the base branch lacks, so the branch
    # starts with an empty commit on top of the default branch.
    info = graphql("""query($o: String!, $r: String!) { repository(owner: $o, name: $r) {
        id defaultBranchRef { name target { oid ... on Commit { tree { oid } } } } } }""",
                   o=OWNER, r=REPO)["repository"]
    base = info["defaultBranchRef"]
    commit = gh("api", f"repos/{repo}/git/commits", "--input", "-", body={
        "message": f"Start #{sub['number']}", "tree": base["target"]["tree"]["oid"],
        "parents": [base["target"]["oid"]]})
    graphql("""mutation($issue: ID!, $repo: ID!, $oid: GitObjectID!, $name: String!) {
        createLinkedBranch(input: {issueId: $issue, repositoryId: $repo, oid: $oid, name: $name}) {
          linkedBranch { id } } }""",
            issue=sub["id"], repo=info["id"], oid=commit["sha"], name=branch)

    pr = gh("api", f"repos/{repo}/pulls", "--input", "-", body={
        "title": sub["title"], "head": branch, "base": base["name"], "draft": True,
        "body": f"Closes #{sub['number']}\n\nPart of #{parent}."})
    # The PR belongs to whoever has the sub-issue, or else the parent's owner.
    assignees = [a["login"] for a in sub["assignees"]["nodes"]] or ([owner] if owner and owner != "Bucket" else [])
    if assignees:
        gh("api", f"repos/{repo}/issues/{pr['number']}/assignees", "--input", "-", body={"assignees": assignees})

    item = graphql("""mutation($project: ID!, $content: ID!) {
        addProjectV2ItemById(input: {projectId: $project, contentId: $content}) { item { id } } }""",
                   project=project["id"], content=pr["node_id"])["addProjectV2ItemById"]["item"]["id"]
    for field, option in ((project["status"], STATUS), (project["owner"], option_id(project["owner"], owner))):
        if option:
            graphql("""mutation($project: ID!, $item: ID!, $field: ID!, $option: String!) {
                updateProjectV2ItemFieldValue(input: {projectId: $project, itemId: $item, fieldId: $field,
                  value: {singleSelectOptionId: $option}}) { projectV2Item { id } } }""",
                    project=project["id"], item=item, field=field["id"], option=option)
    print(f"  opened {pr['html_url']} for #{sub['number']}")
    return pr["number"]


def list_prs(parent, prs):
    """Write or update the parent's comment that lists its sub-issues' PRs."""
    if not prs:
        return
    body = "\n".join([MARKER, "**PRs for the sub-issues**", ""]
                     + [f"- #{sub} → {', '.join(f'#{n}' for n in numbers)}" for sub, numbers in prs])
    if DRY_RUN:
        print(f"  would list on #{parent}: " + "; ".join(f"#{sub} → {numbers}" for sub, numbers in prs))
        return
    repo = f"{OWNER}/{REPO}"
    comments = gh("api", "--paginate", "--slurp", f"repos/{repo}/issues/{parent}/comments")
    old = next((c for page in comments for c in page if c["body"].startswith(MARKER)), None)
    if old is None:
        gh("api", f"repos/{repo}/issues/{parent}/comments", "--input", "-", body={"body": body})
    elif old["body"] != body:
        gh("api", "-X", "PATCH", f"repos/{repo}/issues/comments/{old['id']}", "--input", "-", body={"body": body})


def set_box(issue, ticked):
    """Put the checkbox in the issue's description, ticked or not."""
    body = (issue["body"] or "").replace("\r\n", "\n")
    want = f"- [{'x' if ticked else ' '}] {BOX}"
    lines = body.split("\n")
    new = "\n".join(want if line.strip() in (f"- [ ] {BOX}", f"- [x] {BOX}") else line for line in lines)
    if want not in new:
        new = f"{new.rstrip()}\n\n{want}" if new.strip() else want
    if new == body:
        return
    print(f"  {'ticking' if ticked else 'unticking' if BOX in body else 'adding'} the checkbox on #{issue['number']}")
    if not DRY_RUN:
        gh("api", "-X", "PATCH", f"repos/{OWNER}/{REPO}/issues/{issue['number']}", "--input", "-", body={"body": new})


def main():
    opened = 0
    for project, items in project_items():
        for item in items:
            issue = item["content"] or {}
            if not ISSUE and (item["status"] or {}).get("optionId") == MAKE_SUBISSUES and issue.get("state") == "OPEN" \
                    and issue["repository"]["nameWithOwner"] == f"{OWNER}/{REPO}" and BOX not in (issue["body"] or ""):
                set_box(issue, False)
            wanted = issue.get("number") == ISSUE if ISSUE else (item["status"] or {}).get("optionId") == STATUS
            if not wanted or issue.get("state") != "OPEN" or issue["repository"]["nameWithOwner"] != f"{OWNER}/{REPO}":
                continue
            owner = (item["owner"] or {}).get("name")
            print(f"#{issue['number']} ({owner or 'no owner'})")
            prs = []
            for sub in issue["subIssues"]["nodes"]:
                numbers = [pr["number"] for pr in sub["closedByPullRequestsReferences"]["nodes"]]
                if sub["state"] != "OPEN" or sub["repository"]["nameWithOwner"] != f"{OWNER}/{REPO}" \
                        or numbers or sub["linkedBranches"]["totalCount"]:
                    pass
                elif opened == MAX_PRS:
                    print(f"  #{sub['number']} waits for the next run ({MAX_PRS} PRs per run)")
                else:
                    numbers = [n for n in [open_pr(project, sub, issue["number"], owner)] if n]
                    opened += 1
                if numbers:
                    prs.append((sub["number"], numbers))
            list_prs(issue["number"], prs)
            if ISSUE:
                set_box(issue, False)


if __name__ == "__main__":
    main()
