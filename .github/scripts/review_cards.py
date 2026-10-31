#!/usr/bin/env python3
"""Give each reviewer of a PR in "Generated PRs" a copy of its card.

For each PR on the project board (github.com/orgs/ALPSim/projects/5), or only PR
number PR if that is set:
- its Issue field is set to the issue numbers it closes, or NA if none;
- if its Status is "Generated PRs", every requested or past reviewer without a
  card gets a draft card titled like the PR, in "PRs ready for review", under
  their own Owner name strip (added to the Owner field if missing), showing the
  PR owner's icon and the same Issue value;
- its Reviews field shows one light per reviewer, in name order, and each
  reviewer's card shows their own: 🟢 approved, 🔴 changes requested,
  🟡 waiting or only commented.

Running it again only adds what is missing. With DRY_RUN=1 it only prints what
it would do.
"""
import json
import os
import subprocess

ORG, PROJECT = "ALPSim", 5
OWNER, REPO = "ALPSim", "ALPS"
# Status options, by ID so renaming the columns doesn't break this.
GENERATED = "1d217230"
READY_FOR_REVIEW = "099ec7f0"
PR = int(os.environ.get("PR") or 0)
DRY_RUN = os.environ.get("DRY_RUN") == "1"

ITEMS = """
query($org: String!, $number: Int!, $after: String) {
  organization(login: $org) {
    projectV2(number: $number) {
      id
      status: field(name: "Status") { ... on ProjectV2SingleSelectField { id } }
      owner: field(name: "Owner") { ... on ProjectV2SingleSelectField { id options { id name color description } } }
      issue: field(name: "Issue") { ... on ProjectV2Field { id } }
      reviews: field(name: "Reviews") { ... on ProjectV2Field { id } }
      items(first: 100, after: $after) {
        pageInfo { hasNextPage endCursor }
        nodes {
          id
          status: fieldValueByName(name: "Status") { ... on ProjectV2ItemFieldSingleSelectValue { optionId } }
          issue: fieldValueByName(name: "Issue") { ... on ProjectV2ItemFieldTextValue { text } }
          reviews: fieldValueByName(name: "Reviews") { ... on ProjectV2ItemFieldTextValue { text } }
          content {
            __typename
            ... on DraftIssue { body }
            ... on PullRequest {
              number title url state
              repository { nameWithOwner }
              author { login }
              assignees(first: 10) { nodes { id login } }
              closingIssuesReferences(first: 10) { nodes { number } }
              reviewRequests(first: 20) { nodes { requestedReviewer { ... on User { login } } } }
              latestReviews(first: 20) { nodes { state author { login } } }
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
    return gh("api", "graphql", "--input", "-", body={"query": query, "variables": variables})["data"]


def load_board():
    project, items, after = None, [], None
    while True:
        page = graphql(ITEMS, org=ORG, number=PROJECT, after=after)["organization"]["projectV2"]
        project = project or page
        items += page["items"]["nodes"]
        if not page["items"]["pageInfo"]["hasNextPage"]:
            return project, items
        after = page["items"]["pageInfo"]["endCursor"]


# A re-requested reviewer is waiting again whatever they said before.
LIGHTS = {"APPROVED": "🟢", "CHANGES_REQUESTED": "🔴"}
WAITING = "🟡"


def marker(pr, login):
    return f"<!-- review-card {OWNER}/{REPO}#{pr} @{login} -->"


def set_text(project, item, field, value):
    """Set a text field on a card unless it already holds value."""
    if (item[field] or {}).get("text") != value and not DRY_RUN:
        set_field(project, item["id"], project[field]["id"], {"text": value})


def set_field(project, item, field, value):
    graphql("""mutation($project: ID!, $item: ID!, $field: ID!, $value: ProjectV2FieldValue!) {
        updateProjectV2ItemFieldValue(input: {projectId: $project, itemId: $item, fieldId: $field, value: $value}) {
          projectV2Item { id } } }""", project=project["id"], item=item, field=field, value=value)


def owner_option(project, login):
    """The Owner option for login, added to the field if it isn't there yet."""
    options = project["owner"]["options"]
    found = next((o["id"] for o in options if o["name"].lower() == login.lower()), None)
    if found or DRY_RUN:
        return found
    # The field is rewritten with every option; existing ones keep their IDs so
    # no card loses its Owner.
    keep = [{"id": o["id"], "name": o["name"], "color": o["color"], "description": o["description"]} for o in options]
    field = graphql("""mutation($field: ID!, $options: [ProjectV2SingleSelectFieldOptionInput!]) {
        updateProjectV2Field(input: {fieldId: $field, singleSelectOptions: $options}) {
          projectV2Field { ... on ProjectV2SingleSelectField { options { id name color description } } } } }""",
                    field=project["owner"]["id"],
                    options=keep + [{"name": login, "color": "GRAY", "description": ""}])["updateProjectV2Field"]["projectV2Field"]
    project["owner"]["options"] = field["options"]
    print(f"  added name strip {login}")
    return next(o["id"] for o in field["options"] if o["name"] == login)


def main():
    project, items = load_board()
    cards = {(item["content"] or {}).get("body", "").split("\n", 1)[0]: item for item in items
             if (item["content"] or {}).get("__typename") == "DraftIssue"}
    for item in items:
        pr = item["content"] or {}
        if pr.get("__typename") != "PullRequest" or pr["repository"]["nameWithOwner"] != f"{OWNER}/{REPO}" \
                or (PR and pr["number"] != PR):
            continue
        issue = ", ".join(f"#{i['number']}" for i in pr["closingIssuesReferences"]["nodes"]) or "NA"
        print(f"#{pr['number']} Issue: {issue}")
        set_text(project, item, "issue", issue)
        if pr["state"] != "OPEN" or (item["status"] or {}).get("optionId") != GENERATED:
            continue

        author = (pr["author"] or {}).get("login")
        lights = {r["author"]["login"]: LIGHTS.get(r["state"], WAITING)
                  for r in pr["latestReviews"]["nodes"] if r["author"]}
        lights.update({r["requestedReviewer"]["login"]: WAITING
                       for r in pr["reviewRequests"]["nodes"] if r["requestedReviewer"]})
        lights.pop(author, None)
        reviewers = sorted(lights, key=str.lower)
        print(f"  Reviews: " + " ".join(f"{lights[r]} {r}" for r in reviewers))
        set_text(project, item, "reviews", "".join(lights[r] for r in reviewers))
        # The card shows the PR owner's icon, so its assignee is the PR's.
        owner_ids = [a["id"] for a in pr["assignees"]["nodes"] if a["login"] == author] \
            or [a["id"] for a in pr["assignees"]["nodes"]][:1]
        for login in reviewers:
            if marker(pr["number"], login) in cards:
                set_text(project, cards[marker(pr["number"], login)], "reviews", lights[login])
                continue
            print(f"  card for reviewer {login}")
            if DRY_RUN:
                continue
            option = owner_option(project, login)
            draft = graphql("""mutation($project: ID!, $title: String!, $body: String!, $assignees: [ID!]) {
                addProjectV2DraftIssue(input: {projectId: $project, title: $title, body: $body, assigneeIds: $assignees}) {
                  projectItem { id } } }""", project=project["id"], title=pr["title"], assignees=owner_ids,
                            body=f"{marker(pr['number'], login)}\nReview {pr['url']}")["addProjectV2DraftIssue"]["projectItem"]["id"]
            set_field(project, draft, project["status"]["id"], {"singleSelectOptionId": READY_FOR_REVIEW})
            set_field(project, draft, project["owner"]["id"], {"singleSelectOptionId": option})
            set_field(project, draft, project["issue"]["id"], {"text": issue})
            set_field(project, draft, project["reviews"]["id"], {"text": lights[login]})


if __name__ == "__main__":
    main()
