#!/usr/bin/env python3
"""Put every open PR on the board under its author, and give each reviewer a
copy of its card.

Every open PR in the repository, or only PR number PR if that is set, is added
to the project board (github.com/orgs/ALPSim/projects/5) if it isn't there yet,
gets Status "Generated PRs", and gets its author as Owner unless it already has
an owner other than Bucket. Then, for each PR on the board:
- its Issue field is set to the issue numbers it closes, or NA if none;
- while it is open, every requested or past reviewer gets a draft card titled
  "ALPS #<PR> <title>" in the last column ("PR review assignments"), under
  their own Owner name strip (added to the Owner field if missing), showing the
  PR author's icon, the same Issue value and the PR's state;
- its Reviews field shows one light per reviewer, in name order, and each
  reviewer's card shows their own: 🟢 approved, 🔴 changes requested,
  🟡 waiting or only commented;
- once it is merged or closed, or a reviewer is removed before reviewing, those
  reviewer cards are removed. The PR's own card goes to Done through the
  board's built-in workflows.

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

CARD = """
fragment card on ProjectV2Item {
    id
    status: fieldValueByName(name: "Status") { ... on ProjectV2ItemFieldSingleSelectValue { optionId } }
    owner: fieldValueByName(name: "Owner") { ... on ProjectV2ItemFieldSingleSelectValue { name } }
    issue: fieldValueByName(name: "Issue") { ... on ProjectV2ItemFieldTextValue { text } }
    reviews: fieldValueByName(name: "Reviews") { ... on ProjectV2ItemFieldTextValue { text } }
    prstate: fieldValueByName(name: "PR state") { ... on ProjectV2ItemFieldTextValue { text } }
    content {
      __typename
      ... on DraftIssue { id title body }
      ... on PullRequest {
        id number title url state isDraft
        repository { nameWithOwner }
        author { login ... on User { id } }
        closingIssuesReferences(first: 10) { nodes { number } }
        reviewRequests(first: 20) { nodes { requestedReviewer { ... on User { login } } } }
        latestReviews(first: 20) { nodes { state author { login } } }
      }
    }
}
"""

ITEMS = """
query($org: String!, $number: Int!, $after: String) {
  organization(login: $org) {
    projectV2(number: $number) {
      id
      status: field(name: "Status") { ... on ProjectV2SingleSelectField { id } }
      owner: field(name: "Owner") { ... on ProjectV2SingleSelectField { id options { id name color description } } }
      issue: field(name: "Issue") { ... on ProjectV2Field { id } }
      reviews: field(name: "Reviews") { ... on ProjectV2Field { id } }
      prstate: field(name: "PR state") { ... on ProjectV2Field { id } }
      items(first: 100, after: $after) {
        pageInfo { hasNextPage endCursor }
        nodes { ...card }
      }
    }
  }
}
""" + CARD


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


OPEN_PRS = """
query($owner: String!, $repo: String!, $after: String) {
  repository(owner: $owner, name: $repo) {
    pullRequests(states: OPEN, first: 100, after: $after) {
      pageInfo { hasNextPage endCursor }
      nodes { id number }
    }
  }
}
"""


def add_open_prs(project, items):
    """Add the open PRs that aren't on the board yet, and return their cards.

    The board's item list lags behind for a minute or so, so new cards are
    fetched by ID instead."""
    on_board = {(item["content"] or {}).get("id") for item in items}
    added, after = [], None
    while True:
        page = graphql(OPEN_PRS, owner=OWNER, repo=REPO, after=after)["repository"]["pullRequests"]
        for pr in page["nodes"]:
            if pr["id"] in on_board or (PR and pr["number"] != PR):
                continue
            print(f"#{pr['number']} added to the board")
            if not DRY_RUN:
                item = graphql("""mutation($project: ID!, $pr: ID!) {
                    addProjectV2ItemById(input: {projectId: $project, contentId: $pr}) { item { id } } }""",
                               project=project["id"], pr=pr["id"])["addProjectV2ItemById"]["item"]["id"]
                added.append(graphql("query($id: ID!) { node(id: $id) { ...card } }" + CARD, id=item)["node"])
        if not page["pageInfo"]["hasNextPage"]:
            return added
        after = page["pageInfo"]["endCursor"]


def add_card(project, title, body, author_id):
    """Add a draft card showing the PR author's icon, or none if they can't be assigned."""
    query = """mutation($project: ID!, $title: String!, $body: String!, $assignees: [ID!]) {
        addProjectV2DraftIssue(input: {projectId: $project, title: $title, body: $body, assigneeIds: $assignees}) {
          projectItem { id } } }"""
    try:
        data = graphql(query, project=project["id"], title=title, body=body, assignees=[author_id] if author_id else [])
    except subprocess.CalledProcessError:
        # Someone from outside the org, for instance.
        data = graphql(query, project=project["id"], title=title, body=body, assignees=[])
    return data["addProjectV2DraftIssue"]["projectItem"]["id"]


def remove(project, card, why):
    print(f"  removing {why}")
    if not DRY_RUN:
        graphql("""mutation($project: ID!, $item: ID!) {
            deleteProjectV2Item(input: {projectId: $project, itemId: $item}) { deletedItemId } }""",
                project=project["id"], item=card["id"])


def pr_state(pr):
    return {"MERGED": "Merged", "CLOSED": "Closed"}.get(pr["state"]) or ("Draft PR" if pr["isDraft"] else "Open")


def main():
    project, items = load_board()
    items += add_open_prs(project, items)
    # A card made by a run a minute earlier can be missing from the item list
    # and get made twice; keep the first of each and remove the rest.
    cards = {}
    for item in items:
        if (item["content"] or {}).get("__typename") != "DraftIssue":
            continue
        key = item["content"]["body"].split("\n", 1)[0]
        if key.startswith("<!-- review-card ") and key in cards:
            remove(project, item, f"duplicate card {key}")
            continue
        cards[key] = item
    seen = set()
    for item in items:
        pr = item["content"] or {}
        if pr.get("__typename") != "PullRequest" or pr["repository"]["nameWithOwner"] != f"{OWNER}/{REPO}" \
                or (PR and pr["number"] != PR):
            continue
        issue = ", ".join(f"#{i['number']}" for i in pr["closingIssuesReferences"]["nodes"]) or "NA"
        print(f"#{pr['number']} Issue: {issue}")
        set_text(project, item, "issue", issue)
        seen.add(pr["number"])
        mine = {key: card for key, card in cards.items()
                if key.startswith(f"<!-- review-card {OWNER}/{REPO}#{pr['number']} @")}
        # A merged or closed PR leaves no cards in the reviewers' column.
        if pr["state"] != "OPEN":
            for key, card in mine.items():
                remove(project, card, f"{key} (PR {pr['state'].lower()})")
            continue

        # An open PR always sits in "Generated PRs", under its author's name.
        author = (pr["author"] or {}).get("login")
        if (item["status"] or {}).get("optionId") != GENERATED:
            print(f"  Status: Generated PRs")
            if not DRY_RUN:
                set_field(project, item["id"], project["status"]["id"], {"singleSelectOptionId": GENERATED})
        if author and (item["owner"] or {}).get("name") in (None, "Bucket"):
            print(f"  Owner: {author}")
            if not DRY_RUN:
                set_field(project, item["id"], project["owner"]["id"], {"singleSelectOptionId": owner_option(project, author)})

        lights = {r["author"]["login"]: LIGHTS.get(r["state"], WAITING)
                  for r in pr["latestReviews"]["nodes"] if r["author"]}
        # Teams have no login and no name strip, so they get no card.
        lights.update({r["requestedReviewer"]["login"]: WAITING
                       for r in pr["reviewRequests"]["nodes"] if (r["requestedReviewer"] or {}).get("login")})
        lights.pop(author, None)
        reviewers = sorted(lights, key=str.lower)
        print(f"  Reviews: " + " ".join(f"{lights[r]} {r}" for r in reviewers))
        set_text(project, item, "reviews", "".join(lights[r] for r in reviewers))
        # A reviewer removed before reviewing loses their card.
        for login in {key.rsplit("@", 1)[1][:-4] for key in mine} - set(reviewers):
            remove(project, mine[marker(pr["number"], login)], f"{login}'s card (no longer a reviewer)")
        # The card shows the PR author's icon, so its assignee is the author.
        author_id = (pr["author"] or {}).get("id")
        # A draft card can't show the PR's number in its header, so the title starts with it.
        title = f"{REPO} #{pr['number']} {pr['title']}"
        # The reviewers' cards can't show the PR's own state icon, so they spell it out.
        state = pr_state(pr)
        for login in reviewers:
            if marker(pr["number"], login) in cards:
                card = cards[marker(pr["number"], login)]
                set_text(project, card, "reviews", lights[login])
                set_text(project, card, "prstate", state)
                if card["content"]["title"] != title:
                    print(f"  title of {login}'s card: {title}")
                    if not DRY_RUN:
                        graphql("""mutation($draft: ID!, $title: String!) {
                            updateProjectV2DraftIssue(input: {draftIssueId: $draft, title: $title}) { draftIssue { id } } }""",
                                draft=card["content"]["id"], title=title)
                continue
            print(f"  card for reviewer {login}")
            if DRY_RUN:
                continue
            option = owner_option(project, login)
            draft = add_card(project, title, f"{marker(pr['number'], login)}\nReview {pr['url']}", author_id)
            set_field(project, draft, project["status"]["id"], {"singleSelectOptionId": READY_FOR_REVIEW})
            set_field(project, draft, project["owner"]["id"], {"singleSelectOptionId": option})
            set_field(project, draft, project["issue"]["id"], {"text": issue})
            set_field(project, draft, project["reviews"]["id"], {"text": lights[login]})
            set_field(project, draft, project["prstate"]["id"], {"text": state})

    # Reviewer cards whose PR was taken off the board go once it is merged or closed.
    prefix = f"<!-- review-card {OWNER}/{REPO}#"
    for key, card in cards.items():
        if not key.startswith(prefix):
            continue
        number = int(key[len(prefix):].split(" ", 1)[0])
        if number in seen or (PR and number != PR):
            continue
        pr = graphql("""query($o: String!, $r: String!, $n: Int!) {
            repository(owner: $o, name: $r) { pullRequest(number: $n) { state } } }""",
                     o=OWNER, r=REPO, n=number)["repository"]["pullRequest"]
        if pr["state"] != "OPEN":
            remove(project, card, f"{key} (PR {pr['state'].lower()})")


if __name__ == "__main__":
    main()
