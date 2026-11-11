#!/usr/bin/env python3
"""List the website pages a change may affect, and check the page map.

.github/website-pages.json maps source paths here to the pages of the ALPS
website (ALPSim/alpsim.github.io) that describe them. Usage in CI::

    python3 .github/scripts/website_pages.py --check
    python3 .github/scripts/website_pages.py --base "$BASE" --head HEAD

--check fails if the map is malformed or a source pattern matches no tracked
file, so renaming a program or tutorial also updates its entry. The site
checks the other half: that every page named in the map exists.

--base/--head lists the pages whose sources changed in base...head, leaving
out the map's ignored paths (tests), in the job summary when GITHUB_STEP_SUMMARY is set and on standard output otherwise. It
only reports; whether a page needs changing is for the author to decide.

Standard library only; keep compatible with Python 3.9.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MAP = ROOT / ".github/website-pages.json"
SITE = "https://github.com/ALPSim/alpsim.github.io/blob/main/content/en/"

sys.path.insert(0, str(ROOT / ".github/scripts/ci"))
from fingerprint import glob_to_regex  # noqa: E402


def git_lines(*args):
    out = subprocess.run(("git",) + args, cwd=ROOT, stdout=subprocess.PIPE, check=True).stdout
    return [line for line in out.decode("utf-8", "surrogateescape").split("\0") if line]


def load(path=MAP):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def problems(page_map, files):
    """What is wrong with the map, given the tracked files."""
    found = []
    entries = page_map.get("entries")
    unmapped = page_map.get("unmapped", {})
    if not isinstance(entries, list) or not entries:
        return ["the map has no entries"]
    mapped = set()
    for i, entry in enumerate(entries):
        name = entry.get("name") or "entry %d" % (i + 1)
        sources, pages = entry.get("sources"), entry.get("pages")
        if not sources or not pages:
            found.append("%s: needs both sources and pages" % name)
            continue
        for page in pages:
            if page.startswith("/") or not page.endswith(".md"):
                found.append("%s: page %r should be a .md path under content/en" % (name, page))
            mapped.add(page)
        for pattern in sources:
            regex = glob_to_regex(pattern)
            if not any(regex.match(f) for f in files):
                found.append("%s: %r matches no tracked file" % (name, pattern))
    for page in sorted(mapped & set(unmapped)):
        found.append("%s is both mapped and listed as unmapped" % page)
    return found


def affected(page_map, changed):
    """[(entry name, changed sources, pages)] for the entries whose sources changed."""
    ignored = [glob_to_regex(p) for p in page_map.get("ignored", [])]
    changed = [f for f in changed if not any(r.match(f) for r in ignored)]
    found = []
    for entry in page_map["entries"]:
        regexes = [glob_to_regex(p) for p in entry["sources"]]
        hits = sorted(f for f in changed if any(r.match(f) for r in regexes))
        if hits:
            found.append((entry["name"], hits, entry["pages"]))
    return found


def report(found):
    if not found:
        return "No website pages describe the files this PR changes.\n"
    pages = {page for _, _, entry_pages in found for page in entry_pages}
    lines = ["### Website pages to check", "",
             "This PR changes code that %d page(s) of the ALPS website describe "
             "(from `.github/website-pages.json`). If one is now out of date, update it "
             "in [alpsim.github.io](https://github.com/ALPSim/alpsim.github.io) "
             "(en, zh-cn and ja) or open an issue there." % len(pages), "",
             "| Area | Changed here | Pages |", "|---|---|---|"]
    for name, hits, entry_pages in found:
        shown = ", ".join("`%s`" % h for h in hits[:3])
        if len(hits) > 3:
            shown += " and %d more" % (len(hits) - 3)
        links = ", ".join("[%s](%s%s)" % (page[:-3], SITE, page) for page in entry_pages)
        lines.append("| %s | %s | %s |" % (name, shown, links))
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check", action="store_true", help="check the map against the tracked files")
    parser.add_argument("--base", help="list the pages affected by base...head")
    parser.add_argument("--head", default="HEAD")
    args = parser.parse_args(argv)
    if not args.check and not args.base:
        parser.error("give --check or --base")
    page_map = load()
    if args.check:
        found = problems(page_map, git_lines("ls-files", "-z"))
        for problem in found:
            print("::error file=.github/website-pages.json::%s" % problem)
        if found:
            return 1
        print("website-pages.json: %d entries, all sources match tracked files" % len(page_map["entries"]))
    if args.base:
        try:
            changed = git_lines("diff", "--name-only", "-z", "--no-renames", "%s...%s" % (args.base, args.head))
        except subprocess.CalledProcessError:
            # The listing only informs, so a missing base never fails the PR.
            print("::warning title=Website pages to check::could not diff %s...%s" % (args.base, args.head))
            return 0
        found = affected(page_map, changed)
        text = report(found)
        summary = os.environ.get("GITHUB_STEP_SUMMARY")
        if summary:
            with open(summary, "a", encoding="utf-8") as f:
                f.write(text)
            if found:
                print("::notice title=Website pages to check::%d area(s) of the website describe code"
                      " this PR changes; see the job summary" % len(found))
        else:
            sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
