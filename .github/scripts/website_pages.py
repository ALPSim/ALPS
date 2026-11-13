#!/usr/bin/env python3
"""List the website pages that describe the files a change touches.

The ALPS website (ALPSim/alpsim.github.io) keeps .github/alps-source-map.json,
which maps glob patterns over this repository to the pages that describe those
files. Usage in CI::

    python3 .github/scripts/website_pages.py --base "$BASE" --head HEAD

It reads the map from the site's main branch (or --map), lists the pages whose
sources changed in base...head, leaving out the map's ignored paths (tests),
and warns about source patterns that no longer match a file at head. It writes
to the job summary when GITHUB_STEP_SUMMARY is set and to standard output
otherwise. It never fails: whether a page needs changing, and the map an
update, is for the author to decide.

Standard library only; keep compatible with Python 3.9.
"""
import argparse
import json
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SITE = "https://github.com/ALPSim/alpsim.github.io/blob/main/"
MAP_PATH = ".github/alps-source-map.json"
MAP_URL = "https://raw.githubusercontent.com/ALPSim/alpsim.github.io/main/" + MAP_PATH

sys.path.insert(0, str(ROOT / ".github/scripts/ci"))
from fingerprint import glob_to_regex  # noqa: E402


def git_lines(*args):
    out = subprocess.run(("git",) + args, cwd=ROOT, stdout=subprocess.PIPE, check=True).stdout
    return [line for line in out.decode("utf-8", "surrogateescape").split("\0") if line]


def load(path=None):
    if path:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    with urllib.request.urlopen(MAP_URL, timeout=30) as response:
        return json.loads(response.read().decode("utf-8"))


def stale(page_map, files):
    """[(entry name, pattern)] for the source patterns that match none of files."""
    return [(entry["name"], pattern) for entry in page_map["entries"] for pattern in entry["sources"]
            if not any(glob_to_regex(pattern).match(f) for f in files)]


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


def report(found, gone=()):
    lines = []
    if found:
        pages = {page for _, _, entry_pages in found for page in entry_pages}
        lines += ["### Website pages to check", "",
                  "This PR changes files that %d page(s) of the ALPS website describe "
                  "(from [`%s`](%s%s)). If one is now out of date, update it in "
                  "[alpsim.github.io](https://github.com/ALPSim/alpsim.github.io) (en, zh-cn and ja) "
                  "or open an issue there." % (len(pages), MAP_PATH, SITE, MAP_PATH), "",
                  "| Area | Changed here | Pages |", "|---|---|---|"]
        for name, hits, entry_pages in found:
            shown = ", ".join("`%s`" % h for h in hits[:3])
            if len(hits) > 3:
                shown += " and %d more" % (len(hits) - 3)
            links = ", ".join("[%s](%scontent/en/%s)" % (page[:-3], SITE, page) for page in entry_pages)
            lines.append("| %s | %s | %s |" % (name, shown, links))
    else:
        lines.append("No website pages describe the files this PR changes.")
    if gone:
        lines += ["", "These patterns in [`%s`](%s%s) match no file here any more; if this PR moved or "
                  "removed those files, update the map in alpsim.github.io:" % (MAP_PATH, SITE, MAP_PATH), ""]
        lines += ["- %s: `%s`" % (name, pattern) for name, pattern in gone]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--base", required=True, help="list the pages affected by base...head")
    parser.add_argument("--head", default="HEAD")
    parser.add_argument("--map", help="a local alps-source-map.json instead of the site's")
    args = parser.parse_args(argv)
    # The listing only informs, so nothing here fails the PR.
    try:
        page_map = load(args.map)
        changed = git_lines("diff", "--name-only", "-z", "--no-renames", "%s...%s" % (args.base, args.head))
        files = git_lines("ls-tree", "-r", "-z", "--name-only", args.head)
    except (OSError, ValueError, subprocess.CalledProcessError) as e:
        print("::warning title=Website pages to check::no listing: %s" % e)
        return 0
    found = affected(page_map, changed)
    gone = stale(page_map, files)
    text = report(found, gone)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as f:
            f.write(text)
        if found:
            print("::notice title=Website pages to check::%d area(s) of the website describe files"
                  " this PR changes; see the job summary" % len(found))
        for name, pattern in gone:
            print("::warning title=Website source map::%s: %r matches no file; update %s in"
                  " alpsim.github.io if this PR moved it" % (name, pattern, MAP_PATH))
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
