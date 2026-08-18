#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Replace ALPS license boilerplate with project and SPDX notices.

The repository is MIT-licensed (see LICENSE.txt), but many files still carry
either boilerplate referring to the superseded "ALPS Library License" or
"ALPS Application License" (see issue #108), or a full copy of the MIT text
inside an ALPS Project header.  This script rewrites those blocks mechanically
so the resulting sweep can be reviewed by re-running the script instead of
reading every changed file:

    python3 script/relicense.py            # dry run: list affected files
    python3 script/relicense.py --write    # rewrite files in place
    python3 script/relicense.py --check    # exit 1 if rewrites remain

The rewrite rule: a recognized contiguous license comment block ending at the
line matching "DEALINGS IN THE SOFTWARE." is replaced by two concise lines:

    <comment-leader> ALPS Project: https://alps.comp-phys.org/
    <comment-leader> SPDX-License-Identifier: MIT

Recognized blocks start with either "This software is part of the ALPS ..." or
the MIT permission grant.  Full MIT blocks are rewritten only when a nearby
comment line identifies the header as part of the ALPS Project; unclassified
and third-party notices are left intact.  The comment leader and the right
edge of boxed headers are preserved, as are copyright attribution lines above
the block.  Legacy license phrases in any other form are never modified; they
are reported for manual attention.

The script is idempotent and only inspects text files; it operates on
raw bytes so files with non-UTF-8 author names are passed through
unchanged apart from the replaced block.
"""

import argparse
import re
import sys
from pathlib import Path

COMMENT_LEADER = rb"[ \t]*(?:#\*?|//+|\*|!)"
LEGACY_START_RE = re.compile(
    rb"^(?P<lead>" + COMMENT_LEADER + rb")[ \t]+This software is part of the ALPS",
    re.IGNORECASE | re.MULTILINE,
)
MIT_START_RE = re.compile(
    rb"^(?P<lead>" + COMMENT_LEADER
    + rb")[ \t]+Permission is hereby granted, free of charge, "
    + rb"to any person[ \t]+obtaining",
    re.IGNORECASE | re.MULTILINE,
)
PROJECT_HEADER_RE = re.compile(
    rb"^" + COMMENT_LEADER
    + rb"[ \t]+ALPS(?:[ \t]+[A-Za-z]+)?[ \t]+Project\b",
    re.IGNORECASE | re.MULTILINE,
)
BOX_PAD_RE = re.compile(rb"[ \t]{2,}(?P<suffix>\*#|[*#])[ \t]*$")
END_RE = re.compile(rb"DEALINGS IN THE SOFTWARE", re.IGNORECASE)
LEGACY_RE = re.compile(
    rb"ALPS[ \t]+(?:Librar(?:y|ies)|Applications?)[ \t]+Licen[cs]e",
    re.IGNORECASE,
)
SPDX = b"SPDX-License-Identifier: MIT"
PROJECT = b"ALPS Project: https://alps.comp-phys.org/"

# The recognized legacy and MIT blocks are at most 21 lines.
MAX_BLOCK_LINES = 25
PROJECT_LOOKBACK_LINES = 40

SKIP_DIRS = {".git", ".hg", "build", "__pycache__"}
SELF = Path(__file__).resolve()


def rewrite(data):
    """Return (new_data, replaced_blocks) for one file's bytes."""
    lines = data.splitlines(keepends=True)
    out = []
    replaced = 0
    i = 0
    while i < len(lines):
        m = LEGACY_START_RE.match(lines[i])
        if not m:
            m = MIT_START_RE.match(lines[i])
            header = b"".join(lines[max(0, i - PROJECT_LOOKBACK_LINES):i])
            if m and not PROJECT_HEADER_RE.search(header):
                m = None
        if m:
            for j in range(i, min(i + MAX_BLOCK_LINES, len(lines))):
                if END_RE.search(lines[j]):
                    eol = b"\r\n" if lines[j].endswith(b"\r\n") else b"\n"
                    # Preserve the right edge of boxed comment blocks.
                    content = lines[i].rstrip(b"\r\n")
                    pad = BOX_PAD_RE.search(content)
                    for notice in (PROJECT, SPDX):
                        replacement = m.group("lead") + b" " + notice
                        suffix = pad.group("suffix") if pad else b""
                        if (pad and len(replacement)
                                < len(content.rstrip()) - len(suffix)):
                            width = len(content.rstrip())
                            replacement += (
                                b" " * (width - len(suffix) - len(replacement))
                                + suffix
                            )
                        out.append(replacement + eol)
                    replaced += 1
                    i = j + 1
                    break
            else:
                out.append(lines[i])
                i += 1
        else:
            out.append(lines[i])
            i += 1
    return b"".join(out), replaced


def candidate_files(paths):
    for root in paths:
        root = Path(root)
        files = [root] if root.is_file() else sorted(
            p for p in root.rglob("*")
            if p.is_file() and not p.is_symlink()
            and not (SKIP_DIRS & set(p.parts))
        )
        for path in files:
            if path.resolve() == SELF:
                continue
            data = path.read_bytes()
            if b"\0" in data[:8192]:  # binary
                continue
            if (LEGACY_RE.search(data) or LEGACY_START_RE.search(data)
                    or MIT_START_RE.search(data)):
                yield path, data


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="*", default=["."],
                        help="files or directories to process (default: .)")
    parser.add_argument("--write", action="store_true",
                        help="rewrite files in place (default: dry run)")
    parser.add_argument("--check", action="store_true",
                        help="report rewriteable or unhandled ALPS license "
                             "blocks and exit 1 if any are found")
    args = parser.parse_args()

    rewritten, pending, manual = [], [], []
    for path, data in candidate_files(args.paths or ["."]):
        new_data, replaced = rewrite(data)
        if args.check:
            if replaced:
                pending.append((path, replaced))
            if LEGACY_RE.search(new_data):
                manual.append(path)
            continue
        if replaced:
            rewritten.append((path, replaced))
            if args.write:
                path.write_bytes(new_data)
        if LEGACY_RE.search(new_data):
            manual.append(path)

    if args.check:
        for path, replaced in pending:
            print(f"rewriteable ALPS license block remains: {path} "
                  f"({replaced} block(s))")
        for path in manual:
            print(f"legacy license phrase remains: {path}")
        failures = set(path for path, _ in pending) | set(manual)
        print(f"{len(failures)} file(s) with rewriteable or unhandled "
              "ALPS license blocks")
        return 1 if failures else 0

    verb = "rewrote" if args.write else "would rewrite"
    for path, replaced in rewritten:
        print(f"{verb} {path} ({replaced} block(s))")
    print(f"{verb} {len(rewritten)} file(s)")
    for path in manual:
        print(f"needs manual attention (legacy phrase outside the standard "
              f"block): {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
