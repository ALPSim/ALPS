"""Append CTest and pytest JUnit results to the current Actions job summary."""

import argparse
import glob
import html
import os
from pathlib import Path
import xml.etree.ElementTree as ET


def summarize(paths):
    lines = ["| Report | Passed | Failed | Skipped |", "| --- | ---: | ---: | ---: |"]
    found = False
    for path in paths:
        cases = list(ET.parse(path).getroot().iter("testcase"))
        failed = sum(case.find("failure") is not None or case.find("error") is not None for case in cases)
        skipped = sum(case.find("skipped") is not None for case in cases)
        label = html.escape(Path(path).as_posix()).replace("|", "&#124;")
        lines.append(f"| {label} | {len(cases) - failed - skipped} | {failed} | {skipped} |")
        found = True
    return "\n".join(lines) + "\n" if found else "No test reports were produced; check the preceding steps.\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reports", nargs="+")
    args = parser.parse_args()
    paths = sorted({path for pattern in args.reports for path in glob.glob(pattern, recursive=True)})
    result = summarize(paths)
    print(result)
    if summary := os.environ.get("GITHUB_STEP_SUMMARY"):
        with Path(summary).open("a", encoding="utf-8") as stream:
            stream.write(result)


if __name__ == "__main__":
    main()
