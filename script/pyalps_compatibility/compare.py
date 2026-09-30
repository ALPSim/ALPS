"""Compare independently collected records, failing on unclassified changes."""

import argparse
import json
from pathlib import Path


def compare(legacy, candidate, expected):
    for key in ("python", "numpy"):
        if legacy["environment"][key] != candidate["environment"][key]:
            raise ValueError(
                f"Compare matching {key} versions, not different environments"
            )
    before, after = legacy["records"], candidate["records"]
    differences = {
        name: {"legacy": before.get(name), "nanobind": after.get(name)}
        for name in sorted(before.keys() | after.keys())
        if before.get(name) != after.get(name)
    }
    classified = {
        name: {key: entry[key] for key in ("legacy", "nanobind")}
        for name, entry in expected.items()
    }
    unexpected = {
        name: entry
        for name, entry in differences.items()
        if classified.get(name) != entry
    }
    stale = sorted(classified.keys() - differences.keys())
    return {
        "records": len(before.keys() | after.keys()),
        "identical": len(before.keys() | after.keys()) - len(differences),
        "classified_differences": len(differences) - len(unexpected),
        "unexpected": unexpected,
        "stale_classifications": stale,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("legacy", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--expected",
        type=Path,
        default=Path(__file__).with_name("expected_differences.json"),
    )
    args = parser.parse_args()
    result = compare(
        *(
            json.loads(path.read_text())
            for path in (args.legacy, args.candidate, args.expected)
        )
    )
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    return bool(result["unexpected"] or result["stale_classifications"])


if __name__ == "__main__":
    raise SystemExit(main())
