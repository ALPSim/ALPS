"""Time CI build phases and report actual compiler-cache hits for this run."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import time


DIRECTORY = Path(os.environ.get("ALPS_CI_METRICS", "_build/metrics"))


def cache_stats():
    try:
        result = subprocess.run(["ccache", "--print-stats"], capture_output=True, text=True)
    except OSError:
        return {}
    if result.returncode:
        return {}
    return {key: int(value) for key, value in (line.split() for line in result.stdout.splitlines())}


def run(label, command, timeout=None):
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    before = cache_stats()
    start = time.monotonic()
    code = 1
    try:
        code = subprocess.run(command, timeout=timeout).returncode
    except subprocess.TimeoutExpired:
        print(f"{label} exceeded its {timeout:g}-second timeout", flush=True)
        code = 124
    except OSError as error:
        print(f"{label}: {error}", flush=True)
        code = 127
    finally:
        elapsed = time.monotonic() - start
        after = cache_stats()
        delta = {key: value - before.get(key, 0) for key, value in after.items()}
        record = {"phase": label, "seconds": elapsed, "exit_code": code,
                  "hits": delta.get("direct_cache_hit", 0) + delta.get("preprocessed_cache_hit", 0),
                  "misses": delta.get("cache_miss", 0)}
        with (DIRECTORY / "phases.jsonl").open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(record) + "\n")
        if after:
            stats = subprocess.run(["ccache", "--show-stats", "--verbose"], capture_output=True, text=True)
            if stats.returncode:
                # The manylinux system package reports full statistics without
                # supporting the newer --verbose option.
                stats = subprocess.run(["ccache", "--show-stats"], capture_output=True, text=True)
            if stats.returncode == 0:
                (DIRECTORY / "ccache.txt").write_text(stats.stdout, encoding="utf-8")
        print(f"{label}: {elapsed:.1f}s, exit {code}, {record['hits']} cache hits, {record['misses']} misses", flush=True)
    return code


def summarize(directory, restored_key=""):
    lines = ["### Build timings and compiler cache", "",
             f"Compiler cache restored: **{'yes' if restored_key else 'no (cold cache)'}**.", "",
             "| Phase | Seconds | Exit code | Cache hits | Cache misses |",
             "| --- | ---: | ---: | ---: | ---: |"]
    phases = directory / "phases.jsonl"
    if phases.exists():
        for line in phases.read_text(encoding="utf-8").splitlines():
            row = json.loads(line)
            phase = row["phase"].replace("|", "&#124;").replace("\n", " ")
            lines.append(f"| {phase} | {row['seconds']:.1f} | {row['exit_code']} | {row['hits']} | {row['misses']} |")
    else:
        lines.extend(["", "No timed phases completed; check the setup steps."])
    stats = directory / "ccache.txt"
    if stats.exists():
        lines.extend(["", "```text", stats.read_text(encoding="utf-8").rstrip(), "```"])
    return "\n".join(lines) + "\n"


def cache_changed(directory):
    """A restored cache needs no new immutable snapshot on an all-hit run."""
    phases = directory / "phases.jsonl"
    return phases.exists() and any(
        json.loads(line).get("misses", 0) > 0 for line in phases.read_text().splitlines()
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="operation", required=True)
    execute = commands.add_parser("run")
    execute.add_argument("--timeout", type=float)
    execute.add_argument("label")
    execute.add_argument("command", nargs=argparse.REMAINDER)
    commands.add_parser("summary")
    commands.add_parser("cache-policy")
    args = parser.parse_args()
    if args.operation == "cache-policy":
        value = str(cache_changed(DIRECTORY)).lower()
        print(f"Compiler cache contains new objects: {value}")
        if output := os.environ.get("GITHUB_OUTPUT"):
            with Path(output).open("a", encoding="utf-8") as stream:
                stream.write(f"changed={value}\n")
        return
    if args.operation == "run":
        if not args.command:
            parser.error("run requires a command")
        raise SystemExit(run(args.label, args.command, args.timeout))
    result = summarize(DIRECTORY, os.environ.get("CCACHE_RESTORED_KEY", ""))
    print(result)
    if summary := os.environ.get("GITHUB_STEP_SUMMARY"):
        with Path(summary).open("a", encoding="utf-8") as stream:
            stream.write(result)


if __name__ == "__main__":
    main()
