#!/usr/bin/env python3
"""Run independently linked providers against shared, temporary archive fixtures."""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alps", required=True, type=Path)
    parser.add_argument("--alpscore", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expect", type=Path, help="Compare measurements with an earlier JSON report")
    args = parser.parse_args()
    programs = {"ALPS": args.alps.resolve(), "ALPSCore": args.alpscore.resolve()}
    measurements = {}
    diagnostics = []

    def run(provider, command, *paths):
        result = subprocess.run(
            [str(programs[provider]), command, *map(str, paths)],
            text=True, capture_output=True, check=True, timeout=60,
        )
        if result.stderr:
            diagnostics.append(f"{provider} {command} {paths}\n{result.stderr}")
        records = {}
        for line in result.stdout.splitlines():
            key, value = line.split("\t", 1)
            if key in records:
                raise ValueError(f"Duplicate probe key: {key}")
            records[key] = value
        if not records:
            raise ValueError(f"Empty probe output: {provider} {command}")
        return records

    # Each command is a fresh process, so overlapping provider symbols cannot mix.
    with tempfile.TemporaryDirectory(prefix="alps-reconciliation-") as scratch:
        scratch = Path(scratch)
        for provider in programs:
            measurements[f"{provider}/params-semantics"] = run(provider, "semantics")
            measurements[f"{provider}/archive-semantics"] = run(
                provider, "archive-semantics", scratch / f"{provider}-semantics.h5"
            )
            for kind in ("archive", "params", "extended-params"):
                path = scratch / f"{provider}-{kind}.h5"
                measurements[f"{provider}/write-{kind}"] = run(provider, f"write-{kind}", path)
                for reader in programs:
                    measurements[f"{provider}->{reader}/{kind}"] = run(reader, f"read-{kind}", path)
                if kind == "params":
                    measurements[f"{provider}->ALPSCore/dictionary"] = run("ALPSCore", "read-dictionary", path)

    report = {
        "environment": {"system": platform.system(), "machine": platform.machine()},
        "probe_sha256": {key: hashlib.sha256(path.read_bytes()).hexdigest() for key, path in programs.items()},
        "measurements": measurements,
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    args.output.with_suffix(".log").write_text("\n".join(diagnostics))
    failures = []
    common_params = {"load": "ok", "integer": "42", "real": "1.25", "text": "example", "boolean": "1", "vector": "ok"}
    core_extended = dict(common_params, unsigned="42", float="1.25", wide=(
        "1099511627776" if int(measurements["ALPSCore/params-semantics"]["native_long_bits"]) >= 64 else "42"
    ))
    for provider in programs:
        if measurements[f"{provider}->{provider}/params"] != common_params:
            failures.append(f"Params self-check failed: {provider}")
    if measurements["ALPSCore->ALPSCore/extended-params"] != core_extended:
        failures.append("Extended params self-check failed: ALPSCore")
    for key, values in measurements.items():
        if "/write-" in key:
            if any(value != "ok" for value in values.values()):
                failures.append(f"Provider self-check failed: {key}")
        if key.endswith("/archive"):
            if any(value != "ok" for name, value in values.items() if not name.endswith("_type_marker")):
                failures.append(f"Typed archive interchange failed: {key}")
    if args.expect:
        expected = json.loads(args.expect.read_text())["measurements"]
        for key in sorted(expected.keys() | measurements.keys()):
            if expected.get(key) != measurements.get(key):
                failures.append(f"Changed characterization: {key}")
    if failures:
        raise SystemExit("\n".join(failures) + f"\nFull report: {args.output}")
    print(f"Recorded {len(measurements)} probe groups in {args.output}")


if __name__ == "__main__":
    main()
