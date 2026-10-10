"""Run installed-pyalps checks and retain reproducible evidence; no SDK build."""

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--wheelhouse", type=Path)
    parser.add_argument("--packaging", action="store_true")
    parser.add_argument("--downstream", action="store_true")
    parser.add_argument("--applications", action="store_true")
    parser.add_argument("--legacy-python", type=Path)
    parser.add_argument("--legacy-modules", type=Path)
    args = parser.parse_args()
    if bool(args.legacy_python) != bool(args.legacy_modules):
        parser.error("--legacy-python and --legacy-modules must be supplied together")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    import pyalps

    package = Path(pyalps.__file__).resolve().parent
    environment = os.environ.copy()
    if args.downstream:
        environment["PYALPS_TEST_DOWNSTREAM_EXPORT"] = "1"
    manifest = {
        "source_revision": git("rev-parse", "HEAD"),
        "source_tree": git("rev-parse", "HEAD^{tree}"),
        "source_status": git("status", "--short"),
        "dirty_file_hashes": {
            name: sha256(ROOT / name) if (ROOT / name).is_file() else None
            for name in git(
                "ls-files", "--modified", "--others", "--exclude-standard"
            ).splitlines()
        },
        "python": sys.version,
        "platform": platform.platform(),
        "packages": dict(
            sorted(
                (d.metadata["Name"], d.version)
                for d in importlib.metadata.distributions()
            )
        ),
        "installed_package": str(package),
        "installed_extension_hashes": {
            str(path.relative_to(package)): sha256(path)
            for path in sorted(package.rglob("*"))
            if path.suffix in {".so", ".pyd", ".dll", ".dylib"}
        },
        "downstream_enabled": environment.get("PYALPS_TEST_DOWNSTREAM_EXPORT") == "1",
        "build_environment": {
            name: os.environ[name]
            for name in ("CPLUS_INCLUDE_PATH", "ALPS_DIR", "CMAKE_ARGS")
            if name in os.environ
        },
        "steps": [],
    }
    if args.wheelhouse:
        manifest["distribution_hashes"] = {
            path.name: sha256(path)
            for path in sorted(args.wheelhouse.iterdir())
            if path.name.endswith((".whl", ".tar.gz"))
        }
    (output / "source.patch").write_text(git("diff", "HEAD") + "\n")

    def run(name, command, env=None):
        command = [str(arg) for arg in command]
        print(f"Running {name}", flush=True)
        start = time.monotonic()
        with (output / f"{name}.log").open("w") as log:
            try:
                result = subprocess.run(
                    command,
                    cwd=ROOT,
                    env=env or environment,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=300,
                )
                code = result.returncode
            except subprocess.TimeoutExpired:
                code = 124
        manifest["steps"].append(
            {
                "name": name,
                "command": command,
                "returncode": code,
                "seconds": round(time.monotonic() - start, 2),
            }
        )
        if code:
            raise RuntimeError(
                f"{name} failed ({code}); see {output / (name + '.log')}"
            )

    success = False
    try:
        tests = ["tests/pyalps"] + (["tests/packaging"] if args.packaging else [])
        run(
            "pytest",
            [
                sys.executable,
                "-m",
                "pytest",
                "-q",
                *tests,
                "--junitxml",
                output / "pytest.xml",
            ],
        )
        scripts = Path(__file__).resolve().parent / "pyalps_compatibility"
        if args.legacy_python:
            old_env = {**environment, "PYTHONPATH": str(args.legacy_modules.resolve())}
            old_python = args.legacy_python.absolute()
            run(
                "legacy-probes",
                [old_python, scripts / "probe.py", "legacy", output / "legacy.json"],
                old_env,
            )
            run(
                "candidate-probes",
                [
                    sys.executable,
                    scripts / "probe.py",
                    "nanobind",
                    output / "candidate.json",
                ],
            )
            run(
                "comparison",
                [
                    sys.executable,
                    scripts / "compare.py",
                    output / "legacy.json",
                    output / "candidate.json",
                    output / "comparison.json",
                ],
            )
            for writer, reader in (("legacy", "nanobind"), ("nanobind", "legacy")):
                for mode, action in ((writer, "write"), (reader, "read")):
                    run(
                        f"checkpoint-{writer}-{action}",
                        [
                            old_python if mode == "legacy" else sys.executable,
                            scripts / "checkpoint.py",
                            mode,
                            action,
                            output / f"{writer}.h5",
                        ],
                        old_env if mode == "legacy" else environment,
                    )
        if args.applications:
            for app in (
                "spinmc",
                "loop",
                "dirloop_sse",
                "sparsediag",
                "fulldiag",
                "dmrg",
            ):
                run(
                    app,
                    [sys.executable, scripts / "applications.py", app],
                    {
                        **environment,
                        "PYALPS_WORKFLOW_ROOT": str(output / "applications"),
                    },
                )
        success = True
    finally:
        manifest["success"] = success
        if (output / "pytest.xml").exists():
            manifest["test_suites"] = [
                suite.attrib
                for suite in ET.parse(output / "pytest.xml").iter("testsuite")
            ]
        (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Validation passed; evidence: {output}")


if __name__ == "__main__":
    main()
