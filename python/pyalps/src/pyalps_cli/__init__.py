# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Bundled shell commands and ALPS executable selection, without importing pyalps.

Keep this outside pyalps: importing that package also loads the scientific
Python stack, which the native command-line applications do not need.
"""

import os
from pathlib import Path
import json
import signal
import sys


NATIVE_PROGRAMS = (
    "checksign", "dirloop_sse", "dmft", "dmrg", "fulldiag", "fulldiag_evaluate",
    "hirschfye", "hybridization", "interaction", "loop", "qwl", "qwl_evaluate",
    "simplemc", "sparsediag", "spinmc", "spinmc_evaluate", "worm", "worm_evaluate",
    "parameter2xml", "printgraph", "convert2xml", "snap2vtk", "maxent",
)
EXPORTERS = {
    "convert2text": "QMCXML2text.xsl",
    "plot2text": "plot2text.xsl",
    "plot2gp": "plot2gp.xsl",
    "plot2xmgr": "plot2xmgr.xsl",
}


def runtime_directory():
    """Locate CMake-installed resources without importing the scientific stack."""
    package = Path(__file__).resolve().parent.parent / "pyalps"
    if not (package / "runtime.json").is_file():
        # Editable installs keep live Python sources apart from native resources.
        from importlib.metadata import distribution
        package = Path(distribution("pyalps").locate_file("pyalps"))
    return package


def resolve_executable(program):
    """Resolve an explicit path, a bundled program, or a bindings-only SDK tool.

    PATH is deliberately not searched. Only bindings-only builds record an SDK
    fallback; ALPS_BIN_PATH can override that SDK, but never a bundled payload.
    """
    program = os.fspath(program)
    if os.path.dirname(program):
        executable = Path(program).resolve()
    else:
        package = runtime_directory()
        metadata = json.loads((package / "runtime.json").read_text(encoding="utf-8"))
        sdk_bin = metadata.get("sdk_bin", "")
        directory = (os.environ.get("ALPS_BIN_PATH") or sdk_bin) if sdk_bin else package / "bin"
        executable = Path(directory).resolve() / program
    if not executable.is_file():
        raise RuntimeError(
            f"{program}: no executable was found at {executable}. "
            "Use a full executable path, or select an installed SDK with "
            "ALPS_BIN_PATH for a bindings-only build. Bundled wheels require "
            "their own PYALPS_BUNDLE_APPLICATIONS=ON payload."
        )
    return str(executable)


def _run(program):
    package = runtime_directory()
    try:
        executable = Path(resolve_executable(program))
    except RuntimeError as error:
        print(error, file=sys.stderr)
        return 127

    env = os.environ.copy()
    # XML resources can be overridden independently. Auxiliary executables
    # (for example DMFT solvers) must come from the selected installation.
    env.setdefault("ALPS_XML_PATH", str(package / "xml"))
    env["ALPS_BIN_PATH"] = str(executable.parent)
    try:
        # Use the resolved path, never PATH: it may contain this
        # very launcher or an unrelated ALPS install. Replace the process so
        # arguments, streams, exit status and signals reach the native tool.
        # exec preserves ignored signals, including the ones Python ignores
        # at startup. Restore their native defaults as subprocess does.
        for name in ("SIGPIPE", "SIGXFZ", "SIGXFSZ"):
            native_signal = getattr(signal, name, None)
            if native_signal is not None:
                signal.signal(native_signal, signal.SIG_DFL)
        os.execve(str(executable), [str(executable), *sys.argv[1:]], env)
    except OSError as error:
        print(f"{program}: cannot execute {executable}: {error}", file=sys.stderr)
        return 126


def _transform(program, stylesheet):
    # Import the XSLT dependency only for exporters. Native commands remain
    # independent of both lxml and the scientific Python stack.
    import argparse
    from lxml import etree

    parser = argparse.ArgumentParser(prog=program, description="Export ALPS XML to stdout.")
    parser.add_argument("input", help="input XML file, or - for standard input")
    args = parser.parse_args()
    package = runtime_directory()
    xml_dir = Path(os.environ.get("ALPS_XML_PATH", package / "xml"))
    # ALPS XML files can refer to a remote DTD; conversion needs only their
    # contents. Keep relative xsl:include resolution for helpers.xsl.
    xml_parser = etree.XMLParser(load_dtd=False, resolve_entities=False, no_network=True)
    if hasattr(signal, "SIGPIPE"):
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    try:
        source = sys.stdin.buffer if args.input == "-" else args.input
        document = etree.parse(source, xml_parser)
        transform = etree.XSLT(etree.parse(str(xml_dir / stylesheet), xml_parser))
        sys.stdout.buffer.write(bytes(transform(document)))
        sys.stdout.buffer.flush()
    except (OSError, etree.Error) as error:
        print(f"{program}: {error}", file=sys.stderr)
        return 1
    return 0


def main():
    program = Path(sys.argv[0]).name
    if program in EXPORTERS:
        return _transform(program, EXPORTERS[program])
    if program in NATIVE_PROGRAMS:
        return _run(program)
    print(f"{program}: unknown ALPS command", file=sys.stderr)
    return 127
