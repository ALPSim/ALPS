# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Shell entry points for bundled or configured SDK programs, without pyalps.

Keep this outside pyalps: importing that package also loads the scientific
Python stack, which the native command-line applications do not need.
"""

import os
from pathlib import Path
import runpy
import signal
import sys


def _run(program):
    package = Path(__file__).resolve().parent.parent / "pyalps"
    executable = package / "bin" / program
    if not executable.is_file():
        # Only bindings-only builds record an SDK fallback. Load this small
        # generated file directly: importing pyalps would load its extensions.
        config = package / "pyalps_config.py"
        sdk_bin = (
            runpy.run_path(str(config)).get("ALPS_BIN_INSTALL_DIR", "")
            if config.is_file() else ""
        )
        if sdk_bin:
            sdk_bin = os.environ.get("ALPS_BIN_PATH") or sdk_bin
            executable = Path(sdk_bin).resolve() / program
            # An SDK installed into this Python environment may have had its
            # executable replaced by pip's launcher. Never execute ourselves.
            if (
                executable.is_file() and os.path.isfile(sys.argv[0])
                and os.path.samefile(executable, sys.argv[0])
            ):
                print(
                    f"{program}: the configured ALPS SDK executable {executable} "
                    "is this pyalps launcher. Use an SDK installed in a separate prefix.",
                    file=sys.stderr,
                )
                return 127
    if not executable.is_file():
        print(
            f"{program}: no executable was found at {executable}. "
            "Install a wheel built with PYALPS_BUNDLE_APPLICATIONS=ON, or invoke "
            "the executable from your ALPS SDK by its full path.",
            file=sys.stderr,
        )
        return 127

    env = os.environ.copy()
    # Respect explicit resource overrides, just like pyalps.tools.
    env.setdefault("ALPS_XML_PATH", str(package / "xml"))
    env.setdefault("ALPS_BIN_PATH", str(executable.parent))
    try:
        # Use the bundled or configured SDK path, never PATH: it may contain this
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


def checksign():
    return _run("checksign")


def dirloop_sse():
    return _run("dirloop_sse")


def dmft():
    return _run("dmft")


def dmrg():
    return _run("dmrg")


def fulldiag():
    return _run("fulldiag")


def fulldiag_evaluate():
    return _run("fulldiag_evaluate")


def hirschfye():
    return _run("hirschfye")


def hybridization():
    return _run("hybridization")


def interaction():
    return _run("interaction")


def loop():
    return _run("loop")


def qwl():
    return _run("qwl")


def qwl_evaluate():
    return _run("qwl_evaluate")


def simplemc():
    return _run("simplemc")


def sparsediag():
    return _run("sparsediag")


def spinmc():
    return _run("spinmc")


def spinmc_evaluate():
    return _run("spinmc_evaluate")


def worm():
    return _run("worm")


def worm_evaluate():
    return _run("worm_evaluate")


def parameter2xml():
    return _run("parameter2xml")


def printgraph():
    return _run("printgraph")


def convert2xml():
    return _run("convert2xml")


def snap2vtk():
    return _run("snap2vtk")


def maxent():
    return _run("maxent")


def _transform(program, stylesheet):
    # Import the XSLT dependency only for exporters. Native commands remain
    # independent of both lxml and the scientific Python stack.
    import argparse
    from lxml import etree

    parser = argparse.ArgumentParser(prog=program, description="Export ALPS XML to stdout.")
    parser.add_argument("input", help="input XML file, or - for standard input")
    args = parser.parse_args()
    package = Path(__file__).resolve().parent.parent / "pyalps"
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


def convert2text():
    return _transform("convert2text", "QMCXML2text.xsl")


def plot2text():
    return _transform("plot2text", "plot2text.xsl")


def plot2gp():
    return _transform("plot2gp", "plot2gp.xsl")


def plot2xmgr():
    return _transform("plot2xmgr", "plot2xmgr.xsl")
