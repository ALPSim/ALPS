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
    "parameter2xml", "printgraph", "convert2xml", "snap2vtk", "maxent", "alps2fcidump",
)
# The same commands and formats as the SDK's alps-xml (src/tools/xml), but
# rendered with lxml so that wheels need no external xsltproc.
XML_STYLES = {
    "plot": {"text": "plot2text.xsl", "html": "ALPS.xsl", "gnuplot": "plot2gp.xsl",
             "matplotlib": "plot2mpl.xsl", "grace": "plot2xmgr.xsl"},
    "convert": {"text": "QMCXML2text.xsl", "html": "ALPS.xsl"},
}
XML_STYLES["extract"] = XML_STYLES["plot"]
# Historical names, kept as aliases of alps-xml commands.
EXPORTERS = {
    "convert2text": ("convert", "text"),
    "plot2text": ("plot", "text"),
    "plot2gp": ("plot", "gnuplot"),
    "plot2xmgr": ("plot", "grace"),
}
COMMANDS = (*NATIVE_PROGRAMS, "alps-xml", *EXPORTERS)


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


def _render(command, style, inputs, specification=None):
    # Import the XSLT dependency only for XML commands. Native commands remain
    # independent of both lxml and the scientific Python stack.
    from lxml import etree

    xml_dir = Path(os.environ.get("ALPS_XML_PATH", runtime_directory() / "xml")).resolve()
    # ALPS XML files can refer to a remote DTD; conversion needs only their
    # contents. Keep relative xsl:include resolution for helpers.xsl.
    parser = etree.XMLParser(load_dtd=False, resolve_entities=False, no_network=True)

    def parse(source):
        return etree.parse(sys.stdin.buffer if str(source) == "-" else str(source), parser)

    def stylesheet(name):
        return etree.XSLT(etree.parse(str(xml_dir / name), parser))

    if command != "extract":
        return bytes(stylesheet(style)(parse(inputs[0])))
    archive = etree.Element("ARCHIVE")
    for source in inputs:
        root = parse(source).getroot()
        if root.tag == "SIMULATION":
            archive.append(root)
        elif root.tag == "ARCHIVE":
            archive.extend(root.findall("SIMULATION"))
        else:
            raise ValueError(f"{source}: expected SIMULATION or ARCHIVE XML")
    if not len(archive):
        raise ValueError("The input contains no simulations")
    # Reparse each intermediate result, as the SDK command does with files:
    # the stylesheets emit part of their markup as unescaped text.
    generated = etree.fromstring(bytes(stylesheet("plot2xsl.xsl")(parse(specification))), parser)
    # The generated stylesheet has no location; resolve its library includes
    # against the installed resources.
    for include in generated.findall("{http://www.w3.org/1999/XSL/Transform}include"):
        include.set("href", (xml_dir / include.get("href")).as_uri())
    plot = etree.fromstring(bytes(etree.XSLT(generated)(etree.ElementTree(archive))), parser)
    return bytes(stylesheet(style)(etree.ElementTree(plot)))


def _emit(program, command, format, inputs, specification=None, output=None):
    from lxml import etree

    if hasattr(signal, "SIGPIPE"):
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    try:
        result = _render(command, XML_STYLES[command][format], inputs, specification)
        # Publish output only after the transformation succeeded.
        if output is None:
            sys.stdout.buffer.write(result)
            sys.stdout.buffer.flush()
        else:
            output.write_bytes(result)
    except (OSError, ValueError, etree.Error) as error:
        print(f"{program}: {error}", file=sys.stderr)
        return 1
    return 0


def _xml_tool(program):
    import argparse

    parser = argparse.ArgumentParser(
        prog=program, description="Render ALPS plots and simulation XML.")
    commands = parser.add_subparsers(dest="command", required=True)
    for command, help_text in (
        ("plot", "render a plot XML file"),
        ("convert", "render simulation or job XML"),
        ("extract", "extract a plot from simulation or archive XML"),
    ):
        subparser = commands.add_parser(command, help=help_text)
        subparser.add_argument("format", choices=XML_STYLES[command])
        if command == "extract":
            subparser.add_argument("specification", help="plot definition XML")
        subparser.add_argument("input", nargs="+" if command == "extract" else None)
        subparser.add_argument("-o", "--output", type=Path,
                               help="write to a file instead of stdout")
    args = parser.parse_args()
    inputs = args.input if args.command == "extract" else [args.input]
    return _emit(program, args.command, args.format, inputs,
                 getattr(args, "specification", None), args.output)


def _exporter(program):
    import argparse

    command, format = EXPORTERS[program]
    parser = argparse.ArgumentParser(
        prog=program, description=f"Same as 'alps-xml {command} {format}'.")
    parser.add_argument("input", help="input XML file, or - for standard input")
    return _emit(program, command, format, [parser.parse_args().input])


def main():
    program = Path(sys.argv[0]).name
    if program == "alps-xml":
        return _xml_tool(program)
    if program in EXPORTERS:
        return _exporter(program)
    if program in NATIVE_PROGRAMS:
        return _run(program)
    print(f"{program}: unknown ALPS command", file=sys.stderr)
    return 127
