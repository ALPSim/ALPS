#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""A scheduler-based Monte Carlo run stops at SWEEPS and is reproducible.

spinmc runs the same task twice with the same SEED and the default check
times. Each run must take exactly SWEEPS measurements, give the same energy to
the last digit, and end when the work is done instead of at the next check.
Usage: check_scheduler_sweeps.py <spinmc> <directory with lattices.xml>
"""
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

SWEEPS = 20000
PARAMETERS = {
    "LATTICE": "square lattice",
    "L": "4",
    "MODEL": "Ising",
    "J": "1",
    "T": "2.269186",
    "UPDATE": "cluster",
    "THERMALIZATION": "2000",
    "SWEEPS": str(SWEEPS),
    "SEED": "2873",
}


def run(spinmc, xml_path, work):
    work.mkdir()
    parameters = "".join(f'<PARAMETER name="{name}">{value}</PARAMETER>' for name, value in PARAMETERS.items())
    (work / "task.in.xml").write_text(f"<SIMULATION><PARAMETERS>{parameters}</PARAMETERS></SIMULATION>")
    (work / "job.in.xml").write_text(
        '<JOB><OUTPUT file="job.out.xml"/><TASK status="new"><INPUT file="task.in.xml"/>'
        '<OUTPUT file="task.out.xml"/></TASK></JOB>'
    )
    environment = dict(os.environ, ALPS_XML_PATH=xml_path, ALPS_NO_CITATIONS="1")
    # The default check interval is 60 s: the timeout fails a run that waits for it.
    subprocess.run([spinmc, "--write-xml", "job.in.xml"], cwd=work, env=environment, check=True,
                   stdout=subprocess.DEVNULL, timeout=40)
    energy = next(a for a in ET.parse(work / "task.out.xml").iter("SCALAR_AVERAGE") if a.get("name") == "Energy")
    return int(energy.findtext("COUNT")), energy.findtext("MEAN"), energy.findtext("ERROR")


def main():
    spinmc, xml_path = sys.argv[1:]
    with tempfile.TemporaryDirectory() as directory:
        first = run(spinmc, xml_path, Path(directory) / "first")
        second = run(spinmc, xml_path, Path(directory) / "second")
    if first[0] != SWEEPS:
        sys.exit(f"{first[0]} measurements, expected SWEEPS = {SWEEPS}")
    if first != second:
        sys.exit(f"two runs with the same SEED differ: {first} and {second}")


if __name__ == "__main__":
    main()
