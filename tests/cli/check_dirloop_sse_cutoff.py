#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""dirloop_sse grows the operator string when the expansion order nears the cutoff.

A short thermalisation leaves the cutoff too small for the measurement sweeps.
The run used to stop with "nnonzero reached the maximum value"; it now has to
finish and still reproduce the exact energy. The reference is exact
diagonalisation of the S=1/2 Heisenberg ring with L=8 at T=0.5.
Usage: check_dirloop_sse_cutoff.py <dirloop_sse> <directory with lattices.xml>
"""
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

EXACT_ENERGY_DENSITY = -0.343209
TOLERANCE = 0.006
PARAMETERS = {
    "LATTICE": "chain lattice",
    "MODEL": "spin",
    "local_S": "1/2",
    "J": "1",
    "L": "8",
    "T": "0.5",
    "THERMALIZATION": "20",
    "SWEEPS": "50000",
    "SEED": "1",
}


def main():
    dirloop_sse, xml_path = sys.argv[1:]
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        parameters = "".join(
            f'<PARAMETER name="{name}">{value}</PARAMETER>' for name, value in PARAMETERS.items()
        )
        (work / "task.in.xml").write_text(f"<SIMULATION><PARAMETERS>{parameters}</PARAMETERS></SIMULATION>")
        (work / "job.in.xml").write_text(
            '<JOB><OUTPUT file="job.out.xml"/><TASK status="new"><INPUT file="task.in.xml"/>'
            '<OUTPUT file="task.out.xml"/></TASK></JOB>'
        )
        environment = dict(os.environ, ALPS_XML_PATH=xml_path, ALPS_NO_CITATIONS="1")
        subprocess.run([dirloop_sse, "--Tmin", "1", "--Tmax", "1", "--write-xml", "job.in.xml"], cwd=work,
                       env=environment, check=True, stdout=subprocess.DEVNULL)
        averages = ET.parse(work / "task.out.xml").iter("SCALAR_AVERAGE")
        energy = next(float(a.findtext("MEAN")) for a in averages if a.get("name") == "Energy Density")
    if abs(energy - EXACT_ENERGY_DENSITY) > TOLERANCE:
        sys.exit(f"Energy Density {energy}, expected {EXACT_ENERGY_DENSITY} +- {TOLERANCE}")


if __name__ == "__main__":
    main()
