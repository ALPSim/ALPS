#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""spinmc samples a non-symmetric coupling matrix on a two-site Heisenberg dimer.

The reference is the two-sphere integral of exp(beta s1.J.s2) with J below at
T=0.7 (Gauss-Legendre in cos(theta), uniform in phi). The local update used to
give -0.1989.
Usage: check_spinmc_coupling.py <spinmc> <directory with lattices.xml>
"""
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

EXACT_ENERGY_DENSITY = -0.241965
TOLERANCE = 0.008
PARAMETERS = {
    "LATTICE": "dimer",
    "MODEL": "Heisenberg",
    "T": "0.7",
    "UPDATE": "local",
    "THERMALIZATION": "5000",
    "SWEEPS": "100000",
    "SEED": "2873",
    "J": "1 0.5 0 -0.5 1 0 0 0 1",
}


def main():
    spinmc, xml_path = sys.argv[1:]
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
        subprocess.run([spinmc, "--Tmin", "1", "--Tmax", "1", "--write-xml", "job.in.xml"], cwd=work, env=environment, check=True,
                       stdout=subprocess.DEVNULL)
        averages = ET.parse(work / "task.out.xml").iter("SCALAR_AVERAGE")
        energy = next(float(a.findtext("MEAN")) for a in averages if a.get("name") == "Energy Density")
    if abs(energy - EXACT_ENERGY_DENSITY) > TOLERANCE:
        sys.exit(f"Energy Density {energy}, expected {EXACT_ENERGY_DENSITY} +- {TOLERANCE}")


if __name__ == "__main__":
    main()
