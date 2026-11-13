#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Every run of a parallel qwl task finishes and contributes its results.

A qwl run stores its coefficients only when it finishes. Three runs of one
task under MPI must all be merged, and the energy of the 8-site spin-1/2
Heisenberg ring at T=0.1 must agree with exact diagonalization (ground state
-0.456387 per site plus the triplet's thermal weight).
Usage: check_scheduler_qwl_runs.py <qwl> <qwl_evaluate> <directory with lattices.xml> <mpiexec> [flags...]
"""
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

RUNS = 3
EXACT_ENERGY_DENSITY = -0.4554
TOLERANCE = 0.01
PARAMETERS = {
    "LATTICE": "chain lattice",
    "MODEL": "spin",
    "local_S": "1/2",
    "J": "1",
    "L": "8",
    "CUTOFF": "200",
    # the final stage has to outlast several of the scheduler's checks
    "SWEEPS": "600000",
    "SEED": "11",
}


def main():
    qwl, qwl_evaluate, xml_path, *launcher = sys.argv[1:]
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        parameters = "".join(f'<PARAMETER name="{name}">{value}</PARAMETER>' for name, value in PARAMETERS.items())
        (work / "task.in.xml").write_text(f"<SIMULATION><PARAMETERS>{parameters}</PARAMETERS></SIMULATION>")
        (work / "job.in.xml").write_text(
            '<JOB><OUTPUT file="job.out.xml"/><TASK status="new"><INPUT file="task.in.xml"/>'
            '<OUTPUT file="task.out.xml"/></TASK></JOB>'
        )
        environment = dict(os.environ, ALPS_XML_PATH=xml_path, ALPS_NO_CITATIONS="1")
        quiet = dict(cwd=work, env=environment, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([*launcher, qwl, "--mpi", "--Tmin", "1", "--Tmax", "1", "--write-xml", "job.in.xml"], **quiet)
        coefficients = [a for a in ET.parse(work / "task.out.xml").iter("VECTOR_AVERAGE")
                        if a.get("name") == "Coefficients"]
        if not coefficients:
            sys.exit("no run stored its coefficients")
        merged = int(coefficients[0].find("SCALAR_AVERAGE").findtext("COUNT"))
        subprocess.run([qwl_evaluate, "task.out.xml"], **quiet)
        point = next(ET.parse(work / "task.plot.energy.xml").iter("point"))
        temperature, energy = float(point.findtext("x")), float(point.findtext("y"))
    if merged != RUNS:
        sys.exit(f"{merged} of {RUNS} runs stored their coefficients")
    if abs(temperature - 0.1) > 1e-9 or abs(energy - EXACT_ENERGY_DENSITY) > TOLERANCE:
        sys.exit(f"energy density {energy} at T={temperature}, expected {EXACT_ENERGY_DENSITY} +- {TOLERANCE} at T=0.1")


if __name__ == "__main__":
    main()
