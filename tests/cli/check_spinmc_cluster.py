#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""spinmc cluster updates give the exact squared magnetizations of an open chain.

The exact squared uniform and staggered magnetizations are from
classical_chain.py. With cluster updates the staggered one used to equal the
uniform one for the XY model and for the antiferromagnetic Heisenberg model.
The results are compared in units of their error bars (mc_hypothesis), all
eight comparisons sharing one level.
Usage: check_spinmc_cluster.py <spinmc> <directory with lattices.xml>
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import classical_chain
import mc_hypothesis

L = 8
T = 0.5
CASES = [("XY", 1), ("XY", -1), ("Heisenberg", 1), ("Heisenberg", -1)]


def exact(model, coupling):
    k = coupling / T
    return {
        "Magnetization^2": classical_chain.magnetization_squared(model, k, L),
        "Staggered Magnetization^2": classical_chain.magnetization_squared(model, k, L, staggered=True),
    }


def main():
    spinmc, xml_path = sys.argv[1:]
    failures = []
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        tasks = ""
        for number, (model, coupling) in enumerate(CASES, 1):
            parameters = {
                "LATTICE": "open chain lattice", "L": L, "T": T, "MODEL": model, "J": coupling,
                "UPDATE": "cluster", "THERMALIZATION": 5000, "SWEEPS": 200000, "SEED": 2873 + number,
            }
            text = "".join(f'<PARAMETER name="{name}">{value}</PARAMETER>' for name, value in parameters.items())
            (work / f"task{number}.in.xml").write_text(f"<SIMULATION><PARAMETERS>{text}</PARAMETERS></SIMULATION>")
            tasks += (f'<TASK status="new"><INPUT file="task{number}.in.xml"/>'
                      f'<OUTPUT file="task{number}.out.xml"/></TASK>')
        (work / "job.in.xml").write_text(f'<JOB><OUTPUT file="job.out.xml"/>{tasks}</JOB>')
        environment = dict(os.environ, ALPS_XML_PATH=xml_path, ALPS_NO_CITATIONS="1")
        subprocess.run([spinmc, "--Tmin", "1", "--Tmax", "1", "--write-xml", "job.in.xml"], cwd=work, env=environment,
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for number, (model, coupling) in enumerate(CASES, 1):
            averages = mc_hypothesis.read_averages(work / f"task{number}.out.xml")
            failures += mc_hypothesis.compare(f"{model} J={coupling} (SEED={2873 + number})", averages,
                                              exact(model, coupling), family=2 * len(CASES))
    if failures:
        sys.exit("\n".join(failures))


if __name__ == "__main__":
    main()
