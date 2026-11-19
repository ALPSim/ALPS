#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""spinmc reproduces exact results of small systems.

- Non-symmetric coupling on a Heisenberg dimer: the reference is the
  two-sphere integral of exp(beta s1.J.s2) at T=0.7 (Gauss-Legendre in
  cos(theta), uniform in phi). The local update used to give -0.1989.
- Ising chain of 6 sites in a field, with the default update: the reference
  enumerates the 64 states. The default used to be the cluster update, which
  ignores the field and gave the h=0 magnetization 0.

Each average is compared with its exact value in units of its error bar
(mc_hypothesis), all comparisons sharing one level.
Usage: check_spinmc_exact.py <spinmc> <directory with lattices.xml>
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import mc_hypothesis

CASES = {
    "nonsymmetric coupling": (
        {
            "LATTICE": "dimer",
            "MODEL": "Heisenberg",
            "T": "0.7",
            "UPDATE": "local",
            "THERMALIZATION": "5000",
            "SWEEPS": "100000",
            "SEED": "2873",
            "J": "1 0.5 0 -0.5 1 0 0 0 1",
        },
        {"Energy Density": -0.241965},
    ),
    "default update in a field": (
        {
            "LATTICE": "chain lattice",
            "L": "6",
            "MODEL": "Ising",
            "T": "1.5",
            "J": "1",
            "h": "0.3",
            "THERMALIZATION": "5000",
            "SWEEPS": "1000000",
            "SEED": "2873",
        },
        {"Energy": -5.28997, "Magnetization along Field": 0.58619},
    ),
}


def run(spinmc, xml_path, parameters):
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        text = "".join(f'<PARAMETER name="{name}">{value}</PARAMETER>' for name, value in parameters.items())
        (work / "task.in.xml").write_text(f"<SIMULATION><PARAMETERS>{text}</PARAMETERS></SIMULATION>")
        (work / "job.in.xml").write_text(
            '<JOB><OUTPUT file="job.out.xml"/><TASK status="new"><INPUT file="task.in.xml"/>'
            '<OUTPUT file="task.out.xml"/></TASK></JOB>'
        )
        environment = dict(os.environ, ALPS_XML_PATH=xml_path, ALPS_NO_CITATIONS="1")
        subprocess.run([spinmc, "--Tmin", "1", "--Tmax", "1", "--write-xml", "job.in.xml"], cwd=work, env=environment,
                       check=True, stdout=subprocess.DEVNULL)
        return mc_hypothesis.read_averages(work / "task.out.xml")


def main():
    spinmc, xml_path = sys.argv[1:]
    family = sum(len(expected) for _, expected in CASES.values())
    failures = []
    for case, (parameters, expected) in CASES.items():
        averages = run(spinmc, xml_path, parameters)
        failures += mc_hypothesis.compare(f"{case} (SEED={parameters['SEED']})", averages, expected, family=family)
    if failures:
        sys.exit("\n".join(failures))


if __name__ == "__main__":
    main()
