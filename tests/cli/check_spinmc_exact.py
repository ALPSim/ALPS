#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""spinmc reproduces exact results of small systems.

- Non-symmetric coupling on a Heisenberg dimer: the reference is the
  two-sphere integral of exp(beta s1.J.s2) at T=0.7 (Gauss-Legendre in
  cos(theta), uniform in phi). The local update used to give -0.1989.
- Ising chain of 6 sites in a field, with the default update: the reference
  enumerates the 64 states. The default used to be the cluster update, which
  ignores the field and gave the h=0 magnetization 0.
- 3-state Potts models with antiferromagnetic bonds, with the default update:
  a ring of 8 sites (transfer-matrix value) and a 4x2 ladder with
  ferromagnetic legs and antiferromagnetic rungs (enumeration of the 3^8
  states). The default used to be the cluster update, which never joins an
  antiferromagnetic bond: infinite temperature (energy density 1/3) on the
  ring and -3.54 instead of -4.68 on the ladder.
- J' sets the second bond type, as in the model library's spin model: J'=0 on the
  anisotropic square lattice leaves four independent Ising rings of 4 sites,
  E = -4 N (t + t^3) / (1 + t^4) with t = tanh(J/T). J' used to be ignored.
Inputs spinmc cannot run correctly stop with an error instead of a result:
a negative T, both T and beta, a two-component Potts field (only the second
number was used), and cluster updates for an antiferromagnetic Potts model.
Usage: check_spinmc_exact.py <spinmc> <directory with lattices.xml>
"""
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

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
        {"Energy Density": (-0.241965, 0.008)},
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
        {"Energy": (-5.28997, 0.03), "Magnetization along Field": (0.58619, 0.008)},
    ),
    "antiferromagnetic Potts model": (
        {
            "LATTICE": "chain lattice",
            "L": "8",
            "MODEL": "Potts",
            "q": "3",
            "J": "-1",
            "T": "1",
            "THERMALIZATION": "5000",
            "SWEEPS": "200000",
            "SEED": "41",
        },
        {"Energy Density": (0.155324, 0.003)},
    ),
    "Potts model with couplings of both signs": (
        {
            "LATTICE": "ladder",
            "L": "4",
            "MODEL": "Potts",
            "q": "3",
            "J0": "1",
            "J1": "-1",
            "T": "1",
            "THERMALIZATION": "5000",
            "SWEEPS": "200000",
            "SEED": "41",
        },
        {"Energy": (-4.676918, 0.05)},
    ),
    "J' sets the second bond type": (
        {
            "LATTICE": "anisotropic square lattice",
            "L": "4",
            "MODEL": "Ising",
            "T": "1.5",
            "J": "1",
            "J'": "0",
            "THERMALIZATION": "5000",
            "SWEEPS": "200000",
            "SEED": "2873",
        },
        {"Energy": (-11.19957, 0.08)},
    ),
}

ISING_RING = {"LATTICE": "chain lattice", "L": "8", "MODEL": "Ising", "J": "1", "SWEEPS": "100"}
POTTS_RING = dict(ISING_RING, MODEL="Potts", q="3", T="1")
REJECTED = {
    "T": (dict(ISING_RING, T="-1.5"), "temperature T must be positive"),
    "T and beta": (dict(ISING_RING, T="1.5", beta="0.1"), "either T or beta"),
    "Potts field": (dict(POTTS_RING, h="0.5 0"), "one value for h"),
    "Potts cluster update": (dict(POTTS_RING, J="-1", UPDATE="cluster"), "antiferromagnetic coupling"),
    "J' and J1": (dict(ISING_RING, LATTICE="ladder", L="4", T="1.5", **{"J'": "0.5", "J1": "0.5"}), "either J' or J1"),
}


def run(spinmc, xml_path, parameters, check=True):
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        text = "".join(f'<PARAMETER name="{name}">{value}</PARAMETER>' for name, value in parameters.items())
        (work / "task.in.xml").write_text(f"<SIMULATION><PARAMETERS>{text}</PARAMETERS></SIMULATION>")
        (work / "job.in.xml").write_text(
            '<JOB><OUTPUT file="job.out.xml"/><TASK status="new"><INPUT file="task.in.xml"/>'
            '<OUTPUT file="task.out.xml"/></TASK></JOB>'
        )
        environment = dict(os.environ, ALPS_XML_PATH=xml_path, ALPS_NO_CITATIONS="1")
        result = subprocess.run([spinmc, "--Tmin", "1", "--Tmax", "1", "--write-xml", "job.in.xml"], cwd=work,
                                env=environment, check=check, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                                text=True)
        if not check:
            return result
        return {a.get("name"): float(a.findtext("MEAN")) for a in ET.parse(work / "task.out.xml").iter("SCALAR_AVERAGE")}


def main():
    spinmc, xml_path = sys.argv[1:]
    failures = []
    for case, (parameters, expected) in CASES.items():
        averages = run(spinmc, xml_path, parameters)
        for name, (exact, tolerance) in expected.items():
            if abs(averages[name] - exact) > tolerance:
                failures.append(f"{case}: {name} {averages[name]}, expected {exact} +- {tolerance}")
    for case, (parameters, message) in REJECTED.items():
        result = run(spinmc, xml_path, parameters, check=False)
        if result.returncode == 0 or message not in result.stderr:
            failures.append(f"{case}: expected an error containing '{message}', got exit status "
                            f"{result.returncode}:\n{result.stderr}")
    if failures:
        sys.exit("\n".join(failures))


if __name__ == "__main__":
    main()
