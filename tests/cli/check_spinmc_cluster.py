#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""spinmc cluster updates give the exact squared magnetizations of an open chain.

For an open chain of classical unit spins with H = -J sum s_i.s_j the
correlations are <s_i.s_j> = u^|i-j|, with u = I1(K)/I0(K) for the XY model and
u = coth(K) - 1/K for the Heisenberg model, K = J/T. The squared uniform and
staggered magnetizations follow as sums over all pairs. With cluster updates
the staggered one used to equal the uniform one for the XY model and for the
antiferromagnetic Heisenberg model.
Usage: check_spinmc_cluster.py <spinmc> <directory with lattices.xml>
"""
import math
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

L = 8
T = 0.5
TOLERANCE = 0.01
CASES = [("XY", 1), ("XY", -1), ("Heisenberg", 1), ("Heisenberg", -1)]


def bessel_ratio(k):
    """I1(k)/I0(k) from the power series."""
    i0 = sum((k / 2) ** (2 * n) / math.factorial(n) ** 2 for n in range(60))
    i1 = sum((k / 2) ** (2 * n + 1) / (math.factorial(n) * math.factorial(n + 1)) for n in range(60))
    return i1 / i0


def exact(model, coupling):
    k = coupling / T
    u = bessel_ratio(k) if model == "XY" else 1 / math.tanh(k) - 1 / k
    pairs = [(i, j) for i in range(L) for j in range(L)]
    return {
        "Magnetization^2": sum(u ** abs(i - j) for i, j in pairs) / L**2,
        "Staggered Magnetization^2": sum((-u) ** abs(i - j) for i, j in pairs) / L**2,
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
            averages = {a.get("name"): float(a.findtext("MEAN"))
                        for a in ET.parse(work / f"task{number}.out.xml").iter("SCALAR_AVERAGE")}
            for name, expected in exact(model, coupling).items():
                if abs(averages[name] - expected) > TOLERANCE:
                    failures.append(f"{model} J={coupling}: {name} = {averages[name]:.4f}, exact {expected:.4f}")
    if failures:
        sys.exit("\n".join(failures))


if __name__ == "__main__":
    main()
