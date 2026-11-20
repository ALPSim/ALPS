#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""simplemc averages of open classical spin chains agree with the exact ones.

Each case is an open chain of L = 8 classical spins with H = -J sum s_i.s_j at
T = 1, ferromagnetic and antiferromagnetic, for the Ising, XY and Heisenberg
algorithms. The exact values are from classical_chain.py. Every average is
compared with its exact value in units of its own error bar; all comparisons
share the level mc_hypothesis.ALPHA.

Usage: check_simplemc_exact.py <simplemc> <parameter2xml> <directory with lattices.xml>
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import classical_chain
import mc_hypothesis

L = 8
T = 1.0
SWEEPS = 1000000
# algorithm, model, the averages of one spin component and the number of components
MODELS = [("ising", "Ising", "Magnetization Density^2", 1), ("xy", "XY", "Magnetization Density X^2", 2),
          ("heisenberg", "Heisenberg", "Magnetization Density Z^2", 3)]
CASES = [(algorithm, model, component, n, coupling) for algorithm, model, component, n in MODELS for coupling in (1, -1)]


def exact(model, component, n, coupling):
    k = coupling / T
    m2 = classical_chain.magnetization_squared(model, k, L)
    values = {
        "Energy Density": classical_chain.energy_density(model, k, coupling, L),
        "Specific Heat": classical_chain.specific_heat(model, k, L),
        "Magnetization Density^2": m2,
    }
    if n > 1:
        values[component] = m2 / n  # each component carries 1/n of <M^2>
    return values


def main():
    simplemc, parameter2xml, xml_path = sys.argv[1:]
    environment = dict(os.environ, ALPS_XML_PATH=xml_path, ALPS_NO_CITATIONS="1")
    cases = [(f"{model} J={coupling} (SEED={4800 + number})", exact(model, component, n, coupling))
             for number, (_, model, component, n, coupling) in enumerate(CASES, 1)]
    family = sum(len(values) for _, values in cases)
    failures = []
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        tasks = "".join(f'{{ ALGORITHM="{algorithm}"; J={coupling}; SEED={4800 + number} }}\n'
                        for number, (algorithm, _, _, _, coupling) in enumerate(CASES, 1))
        (work / "parm").write_text(f'LATTICE="open chain lattice"\nL={L}\nT={T}\n'
                                   f'THERMALIZATION=10000\nSWEEPS={SWEEPS}\n{tasks}')
        for command in ([parameter2xml, "parm"], [simplemc, "--write-xml", "parm.in.xml"]):
            done = subprocess.run(command, cwd=work, env=environment, stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, text=True)
            if done.returncode != 0:
                sys.exit(f"{' '.join(command)} failed with status {done.returncode}:\n{done.stdout}")
        for number, (label, values) in enumerate(cases, 1):
            averages = mc_hypothesis.read_averages(work / f"parm.task{number}.out.xml")
            failures += mc_hypothesis.compare(label, averages, values, family=family)
    if failures:
        sys.exit("\n".join(failures))


if __name__ == "__main__":
    main()
