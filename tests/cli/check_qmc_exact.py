#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Quantum Monte Carlo averages of spin-1/2 chains agree with exact diagonalization.

Each case is a spin-1/2 XXZ chain at T = 0.5,

    H = sum_<ij> [Jz Sz_i Sz_j + Jxy/2 (S+_i S-_j + S-_i S+_j)] - h sum_i Sz_i,

small enough for exact diagonalization. The exact values are from
exact_spin_chain.py (full diagonalization from the formula, independent of
ALPS's model library), for example `exact_spin_chain.py 8 0.5 1 1 0 periodic`.
Every average is compared with its exact value in units of its own error bar;
all comparisons for one application share the level mc_hypothesis.ALPHA.

The applications report "Susceptibility" as beta <M^2> / L, which is the
susceptibility only when <M> = 0, so it is compared only at h = 0.

With --short (Debug builds, about ten times slower) no run has more than
200000 sweeps, fewer are too few for a converged binning analysis of worm.
The comparison stays valid with larger error bars; the Debug build then
checks that the assertions hold and the averages are not far off, and the
Release builds catch the small biases.

Usage: check_qmc_exact.py <application> <path to application> <path to parameter2xml> <directory with lattices.xml> [--short]
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import mc_hypothesis

HEISENBERG_8 = {  # exact_spin_chain.py 8 0.5 1 1 0 periodic
    "Energy Density": -0.3432093, "Specific Heat": 0.363321, "Magnetization^2": 0.5722053,
    "Staggered Magnetization^2": 4.885447, "Susceptibility": 0.1430513, "Staggered Susceptibility": 0.957736,
}
XXZ_8_FIELD = {  # exact_spin_chain.py 8 0.5 0.5 1 0.3 periodic
    "Energy Density": -0.2741086, "Specific Heat": 0.3236334, "Magnetization Density": 0.06466709,
    "Magnetization^2": 1.136245, "Staggered Magnetization^2": 3.45142, "Staggered Susceptibility": 0.6126951,
}
ISING_LIKE_4_FIELD = {  # exact_spin_chain.py 4 0.5 2 1 0.3 periodic
    "Energy Density": -0.6118493, "Magnetization Density": 0.009665461, "Magnetization^2": 0.07217696,
    "Staggered Magnetization^2": 3.11454,
}

SHORT_SWEEPS = 200000

COMMON = {"LATTICE": '"chain lattice"', "MODEL": '"spin"', "T": 0.5, "THERMALIZATION": 10000}

# application: [(label, parameters, observables and exact values)]
CASES = {
    "loop": [
        # The SSE loop algorithm has no longitudinal field
        ("path integral, XXZ L=8 Jz=0.5 h=0.3",
         {"ALGORITHM": '"loop"', "L": 8, "Jz": 0.5, "Jxy": 1, "h": 0.3, "SWEEPS": 200000, "SEED": 4711},
         XXZ_8_FIELD),
        ("SSE, Heisenberg L=8",
         {"ALGORITHM": '"loop; sse"', "L": 8, "J": 1, "SWEEPS": 200000, "SEED": 4712}, HEISENBERG_8),
    ],
    "dirloop_sse": [
        ("Heisenberg L=8", {"L": 8, "J": 1, "SWEEPS": 200000, "SEED": 4721},
         {name: HEISENBERG_8[name] for name in
          ("Energy Density", "Magnetization^2", "Staggered Magnetization^2", "Susceptibility")}),
        ("XXZ L=4 Jz=2 h=0.3", {"L": 4, "Jz": 2, "Jxy": 1, "h": 0.3, "SWEEPS": 200000, "SEED": 4722},
         ISING_LIKE_4_FIELD),
    ],
}


def run(application, binary, parameter2xml, work, parameters, short=False):
    if short:
        parameters = {**parameters, "SWEEPS": min(parameters["SWEEPS"], SHORT_SWEEPS)}
    text = "".join(f"{name}={value}\n" for name, value in {**COMMON, **parameters}.items())
    (work / "parm").write_text(text)
    commands = [[parameter2xml, "parm"]]
    if application == "loop":
        commands.append([binary, "--write-xml", "parm.in.xml"])
    else:
        commands.append([binary, "--Tmin", "1", "--Tmax", "1", "--write-xml", "parm.in.xml"])
    for command in commands:
        done = subprocess.run(command, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if done.returncode != 0:
            sys.exit(f"{' '.join(map(str, command))} failed with status {done.returncode}:\n{done.stdout}")
    return mc_hypothesis.read_averages(work / "parm.task1.out.xml")


def main():
    arguments = [argument for argument in sys.argv[1:] if argument]  # CMake passes an empty one outside Debug
    short = "--short" in arguments
    application, binary, parameter2xml, xml_path = [argument for argument in arguments if argument != "--short"]
    os.environ.update(ALPS_XML_PATH=xml_path, ALPS_NO_CITATIONS="1")
    cases = CASES[application]
    family = sum(len(exact) for _, _, exact in cases)
    failures = []
    for label, parameters, exact in cases:
        with tempfile.TemporaryDirectory() as directory:
            averages = run(application, binary, parameter2xml, Path(directory), parameters, short)
        failures += mc_hypothesis.compare(f"{application} {label} (SEED={parameters['SEED']})",
                                          averages, exact, family=family)
    if failures:
        sys.exit("\n".join(failures))


if __name__ == "__main__":
    main()
