#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exchange Monte Carlo of the Ising chain under MPI gives the exact chain values.

Runs `exchange --mpi` with the replicas spread over the ranks
(parallel_exchange_worker) and with two ranks per replica group
(multiple_parallel_exchange_worker around parallel_ising_worker). Every
printed `Energy` and `Magnetization^2` must match the exact value of the
periodic 16-site chain at one of the replica temperatures, and every
temperature must be matched by some energy.
Usage: check_mpi.py <input directory> <exchange> <launcher> [launcher arguments...]
The launcher arguments must end with the flag that sets the number of ranks.
"""
import math
import re
import subprocess
import sys
from pathlib import Path

RANKS = 4
SITES = 16
INPUTS = ["mpi-exchange.ip", "mpi-exchange-groups.ip"]
# The inputs' ladder: TEMPERATURE_DISTRIBUTION_TYPE = 2 spaces NUM_REPLICAS = 4
# temperatures evenly from T_MIN = 2 to T_MAX = 4.
TEMPERATURES = [2, 8 / 3, 10 / 3, 4]
VALUE = re.compile(r"^(Energy|Magnetization\^2): ([^;\s]+) \+/- ([^;\s]+)")


def exact(temperature):
    """Energy and <M^2> of the periodic Ising chain, H = -sum s_i s_j."""
    t = math.tanh(1 / temperature)
    z = 1 + t**SITES
    energy = -SITES * (t + t ** (SITES - 1)) / z
    m2 = sum((t**abs(i - j) + t ** (SITES - abs(i - j))) / z for i in range(SITES) for j in range(SITES))
    return {"Energy": energy, "Magnetization^2": m2}


def check(output, name):
    references = [exact(t) for t in TEMPERATURES]
    failures = []
    matched = set()
    count = 0
    for line in output:
        m = VALUE.match(line)
        if not m:
            continue
        count += 1
        quantity, mean, error = m.group(1), float(m.group(2)), float(m.group(3))
        tolerance = 5 * error + 1e-3 * abs(mean)
        hits = [i for i, r in enumerate(references) if abs(mean - r[quantity]) <= tolerance]
        if not hits:
            nearest = min(references, key=lambda r: abs(mean - r[quantity]))[quantity]
            failures.append(f"{name}: {quantity} {mean} +/- {error} matches no temperature (nearest exact {nearest:.4f})")
        elif quantity == "Energy":
            matched.update(hits)
    if count == 0:
        failures.append(f"{name}: no results in the output")
    for i, t in enumerate(TEMPERATURES):
        if i not in matched:
            failures.append(f"{name}: no energy for T = {t} (exact {references[i]['Energy']:.4f})")
    return failures


def main():
    directory, program, *launcher = sys.argv[1:]
    failures = []
    for name in INPUTS:
        text = (Path(directory) / name).read_text()
        run = subprocess.run([*launcher, str(RANKS), program, "--mpi"], input=text, capture_output=True, text=True)
        if run.returncode != 0:
            failures.append(f"{name}: exit status {run.returncode}\n{run.stderr}")
            continue
        failures += check(run.stdout.splitlines(), name)
    if failures:
        sys.exit("\n".join(failures))


if __name__ == "__main__":
    main()
