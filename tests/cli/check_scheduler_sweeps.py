#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""A scheduler-based Monte Carlo run stops at SWEEPS and is reproducible.

spinmc runs a job of three tasks twice with the same SEEDs and the default
check times. Each task must take exactly SWEEPS measurements, give the same
energy to the last digit, and end when the work is done instead of at the next
check. With an MPI launcher, the second run uses MPI with one run per task
(--Nmax 1) and must give the same results as the serial run.
Usage: check_scheduler_sweeps.py <spinmc> <directory with lattices.xml>
                                 [<mpiexec> <numproc flag> <ranks> [preflags...]]
"""
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

SWEEPS = 20000
SEEDS = (2873, 511, 77)
PARAMETERS = {
    "LATTICE": "square lattice",
    "L": "4",
    "MODEL": "Ising",
    "J": "1",
    "T": "2.269186",
    "UPDATE": "cluster",
    "THERMALIZATION": "2000",
    "SWEEPS": str(SWEEPS),
}


def run(command, xml_path, work):
    work.mkdir()
    tasks = ""
    for i, seed in enumerate(SEEDS):
        parameters = "".join(f'<PARAMETER name="{name}">{value}</PARAMETER>'
                             for name, value in dict(PARAMETERS, SEED=seed).items())
        (work / f"task{i}.in.xml").write_text(f"<SIMULATION><PARAMETERS>{parameters}</PARAMETERS></SIMULATION>")
        tasks += f'<TASK status="new"><INPUT file="task{i}.in.xml"/><OUTPUT file="task{i}.out.xml"/></TASK>'
    (work / "job.in.xml").write_text(f'<JOB><OUTPUT file="job.out.xml"/>{tasks}</JOB>')
    environment = dict(os.environ, ALPS_XML_PATH=xml_path, ALPS_NO_CITATIONS="1")
    # The default check interval is 60 s: the timeout fails a run that waits for it.
    subprocess.run(command + ["--write-xml", "job.in.xml"], cwd=work, env=environment, check=True,
                   stdout=subprocess.DEVNULL, timeout=40)
    results = []
    for i in range(len(SEEDS)):
        energy = next(a for a in ET.parse(work / f"task{i}.out.xml").iter("SCALAR_AVERAGE")
                      if a.get("name") == "Energy")
        results.append((int(energy.findtext("COUNT")), energy.findtext("MEAN"), energy.findtext("ERROR")))
    return results


def main():
    spinmc, xml_path, *mpi = sys.argv[1:]
    second_command = [spinmc]
    if mpi:
        mpiexec, numproc_flag, ranks, *preflags = mpi
        second_command = [mpiexec, numproc_flag, ranks, *preflags, spinmc, "--mpi", "--Nmax", "1"]
    with tempfile.TemporaryDirectory() as directory:
        first = run([spinmc], xml_path, Path(directory) / "first")
        second = run(second_command, xml_path, Path(directory) / "second")
    for result in first:
        if result[0] != SWEEPS:
            sys.exit(f"{result[0]} measurements, expected SWEEPS = {SWEEPS}")
    if first != second:
        sys.exit(f"two runs with the same SEEDs differ: {first} and {second}")


if __name__ == "__main__":
    main()
