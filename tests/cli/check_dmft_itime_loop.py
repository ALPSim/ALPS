"""The imaginary-time DMFT loop (OMEGA_LOOP = 0) refuses what it cannot run
(#388, case 3), before it starts a solver.

- It builds the next hybridization function from the Bethe-lattice condition
  Delta(tau) = t^2 G(tau), so it ran the Bethe lattice whatever TWODBS or
  DOSFILE said, and printed "Semicircular DOS".
- A solver that does not read the hybridization function needs an
  imaginary-time Hilbert transform, which does not exist: the run failed after
  the first solver run.

The Matsubara loop (OMEGA_LOOP = 1) runs both, and the allowed imaginary-time
cases still get past the checks. A missing initial G0 file ends those runs
before any solver runs."""
from pathlib import Path
import subprocess
import sys
import tempfile

dmft = sys.argv[1]
BASE = 'N=8\nNMATSUBARA=8\nBETA=1\nU=1\nMU=0\nH=0\nFLAVORS=2\nSITES=1\nCONVERGED=0.01\n' \
       'SYMMETRIZATION=true\nt=0.5\nG0TAU_INPUT="missing-green-function"\n' \
       'G0OMEGA_INPUT="missing-green-function"\n'


def run(lines):
    with tempfile.TemporaryDirectory(prefix="alps-dmft-itime-loop-") as cwd:
        Path(cwd, "dos.dat").write_text("-1 0.5\n0 0.5\n1 0.5\n")
        Path(cwd, "parameters").write_text(BASE + lines)
        return subprocess.run([dmft, "parameters"], cwd=cwd, text=True,
                              capture_output=True, timeout=45)


def refused(lines, message):
    result = run(lines)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    error = result.stderr.strip()
    assert error.startswith(message) and "\n" not in error, output
    # refused before a solver or the lattice is set up
    assert "using external solver" not in result.stdout, output
    assert "Semicircular DOS" not in result.stdout, output


def accepted(lines):
    result = run(lines)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "could not open inital G0 file" in result.stderr, output
    assert "imaginary-time loop" not in result.stderr, output


lattice = "The imaginary-time loop (OMEGA_LOOP = 0) supports only the Bethe lattice"
for solver in ("hybridization", "ct-hyb-segment"):
    for key, value in (("TWODBS", "1"), ("TWODBS", '"hexagonal"'), ("DOSFILE", '"dos.dat"')):
        refused(f'SOLVER="{solver}"\nOMEGA_LOOP=0\n{key}={value}\n',
                f"{lattice} (semicircular DOS), and would ignore {key}; set OMEGA_LOOP = 1 for {key}.")
        # the Matsubara loop has the general Hilbert transform
        accepted(f'SOLVER="{solver}"\nOMEGA_LOOP=1\n{key}={value}\n')
    # the Bethe lattice
    accepted(f'SOLVER="{solver}"\nOMEGA_LOOP=0\n')

for solver in ("some-solver", "/usr/local/bin/some-solver"):
    refused(f'SOLVER="{solver}"\nOMEGA_LOOP=0\n',
            f"SOLVER = {solver} needs the Matsubara loop; set OMEGA_LOOP = 1. The imaginary-time loop")
    # SC_WRITE_DELTA = 1 hands any solver the hybridization function
    accepted(f'SOLVER="{solver}"\nOMEGA_LOOP=0\nSC_WRITE_DELTA=1\n')

result = run('SOLVER="Hybridization"\nOMEGA_LOOP=0\n')
assert result.returncode != 0 and "replaced by a standalone hybridzation solver" in result.stderr, \
    result.stdout + result.stderr
