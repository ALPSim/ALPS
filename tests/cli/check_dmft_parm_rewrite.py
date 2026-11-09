"""Check that a finished dmft run succeeds when it cannot rewrite its input file."""
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile

PARM = '''SOLVER = "Hirsch-Fye"
OMEGA_LOOP = 1
ANTIFERROMAGNET = 1
H = 0
H_INIT = 0.0225
FLAVORS = 2
SITES = 1
SYMMETRIZATION = 0
N = 16
NMATSUBARA = 64
N_MEAS = 100
BETA = 6
U = 3
MU = 0
t = 0.707106781186547
SEED = 0
SWEEPS = 200
THERMALIZATION = 20
MAX_TIME = 20
MAX_IT = 1
CONVERGED = 0.003
TOLERANCE = 0.001
'''


def run(program, read_only):
    with tempfile.TemporaryDirectory(prefix="alps-dmft-rewrite-") as directory:
        parm = Path(directory, "parm")
        parm.write_text(PARM)
        if read_only:
            parm.chmod(stat.S_IRUSR)
            if os.access(parm, os.W_OK):  # root writes to read-only files
                return
        result = subprocess.run([program, "parm"], cwd=directory, capture_output=True,
                                text=True, timeout=120)
        assert result.returncode == 0, result.stdout + result.stderr
        assert Path(directory, "parm.h5").exists()
        warning = "Warning: could not rewrite parameter file 'parm'" in result.stderr
        assert warning == read_only, result.stderr
        assert (parm.read_text() == PARM) == read_only, parm.read_text()


def main():
    run(sys.argv[1], read_only=False)
    run(sys.argv[1], read_only=True)


if __name__ == "__main__":
    main()
