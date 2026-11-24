"""dmft's built-in interaction-expansion solver warns when the expansion
order reached MAX_ORDER (#389): it samples the series cut at that order, which
changes the results (n = 0.581 instead of 0.569 with MAX_ORDER = 4 here). The
default MAX_ORDER is not reached and gives no warning."""
from pathlib import Path
import subprocess
import sys
import tempfile

dmft = sys.argv[1]
BASE = '''SOLVER = "Interaction Expansion";
SITES = 1; FLAVORS = 2; SYMMETRIZATION = 1; ANTIFERROMAGNET = 0;
BETA = 5; U = 2; MU = 0.5; H = 0; t = 1; ALPHA = -0.01;
N = 100; NMATSUBARA = 100; OMEGA_LOOP = 1; MAX_IT = 1; CONVERGED = 0.003;
SWEEPS = 50000; THERMALIZATION = 2000; MAX_TIME = 60; SEED = 7;
'''


def run(extra):
    with tempfile.TemporaryDirectory(prefix="alps-dmft-max-order-") as cwd:
        Path(cwd, "parm").write_text(BASE + extra)
        return subprocess.run([dmft, "parm"], cwd=cwd, text=True, capture_output=True, timeout=110)


capped = run("MAX_ORDER = 4;\n")
assert capped.returncode == 0, capped.stdout + capped.stderr
assert "WARNING: CT-INT: the expansion order was at MAX_ORDER = 4 in " in capped.stderr, capped.stderr
assert "the results are those of the expansion up to order 4, not of the full model." in capped.stderr, capped.stderr

full = run("")
assert full.returncode == 0, full.stdout + full.stderr
assert "MAX_ORDER" not in full.stderr, full.stderr
