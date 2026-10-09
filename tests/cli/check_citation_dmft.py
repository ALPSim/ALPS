"""Check built-in solver selection notices before reading numerical input."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("generate_citations", root / ".github/scripts/generate_citations.py")
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)
_, policy, references, framework = generator.load_catalog(root)
# A missing initial Green function ends each path deterministically after solver
# selection. This keeps the regression independent of Monte Carlo convergence.
for solver, component, omega in (("Hirsch-Fye", "hirschfye", 0),
                                  ("Interaction Expansion", "interaction", 1)):
    for disabled in (False, True):
        env = os.environ.copy()
        env.pop("ALPS_NO_CITATIONS", None)
        if disabled:
            env["ALPS_NO_CITATIONS"] = "1"
        with tempfile.TemporaryDirectory(prefix="alps-dmft-citations-") as cwd:
            Path(cwd, "parameters").write_text(f'''SOLVER="{solver}"
OMEGA_LOOP={omega}
N=8
NMATSUBARA=8
BETA=1
U=1
MU=0
H=0
FLAVORS=2
SITES=1
CONVERGED=0.01
SYMMETRIZATION=1
t=0.5
G0TAU_INPUT="missing-green-function"
G0OMEGA_INPUT="missing-green-function"
''')
            result = subprocess.run([sys.argv[1], "parameters"], cwd=cwd, env=env,
                                    text=True, capture_output=True, timeout=45)
        assert result.returncode != 0, result.stdout
        assert "could not open inital G0 file" in result.stderr, result.stdout + result.stderr
        for key in ("dmft", component):
            expected = generator.notice(policy, references, framework, key)
            assert result.stdout.count(expected) == (0 if disabled else 1), result.stdout
        assert result.stdout.count("Recommended citations for ") == (0 if disabled else 2), result.stdout
