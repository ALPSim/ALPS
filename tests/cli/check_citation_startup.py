"""A tiny real looper run: two stdin tasks must share one startup notice."""
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
expected = generator.notice(policy, references, framework, "looper")
parameters = f'''LATTICE_LIBRARY="{Path(sys.argv[2]) / 'lattices.xml'}"
MODEL_LIBRARY="{Path(sys.argv[2]) / 'models.xml'}"
LATTICE="chain lattice"
MODEL="spin"
local_S=1/2
L=4
J=1
THERMALIZATION=2
SWEEPS=8
SEED=17
ALGORITHM="loop"
{{T=1;}}
{{T=2;}}
'''
for disabled in (False, True):
    env = os.environ.copy()
    env.pop("ALPS_NO_CITATIONS", None)
    if disabled:
        env["ALPS_NO_CITATIONS"] = "1"
    with tempfile.TemporaryDirectory(prefix="alps-citation-calculation-") as cwd:
        result = subprocess.run([sys.argv[1]], cwd=cwd, env=env, input=parameters, text=True,
                                capture_output=True, timeout=45)
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stderr.count(expected) == (0 if disabled else 1), result.stderr
    assert result.stderr.count("Recommended citations for ") == (0 if disabled else 1), result.stderr
    assert "copyright (c)" in result.stderr, result.stderr
    assert "Licensed under the MIT License." in result.stderr, result.stderr
    assert "Recommended citations for " not in result.stdout, result.stdout
    assert result.stdout.startswith("[input parameters]"), result.stdout
    assert result.stdout.count("[results]") == 2, result.stdout
    if disabled:
        assert result.stdout == numerical_stdout, result.stdout
    else:
        numerical_stdout = result.stdout
