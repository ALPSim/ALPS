# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
import re
import subprocess
import sys
from pathlib import Path

SCRIPT = (Path(__file__).resolve().parents[2]
          / 'tutorials' / 'optical-lattice-01-bandstructure' / 'bandstructure.py')


def values(output, label):
    line = next(l for l in output.splitlines() if l.startswith(label))
    return [float(v) for v in re.findall(r'[-+]?\d+\.\d*(?:e[-+]?\d+)?', line.split('=', 1)[1])]


def test_bandstructure_tutorial(tmp_path):
    # Reference values printed by the removed pyalps.dwa.bandstructure for the same lattice.
    out = subprocess.run([sys.executable, str(SCRIPT)], cwd=str(tmp_path),
                         capture_output=True, text=True, check=True).stdout
    assert all(abs(v - 154.89) < 0.01 for v in values(out, 'Er2nK'))
    assert all(abs(v - 4.77051) < 1e-4 for v in values(out, 't [nK]'))
    assert abs(values(out, 'U [nK]')[0] - 38.7018) < 1e-3
    assert all(abs(v - 8.11272) < 1e-4 for v in values(out, 'U/t'))
    assert all(abs(v - 1) < 1e-9 for v in values(out, 'norm'))


if __name__ == '__main__':
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        test_bandstructure_tutorial(Path(d))
