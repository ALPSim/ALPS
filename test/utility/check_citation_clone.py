"""Exercise actual parapack clone saves and relocation, including rank files."""
import subprocess
import sys
import tempfile

with tempfile.TemporaryDirectory(prefix="alps-citation-clone-") as directory:
    command = sys.argv[2:] + [sys.argv[1], "clone-hdf5", directory]
    result = subprocess.run(command, capture_output=True, text=True, timeout=45)
    assert result.returncode == 0, result.stdout + result.stderr
if len(sys.argv) == 2:
    with tempfile.TemporaryDirectory(prefix="alps-citation-aggregate-") as directory:
        result = subprocess.run([sys.argv[1], "aggregate-hdf5", directory],
                                capture_output=True, text=True, timeout=45)
        assert result.returncode == 0, result.stdout + result.stderr
