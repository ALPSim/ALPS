"""Check that p2h5 writes one parameter set and rejects other input with a message."""
from pathlib import Path
import subprocess
import sys
import tempfile


def run(program, *args, text, success):
    result = subprocess.run([program, *map(str, args)], input=text, capture_output=True,
                            text=True, timeout=15)
    assert (result.returncode == 0) == success, (result.returncode, result.stderr)
    return result.stderr


def main():
    program = sys.argv[1]
    with tempfile.TemporaryDirectory(prefix="alps p2h5 ") as directory:
        output = Path(directory) / "parms.h5"
        run(program, output, text="U = 3\nSOLVER = \"hybridization\"\n", success=True)
        assert output.exists()
        output.unlink()
        for text, runs in (("L = 4\n{ T = 1 }\n", "1 run"), ("L = 4\n{ T = 1 }\n{ T = 2 }\n", "2 runs")):
            message = run(program, output, text=text, success=False)
            assert runs in message and "parameter2hdf5" in message, message
            assert not output.exists()
        message = run(program, output, text="U = = 3\n", success=False)
        assert message.startswith("p2h5: "), message
        assert "usage" in run(program, text="", success=False)


if __name__ == "__main__":
    main()
