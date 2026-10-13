"""The native test watchdog preserves subprocess results and enforces its limit."""

from pathlib import Path
import subprocess
import sys

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / ".github/scripts/run_with_timeout.py"


@pytest.mark.parametrize("code", [0, 7])
def test_command_exit_status_and_output_are_preserved(code):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "5", sys.executable, "-c",
         f"print('compiler output'); raise SystemExit({code})"],
        text=True, capture_output=True, timeout=10,
    )
    assert result.returncode == code
    assert "compiler output" in result.stdout


def test_timeout_stops_the_process():
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "0.1", sys.executable, "-c",
         "import time; time.sleep(30)"],
        text=True, capture_output=True, timeout=5,
    )
    assert result.returncode == 124
