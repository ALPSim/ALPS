"""The native test watchdog preserves subprocess results and enforces its limit."""

from pathlib import Path
import subprocess
import sys

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / ".github/scripts/run_with_timeout.py"
pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="CI watchdog uses POSIX process groups")


@pytest.mark.parametrize("code", [0, 7])
def test_command_exit_status_and_output_are_preserved(code):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "5", sys.executable, "-c",
         f"print('compiler output'); raise SystemExit({code})"],
        text=True, capture_output=True, timeout=10,
    )
    assert result.returncode == code
    assert "compiler output" in result.stdout


def test_timeout_stops_nested_processes_and_closes_output():
    # The real wheel runner launches pytest as a child. Its inherited pipes
    # keep capture_output blocked if the watchdog kills only the launcher.
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "1", sys.executable, "-c",
         "import signal, subprocess, sys; signal.signal(signal.SIGTERM, signal.SIG_IGN); "
         "subprocess.run([sys.executable, '-c', "
         "\"import time; print('child ready', flush=True); time.sleep(12)\"])"],
        text=True, capture_output=True, timeout=9,
    )
    assert result.returncode == 124
    assert "child ready" in result.stdout
