"""Exercise the router through the same event/output files used by Actions."""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    ("body", "expected"),
    [
        ("### Method\n\nExact Diagonalization (ED)", "vws100"),
        ("### Method\n\nQuantum Monte Carlo (SSE)", "wistaria"),
        ("### Method\n\nQuantum Monte Carlo (Worm)", "LodePollet"),
        ("### Method\n\nDensity Matrix Renormalization Group (DMRG)", "afeiguin"),
        ("### Method\n\nDynamical Mean Field Theory (DMFT)", "egull"),
        ("### Method\n\nOther", "Ooolab"),
        ("No form fields", "Ooolab"),
        ("### Method\n\n### Details\nQuantum Monte Carlo (Worm)", "Ooolab"),
        ("### Method\n\n$(exit 1)", "Ooolab"),
        (None, "Ooolab"),
    ],
)
def test_event_routing(tmp_path, body, expected):
    event = tmp_path / "event.json"
    output = tmp_path / "output.txt"
    event.write_text(json.dumps({"issue": {"body": body}}), encoding="utf-8")
    script = Path(__file__).resolve().parents[2] / ".github/scripts/route_issue.py"
    subprocess.run(
        [sys.executable, str(script)],
        env={**os.environ, "GITHUB_EVENT_PATH": str(event), "GITHUB_OUTPUT": str(output)},
        check=True,
    )
    assert output.read_text(encoding="utf-8") == f"assignee={expected}\n"
