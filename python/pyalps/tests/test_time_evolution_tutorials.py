# Copyright (C) 2026 by the ALPS collaboration
# Part of the ALPS Project — see LICENSE.txt for full license text.
# SPDX-License-Identifier: MIT

"""Run the time-evolution tutorials (tutorials/time-evolution) as tests.

Each script checks its physics with assertions and prints a closing line;
it runs from a scratch directory with a non-interactive plotting backend.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("pyalps.classic_dmrg")
pytest.importorskip("matplotlib")
pytest.importorskip("scipy")

TUTORIALS = Path(__file__).resolve().parents[3] / "tutorials" / "time-evolution"
SCRIPTS = ["time_evolution_intro", "interaction_quench", "domain_wall", "quench_hardcore_boson"]


@pytest.mark.parametrize("name", SCRIPTS)
def test_time_evolution_tutorial(name, tmp_path):
    script = TUTORIALS / f"{name}.py"
    if not script.is_file():
        pytest.skip(f"{script} is not present in this tree")
    environment = dict(os.environ, MPLBACKEND="Agg")
    completed = subprocess.run([sys.executable, str(script)], cwd=str(tmp_path), env=environment,
                               capture_output=True, text=True, timeout=900)
    assert completed.returncode == 0, (
        f"{name}.py failed\n--- stdout ---\n{completed.stdout}\n--- stderr ---\n{completed.stderr}")
    assert "All physics assertions passed." in completed.stdout, completed.stdout


def test_every_time_evolution_tutorial_is_run():
    if not TUTORIALS.is_dir():
        pytest.skip("tutorials/time-evolution is not present in this tree")
    assert sorted(path.stem for path in TUTORIALS.glob("*.py")) == sorted(SCRIPTS)
