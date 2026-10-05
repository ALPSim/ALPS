"""Exercise the merge gate's selection, including deletion/rename diffs."""

import importlib.util
import json
from pathlib import Path
import subprocess

import pytest


SCRIPT = Path(__file__).resolve().parents[2] / ".github/scripts/ci_policy.py"
SPEC = importlib.util.spec_from_file_location("ci_policy", SCRIPT)
policy = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(policy)


@pytest.mark.parametrize("paths,source,packaging", [
    (["README.md", "tutorials/README.md"], False, "none"),
    (["src/alps/alea/src/alea/observable.C"], True, "none"),
    (["src/alps/numerics/tests/vector_functions.C"], True, "none"),
    (["python/pyalps/src/pyalps/tools.py"], False, "quick"),
    (["tests/pyalps/test_binding_surface.py"], False, "quick"),
    (["src/alps/alea/src/alea/observable.C", "python/pyalps/src/pyalps/tools.py"], True, "quick"),
    (["README.md", "src/alps/alea/include/alps/alea/observable.h"], True, "full"),
    (["cmake/config/config.h.in"], True, "full"),
    (["src/alps/alea/include/alps/alea/simpleobservable.ipp"], True, "full"),
    (["python/pyalps/pyproject.toml"], True, "full"),
    (["src/apps/new/CMakeLists.txt"], True, "full"),
    (["cmake/ALPSConfig.cmake.in"], True, "full"),
    ([".github/workflows/build.yml"], True, "full"),
    (["cmake/vcpkg/vcpkg.json"], True, "full"),
    ([".github/scripts/prepare_dependencies.py"], True, "full"),
    (["future-source-directory/new.cpp"], True, "full"),
    (None, True, "full"),
])
@pytest.mark.parametrize("event", ["pull_request", "merge_group"])
def test_change_selection(paths, source, packaging, event):
    assert policy.select(event, paths) == {"source": source, "packaging": packaging}


@pytest.mark.parametrize("event", ["push", "schedule", "workflow_dispatch"])
def test_non_pr_events_always_validate_packages(event):
    assert policy.select(event, ["README.md"]) == {"source": True, "packaging": "full"}


def test_rename_out_of_source_cannot_hide_deleted_code(tmp_path, monkeypatch):
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=tmp_path, text=True).strip()
    git("init", "--quiet")
    (tmp_path / "src").mkdir()
    (tmp_path / "src/header.hpp").write_text("header\n")
    git("add", ".")
    commit_args = ("-c", "user.name=CI Test", "-c", "user.email=ci@example.invalid",
                   "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "fixture")
    git(*commit_args)
    base = git("rev-parse", "HEAD")
    git("mv", "src/header.hpp", "README.md")
    git(*commit_args)
    monkeypatch.chdir(tmp_path)
    paths = policy.changed_paths("pull_request", {"pull_request": {"base": {"sha": base}}})
    assert set(paths) == {"src/header.hpp", "README.md"}
    assert policy.select("pull_request", paths)["packaging"] == "full"


def test_missing_base_falls_back_to_full_coverage(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert policy.changed_paths("merge_group", {"merge_group": {"base_sha": "0" * 40}}) is None


def test_quick_smoke_matrix_only_uses_built_artifacts():
    matrices = policy.packaging_matrices("quick")
    assert len(matrices["wheel_matrix"]["plat"]) == 1
    assert {p["artifact"] for p in matrices["smoke_matrix"]["plat"]} == {"cibw-wheels-manylinux"}
    assert matrices["smoke_matrix"]["python"] == ["3.11", "3.14"]


def test_actions_output_round_trips_json(tmp_path, monkeypatch):
    output = tmp_path / "output"
    monkeypatch.setenv("GITHUB_OUTPUT", str(output))
    values = policy.select("pull_request", ["README.md"]) | policy.packaging_matrices("none")
    policy.emit(values)
    decoded = dict(line.split("=", 1) for line in output.read_text().splitlines())
    assert decoded["source"] == "false"
    assert decoded["packaging"] == "none"
    assert json.loads(decoded["wheel_matrix"]) == values["wheel_matrix"]


def test_full_packaging_only_builds_and_smokes_supported_platforms():
    matrices = policy.packaging_matrices('full')
    built = {'cibw-wheels-' + p['family'] for p in matrices['wheel_matrix']['plat']}
    smoke = matrices['smoke_matrix']
    assert {p['artifact'] for p in smoke['plat']} <= built
    assert all(p['os'].startswith(('ubuntu-', 'macos-')) for p in smoke['plat'])
    assert not smoke.get('include')
    assert smoke['python'] == ['3.11', '3.12', '3.13', '3.14']


def test_native_windows_is_manual_and_outside_release_gates():
    workflows = SCRIPT.parents[1] / 'workflows'
    native = (workflows / 'windows.yml').read_text()
    triggers = native.split('on:\n', 1)[1].split('permissions:\n', 1)[0]
    assert '  workflow_dispatch:' in triggers
    assert 'workflow_call:' not in triggers
    assert 'pull_request:' not in triggers
    assert 'push:' not in triggers
    assert 'schedule:' not in triggers
    for name in ('build.yml', 'build_wheels.yml'):
        required = (workflows / name).read_text().lower()
        assert 'windows' not in required, f'{name} must not depend on native Windows'
