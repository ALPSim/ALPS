"""Regression coverage for release version drift and stale build artifacts."""

import importlib.util
import io
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import zipfile

from packaging.version import Version
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "script" / "check_release_version.py"
SPEC = importlib.util.spec_from_file_location("check_release_version", SCRIPT)
release = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release)


@pytest.fixture
def versions(tmp_path):
    def write(core="3.0.0", python="3.0.0"):
        (tmp_path / "ALPS_VERSION.txt").write_text(core + "\n")
        project_dir = tmp_path / "bindings/python/pyalps"
        project_dir.mkdir(parents=True, exist_ok=True)
        (project_dir / "pyproject.toml").write_text(
            f'[project]\nname = "pyalps"\nversion = "{python}"\n'
        )
        return tmp_path
    return write


@pytest.mark.parametrize("ref,expected", [
    ("refs/heads/master", "3.0.0"),
    ("refs/tags/v3.0.0", "3.0.0"),
    ("refs/tags/v3.0.0-beta.2", "3.0.0b2"),
    ("refs/tags/v3.0.0-rc.1", "3.0.0rc1"),
])
def test_standalone_dynamic_version(versions, monkeypatch, ref, expected):
    root = versions()
    monkeypatch.delenv("ALPS_VERSION_PRERELEASE", raising=False)
    (root / "bindings/python/pyalps/pyproject.toml").write_text(
        '[project]\nname = "pyalps"\ndynamic = ["version"]\n'
    )
    assert release.check_version(root, ref) == Version(expected)


def test_dynamic_version_rejects_stale_tag_and_conflicting_label(versions, monkeypatch):
    root = versions()
    (root / "bindings/python/pyalps/pyproject.toml").write_text(
        '[project]\nname = "pyalps"\ndynamic = ["version"]\n'
    )
    with pytest.raises(ValueError, match="disagrees with ALPS_VERSION.txt"):
        release.check_version(root, "refs/tags/v2.3.4")
    monkeypatch.setenv("ALPS_VERSION_PRERELEASE", "beta.1")
    with pytest.raises(ValueError, match="disagrees with ALPS_VERSION_PRERELEASE"):
        release.check_version(root, "refs/tags/v3.0.0")


def test_prerelease_sdist_keeps_its_version_without_the_build_environment(tmp_path):
    """Exercise the backend: a user rebuild must keep the published version."""
    repository = SCRIPT.parents[1]
    project = repository / "bindings/python/pyalps"
    core = (repository / "ALPS_VERSION.txt").read_text().strip()
    environment = {**os.environ, "GITHUB_REF": f"refs/tags/v{core}-beta.2"}
    environment.pop("ALPS_VERSION_PRERELEASE", None)
    subprocess.run(
        [sys.executable, "-c",
         "from scikit_build_core.build import build_sdist; "
         "import sys; build_sdist(sys.argv[1])", str(tmp_path)],
        cwd=project, env=environment, check=True, capture_output=True, text=True,
    )
    with tarfile.open(tmp_path / f"pyalps-{core}b2.tar.gz") as archive:
        assert not any("/_vendor/" in name for name in archive.getnames())
        archive.extractall(tmp_path, filter="data")
    environment.pop("GITHUB_REF")
    unpacked = tmp_path / f"pyalps-{core}b2"
    assert (unpacked / "_build_support/runtime_manifest.py").is_file()
    assert (unpacked / "examples/ising/CMakeLists.txt").is_file()
    assert (unpacked / "examples/ising/export2py.cpp").is_file()
    for binding in ("maxent_c", "cthyb", "ctint"):
        assert (unpacked / "cpp" / f"{binding}.cpp").is_file()
    completed = subprocess.run(
        [sys.executable, "-c",
         "from scikit_build_core.build import prepare_metadata_for_build_wheel; "
         "print(prepare_metadata_for_build_wheel('metadata'))"],
        cwd=unpacked, env=environment, check=True, capture_output=True, text=True,
    )
    assert f"pyalps-{core}b2.dist-info" in completed.stdout


@pytest.mark.parametrize("ref", [
    "", "refs/heads/master", "refs/heads/v4.0.0", "refs/pull/142/merge", "refs/tags/v3.0.0"
])
def test_matching_versions(versions, ref):
    assert release.check_version(versions(), ref) == Version("3.0.0")


def test_ci_rejects_original_release_using_github_ref(versions):
    result = subprocess.run(
        [sys.executable, release.__file__, "--root", str(versions("2.3.4", "2.3.4b1"))],
        env={**os.environ, "GITHUB_REF": "refs/tags/v3.0.0"},
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "Release tag v3.0.0 disagrees" in result.stderr


@pytest.mark.parametrize("ref", [
    "", "refs/heads/master", "refs/pull/142/merge", "refs/tags/v3.0.0"
])
def test_sdk_and_python_must_agree_even_on_branches(versions, ref):
    with pytest.raises(ValueError, match="disagrees with ALPS_VERSION.txt"):
        release.check_version(versions("2.3.4", "3.0.0"), ref)


@pytest.mark.parametrize("label,suffix", [
    ("alpha.1", "a1"), ("beta.2", "b2"), ("rc.3", "rc3"), ("dev.4", ".dev4")
])
def test_prerelease_tags_match_pep440_versions(versions, label, suffix):
    python = "3.0.0" + suffix
    version = release.check_version(versions(python=python), "refs/tags/v3.0.0-" + label)
    assert version == Version(python)


@pytest.mark.parametrize("python,tag", [
    ("3.0.0b1", "v3.0.0"), ("3.0.0", "v3.0.0-beta.1"), ("3.0.0b1", "v3.0.0-beta.2")
])
def test_prerelease_cannot_be_published_as_final_or_different_prerelease(versions, python, tag):
    with pytest.raises(ValueError, match="disagrees with pyproject.toml"):
        release.check_version(versions(python=python), "refs/tags/" + tag)


@pytest.mark.parametrize("tag", [
    "v3.0", "v3.0.0-beta", "v3.0.0-final", "v3.0.0_1", "v03.0.0", "v3.0.0+local"
])
def test_malformed_tags_fail_closed(versions, tag):
    with pytest.raises(ValueError, match="Invalid release tag"):
        release.check_version(versions(), "refs/tags/" + tag)


@pytest.mark.parametrize("core", ["3.0", "v3.0.0", "3.0.0-beta.1", "3.0.0\n2.3.4", "03.0.0"])
def test_invalid_numeric_core(versions, core):
    with pytest.raises(ValueError, match="must contain MAJOR.MINOR.PATCH"):
        release.check_version(versions(core=core), "")


@pytest.fixture
def artifact(tmp_path):
    def write(kind, version="3.0.0", metadata_version=None, name="pyalps"):
        metadata = (
            f"Metadata-Version: 2.1\nName: {name}\nVersion: {metadata_version or version}\n"
        ).encode()
        if kind == "wheel":
            path = tmp_path / f"{name}-{version}-cp313-cp313-manylinux_2_28_x86_64.whl"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr(f"{name}-{version}.dist-info/METADATA", metadata)
        else:
            path = tmp_path / f"{name}-{version}.tar.gz"
            with tarfile.open(path, "w:gz") as archive:
                member = tarfile.TarInfo(f"{name}-{version}/PKG-INFO")
                member.size = len(metadata)
                archive.addfile(member, io.BytesIO(metadata))
    return write


@pytest.mark.parametrize("version", ["3.0.0", "3.0.0b1"])
def test_matching_wheel_and_sdist(tmp_path, artifact, version):
    artifact("wheel", version)
    artifact("sdist", version)
    release.check_distributions(tmp_path, Version(version))


def test_stale_artifact_rejects_the_entire_batch(tmp_path, artifact):
    artifact("wheel")
    artifact("sdist", "2.3.4b1")
    with pytest.raises(ValueError, match="expected a pyalps 3.0.0"):
        release.check_distributions(tmp_path, Version("3.0.0"))


@pytest.mark.parametrize("kind", ["wheel", "sdist"])
def test_renaming_an_old_artifact_does_not_fix_its_version(tmp_path, artifact, kind):
    artifact(kind, metadata_version="2.3.4b1")
    with pytest.raises(ValueError, match="metadata does not describe pyalps 3.0.0"):
        release.check_distributions(tmp_path, Version("3.0.0"))


def test_other_project_is_rejected(tmp_path, artifact):
    artifact("wheel", name="other")
    with pytest.raises(ValueError, match="expected a pyalps 3.0.0"):
        release.check_distributions(tmp_path, Version("3.0.0"))


def test_empty_dist_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="No distributions"):
        release.check_distributions(tmp_path, Version("3.0.0"))


def test_unexpected_file_is_rejected(tmp_path):
    (tmp_path / "README.txt").write_text("not a distribution")
    with pytest.raises(ValueError, match="Unexpected distribution"):
        release.check_distributions(tmp_path, Version("3.0.0"))
