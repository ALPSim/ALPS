"""Regression coverage for release version drift and stale build artifacts."""

import importlib.util
import io
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest
import zipfile

from packaging.version import Version


SCRIPT = Path(__file__).resolve().parents[2] / "script" / "check_release_version.py"
SPEC = importlib.util.spec_from_file_location("check_release_version", SCRIPT)
release = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(release)


class ReleaseVersionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_versions(self, core="3.0.0", python="3.0.0"):
        (self.root / "ALPS_VERSION.txt").write_text(core + "\n")
        (self.root / "pyproject.toml").write_text(
            f'[project]\nname = "pyalps"\nversion = "{python}"\n'
        )

    def test_final_release(self):
        self.write_versions()
        self.assertEqual(
            release.check_version(self.root, "refs/tags/v3.0.0"), Version("3.0.0")
        )

    def test_v3_release_regression(self):
        self.write_versions("2.3.4", "2.3.4b1")
        with self.assertRaisesRegex(ValueError, "Release tag v3.0.0 disagrees"):
            release.check_version(self.root, "refs/tags/v3.0.0")

    def test_ci_rejects_original_release_using_github_ref(self):
        self.write_versions("2.3.4", "2.3.4b1")
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--root", str(self.root)],
            env={**os.environ, "GITHUB_REF": "refs/tags/v3.0.0"},
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 1)
        self.assertIn("Release tag v3.0.0 disagrees", result.stderr)

    def test_sdk_and_python_must_agree_even_on_branches(self):
        for ref in ("", "refs/heads/master", "refs/pull/142/merge", "refs/tags/v3.0.0"):
            with self.subTest(ref=ref):
                self.write_versions("2.3.4", "3.0.0")
                with self.assertRaisesRegex(ValueError, "disagrees with ALPS_VERSION.txt"):
                    release.check_version(self.root, ref)

    def test_non_tag_refs_are_not_release_tags(self):
        self.write_versions()
        for ref in ("", "refs/heads/master", "refs/heads/v4.0.0", "refs/pull/142/merge"):
            with self.subTest(ref=ref):
                self.assertEqual(release.check_version(self.root, ref), Version("3.0.0"))

    def test_prerelease_tags_match_pep440_versions(self):
        for label, suffix in (
            ("alpha.1", "a1"), ("beta.2", "b2"), ("rc.3", "rc3"), ("dev.4", ".dev4")
        ):
            with self.subTest(label=label):
                self.write_versions(python="3.0.0" + suffix)
                self.assertEqual(
                    release.check_version(self.root, "refs/tags/v3.0.0-" + label),
                    Version("3.0.0" + suffix),
                )

    def test_prerelease_cannot_be_published_as_final_or_different_prerelease(self):
        for python, tag in (
            ("3.0.0b1", "v3.0.0"),
            ("3.0.0", "v3.0.0-beta.1"),
            ("3.0.0b1", "v3.0.0-beta.2"),
        ):
            with self.subTest(python=python, tag=tag):
                self.write_versions(python=python)
                with self.assertRaisesRegex(ValueError, "disagrees with pyproject.toml"):
                    release.check_version(self.root, "refs/tags/" + tag)

    def test_malformed_tags_fail_closed(self):
        self.write_versions()
        for tag in (
            "v3.0", "v3.0.0-beta", "v3.0.0-final", "v3.0.0_1", "v03.0.0", "v3.0.0+local"
        ):
            with self.subTest(tag=tag):
                with self.assertRaisesRegex(ValueError, "Invalid release tag"):
                    release.check_version(self.root, "refs/tags/" + tag)

    def test_invalid_numeric_core(self):
        for core in ("3.0", "v3.0.0", "3.0.0-beta.1", "3.0.0\n2.3.4", "03.0.0"):
            with self.subTest(core=core):
                self.write_versions(core=core)
                with self.assertRaisesRegex(ValueError, "must contain MAJOR.MINOR.PATCH"):
                    release.check_version(self.root, "")


class DistributionVersionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dist = Path(self.temp.name)

    def artifact(self, kind, version="3.0.0", metadata_version=None, name="pyalps"):
        metadata = (
            f"Metadata-Version: 2.1\nName: {name}\n"
            f"Version: {metadata_version or version}\n"
        ).encode()
        if kind == "wheel":
            path = self.dist / f"{name}-{version}-cp313-cp313-manylinux_2_28_x86_64.whl"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr(f"{name}-{version}.dist-info/METADATA", metadata)
        else:
            path = self.dist / f"{name}-{version}.tar.gz"
            with tarfile.open(path, "w:gz") as archive:
                member = tarfile.TarInfo(f"{name}-{version}/PKG-INFO")
                member.size = len(metadata)
                archive.addfile(member, io.BytesIO(metadata))
        return path

    def test_matching_wheel_and_sdist(self):
        for version in ("3.0.0", "3.0.0b1"):
            with self.subTest(version=version):
                paths = [self.artifact(kind, version) for kind in ("wheel", "sdist")]
                release.check_distributions(self.dist, Version(version))
                for path in paths:
                    path.unlink()

    def test_stale_artifact_rejects_the_entire_batch(self):
        self.artifact("wheel")
        self.artifact("sdist", "2.3.4b1")
        with self.assertRaisesRegex(ValueError, "expected a pyalps 3.0.0"):
            release.check_distributions(self.dist, Version("3.0.0"))

    def test_renaming_an_old_artifact_does_not_fix_its_version(self):
        for kind in ("wheel", "sdist"):
            with self.subTest(kind=kind):
                path = self.artifact(kind, metadata_version="2.3.4b1")
                with self.assertRaisesRegex(
                    ValueError, "metadata does not describe pyalps 3.0.0"
                ):
                    release.check_distributions(self.dist, Version("3.0.0"))
                path.unlink()

    def test_other_project_is_rejected(self):
        self.artifact("wheel", name="other")
        with self.assertRaisesRegex(ValueError, "expected a pyalps 3.0.0"):
            release.check_distributions(self.dist, Version("3.0.0"))

    def test_empty_dist_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "No distributions"):
            release.check_distributions(self.dist, Version("3.0.0"))

    def test_unexpected_file_is_rejected(self):
        (self.dist / "README.txt").write_text("not a distribution")
        with self.assertRaisesRegex(ValueError, "Unexpected distribution"):
            release.check_distributions(self.dist, Version("3.0.0"))


if __name__ == "__main__":
    unittest.main()
