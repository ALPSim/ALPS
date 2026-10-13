"""Check public MaxEnt CLI exit and diagnostic behavior using the real binary."""

import argparse
from pathlib import Path
import subprocess
import tempfile
import unittest


class MaxEntCLIContract(unittest.TestCase):
    executable: Path

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="alps-maxent-cli-")
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)

    def run_cli(self, *arguments):
        before = set(self.directory.iterdir())
        result = subprocess.run(
            [str(self.executable), *map(str, arguments)],
            cwd=self.directory,
            capture_output=True,
            encoding="utf-8",
            errors="replace",
            timeout=30,
        )
        self.assertEqual(
            set(self.directory.iterdir()), before,
            "Help and rejected input must not create scientific output or crash artifacts",
        )
        return result

    def assert_help(self, result):
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--input-file", result.stdout)
        self.assertIn("--output-file", result.stdout)
        self.assertEqual(result.stderr, "")

    def assert_failure(self, result, diagnostic):
        self.assertEqual(
            result.returncode, 1,
            f"Expected a reported failure, not an abort: {result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )
        self.assertIn("maxent:", result.stderr)
        self.assertIn(diagnostic, result.stderr)

    def test_help_without_input(self):
        self.assert_help(self.run_cli("--help"))

    def test_information_queries_reject_input_before_loading(self):
        missing = self.directory / "not present.in.h5"
        for query in ("--help", "--license", "--citations"):
            with self.subTest(query=query):
                result = self.run_cli(query, "--input-file", missing)
                self.assert_failure(result, "without calculation arguments")
                self.assertNotIn("Recommended citations for ", result.stdout)

    def test_missing_input(self):
        self.assert_failure(self.run_cli(), "No job file specified")

    def test_unknown_option(self):
        option = "--not-a-maxent-option"
        self.assert_failure(self.run_cli(option), option)

    def test_missing_file(self):
        missing = self.directory / "not present.in.h5"
        self.assert_failure(self.run_cli(missing), missing.name)

    def test_invalid_archive(self):
        invalid = self.directory / "invalid.in.h5"
        invalid.write_text("This is not an HDF5 archive.\n", encoding="utf-8")
        self.assert_failure(self.run_cli("--input-file", invalid), invalid.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("maxent", type=Path, help="path to the MaxEnt executable")
    arguments = parser.parse_args()
    MaxEntCLIContract.executable = arguments.maxent.resolve(strict=True)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(MaxEntCLIContract)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
