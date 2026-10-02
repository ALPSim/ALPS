"""Saved citation reader/export contract; no compiled Python bindings required."""
import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[2]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


generator = module("citation_generator", ROOT / "script/generate_citations.py")
citations = module("saved_citations", ROOT / "python/pyalps/src/pyalps/citations.py")
BASE = "/provenance/alps/citations"


class Archive:
    """ALPS archive read interface, including its null-dataspace empty vectors."""
    def __init__(self, records=()):
        self.data = {BASE + "/schema_version": 1}
        for record in records:
            path = BASE + "/records/" + record["id"]
            self.data[path + "/complete"] = 1
            self.data.update({path + "/" + key: copy.deepcopy(value)
                              for key, value in record.items() if key != "id"})

    def __getitem__(self, path):
        return self.data[path]

    def is_group(self, path):
        return any(key.startswith(path + "/") for key in self.data)

    def is_data(self, path):
        return path in self.data

    def is_null(self, path):
        return self.data[path] == []

    def list_children(self, path):
        return sorted({key[len(path) + 1:].split("/")[0] for key in self.data
                       if key.startswith(path + "/")})


class SavedCitationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.catalog = generator.load_catalog(ROOT)

    def snapshot(self, component="looper", activity="calculation", version="test-build"):
        cff, policy, references, framework = self.catalog
        return generator.snapshot(cff, policy, references, framework, component, version, activity)

    def test_roundtrip_and_read_only_history(self):
        old = self.snapshot(version="historical-build")
        new = self.snapshot(activity="analysis", version="new-build")
        source = Archive([old, new])
        before = copy.deepcopy(source.data)
        self.assertEqual(citations.read_citations(source), sorted([old, new], key=lambda s: s["id"]))
        self.assertEqual(source.data, before)

    def test_legacy_empty_roles_and_incomplete_append(self):
        legacy = Archive()
        legacy.data = {"/parameters/SEED": 17}
        self.assertEqual(citations.read_citations(legacy), [])
        record = self.snapshot("framework", "unspecified")
        source = Archive([record])
        self.assertEqual(citations.read_citations(source), [record])
        source.data[BASE + "/records/partial/component"] = "looper"
        self.assertEqual(citations.read_citations(source), [record])

    def test_corruption_and_unknown_schema_are_rejected(self):
        record = self.snapshot()
        for field in ("notice", "catalog_sha256", "software_version", "algorithm"):
            source = Archive([record])
            path = BASE + "/records/" + record["id"] + "/" + field
            source.data[path] = ["changed"] if field == "algorithm" else "changed"
            with self.subTest(field=field), self.assertRaises(ValueError):
                citations.read_citations(source)
        for data in ({BASE: "dataset"}, {BASE + "/schema_version": 2},
                     {BASE + "/schema_version": 1.0}, {BASE + "/schema_version": "1"},
                     {BASE + "/schema_version": 1, BASE + "/records": "dataset"}):
            source = Archive()
            source.data = data
            with self.assertRaises(ValueError):
                citations.read_citations(source)
        source = Archive([record])
        del source.data[BASE + "/records/" + record["id"] + "/notice"]
        with self.assertRaises(ValueError):
            citations.read_citations(source)

    def test_completion_marker_is_a_scalar_zero_or_one(self):
        record = self.snapshot()
        path = BASE + "/records/" + record["id"] + "/complete"
        for value in (2, -1, 1.0, "1", [1]):
            source = Archive([record])
            source.data[path] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                citations.read_citations(source)
        source.data[path] = 0
        self.assertEqual(citations.read_citations(source), [])

    def test_export_cff_requires_unambiguous_selection_and_preserves_history(self):
        calculation = self.snapshot()
        analysis = self.snapshot(activity="analysis")
        source = Archive([calculation, analysis])
        with tempfile.TemporaryDirectory(prefix="alps-cff-export-") as directory:
            path = Path(directory) / "CITATION.cff"
            self.assertEqual(citations.export_citations(source, path), path)
            self.assertEqual(path.read_text(encoding="utf-8"), calculation["bibliography_cff"])
            generator.validate_schema(yaml.safe_load(path.read_text(encoding="utf-8")), "cff-1.2.0.schema.json")
            with self.assertRaises(FileExistsError):
                citations.export_citations(source, path)
            source = Archive([calculation, self.snapshot("interaction")])
            with self.assertRaisesRegex(ValueError, "select a snapshot_id"):
                citations.export_citations(source, Path(directory) / "ambiguous.cff")
            with self.assertRaises(ValueError):
                citations.export_citations(source, Path(directory) / "missing.cff", "missing")
            selected = Path(directory) / "selected.cff"
            citations.export_citations(source, selected, calculation["id"])
            self.assertEqual(selected.read_text(encoding="utf-8"), calculation["bibliography_cff"])
            with self.assertRaises(ValueError):
                citations.export_citations(Archive(), Path(directory) / "empty.cff")


if __name__ == "__main__":
    unittest.main()
