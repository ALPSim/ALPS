"""Policy, generation, and packaging regressions without a full ALPS build."""

import copy
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import jsonschema
import yaml

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("generate_citations", ROOT / "script/generate_citations.py")
generator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(generator)


class CitationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="alps-citations-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for filename in ("CITATION.cff", "CITATIONS.yaml"):
            shutil.copyfile(ROOT / filename, self.root / filename)
        self.cff, self.policy, self.references, self.framework = generator.load_catalog(self.root)

    def save(self):
        for filename, value in (("CITATION.cff", self.cff), ("CITATIONS.yaml", self.policy)):
            (self.root / filename).write_text(yaml.safe_dump(value, allow_unicode=True), encoding="utf-8")

    def test_paper_mappings_and_solver_separation(self):
        dmrg = generator.select_references(self.policy, self.framework, "dmrg")
        self.assertEqual(dmrg["algorithm"], ["white1992", "white1993", "schollwock2005", "hallberg2006"])
        self.assertEqual(dmrg["implementation"], ["feiguin2013", "troyer1998scheduler"])
        qwl = generator.select_references(self.policy, self.framework, "qwl")
        self.assertEqual(qwl["algorithm"], ["wang2001prl", "wang2001pre", "troyer2003qwl"])
        self.assertIn(self.framework, qwl["implementation"])
        self.assertEqual(qwl["framework"], [self.framework])
        expected = {"interaction": ["rubtsov2005", "gull2011"],
                    "hybridization": ["werner2006", "gull2011"],
                    "hirschfye": ["hirschfye1986"]}
        for component, algorithms in expected.items():
            self.assertEqual(generator.select_references(self.policy, self.framework, component)["algorithm"], algorithms)

    def test_framework_deduplicated_but_both_roles_retained(self):
        text = generator.notice(self.policy, self.references, self.framework, "qwl")
        self.assertEqual(" ".join(text.split()).count(self.references[self.framework]["title"]), 1)
        self.assertIn("Implementation: [4], [5]", text)
        self.assertIn("Framework: [4]", text)
        self.assertNotIn("P05001", text)

    def test_valid_cff_rejects_custom_policy_fields(self):
        self.cff["components"] = self.policy["components"]
        self.save()
        with self.assertRaises(jsonschema.ValidationError):
            generator.load_catalog(self.root)

    def test_policy_schema_rejects_unknown_fields_and_missing_roles(self):
        for mutation in (lambda entry: entry.update(algoritm=[]), lambda entry: entry.pop("algorithm")):
            with self.subTest(mutation=mutation):
                policy = copy.deepcopy(self.policy)
                mutation(policy["components"]["dmrg"])
                with self.assertRaises(jsonschema.ValidationError):
                    generator.validate_schema(policy, "policy.schema.json")

    def test_unknown_reference_rejected(self):
        self.policy["components"]["dmrg"]["algorithm"].append("typo1992")
        self.save()
        with self.assertRaisesRegex(ValueError, "unknown reference typo1992"):
            generator.load_catalog(self.root)

    def test_duplicate_reference_id_and_doi_rejected(self):
        reference = copy.deepcopy(self.cff["references"][0])
        reference["title"] += " (duplicate record)"
        self.cff["references"].append(reference)
        self.save()
        with self.assertRaisesRegex(ValueError, "Duplicate reference key"):
            generator.load_catalog(self.root)
        reference["identifiers"][0]["value"] = "different_key"
        self.save()
        with self.assertRaisesRegex(ValueError, "Duplicate DOI"):
            generator.load_catalog(self.root)

    def test_duplicate_yaml_key_rejected(self):
        with (self.root / "CITATIONS.yaml").open("a") as stream:
            stream.write("\nschema_version: 1\n")
        with self.assertRaisesRegex(ValueError, "Duplicate YAML key"):
            generator.load_catalog(self.root)

    def test_unknown_component_and_cycles_rejected(self):
        self.policy["components"]["dmrg"]["uses"] = ["missing"]
        self.save()
        with self.assertRaisesRegex(ValueError, "unknown component"):
            generator.load_catalog(self.root)
        self.policy["components"]["dmrg"]["uses"] = ["scheduler"]
        self.policy["components"]["scheduler"]["uses"] = ["dmrg"]
        self.save()
        with self.assertRaisesRegex(ValueError, "Cyclic component dependency"):
            generator.load_catalog(self.root)

    def test_framework_switch_comes_from_cff(self):
        preferred = self.cff["preferred-citation"]
        old = copy.deepcopy(preferred)
        preferred["identifiers"][0]["value"] = "future_release"
        preferred["title"] = "A future release"
        self.cff["references"].append(old)
        self.save()
        _, policy, references, framework = generator.load_catalog(self.root)
        self.assertEqual(framework, "future_release")
        self.assertIn("A future release", generator.notice(policy, references, framework, "looper"))

    def test_checked_in_document_is_current(self):
        expected = generator.markdown(self.policy, self.references, self.framework)
        self.assertEqual((ROOT / "CITATION.md").read_text(encoding="utf-8"), expected)

    def test_saved_snapshots_are_valid_cff_and_match_cli_selection(self):
        for component in self.policy["components"]:
            with self.subTest(component=component):
                record = generator.snapshot(self.cff, self.policy, self.references, self.framework,
                                            component, "3.0.0-test")
                bibliography = yaml.safe_load(record["bibliography_cff"])
                generator.validate_schema(bibliography, "cff-1.2.0.schema.json")
                available = {generator.reference_key(ref) for ref in
                             [bibliography["preferred-citation"]] + bibliography.get("references", [])}
                selected = {key for role in generator.ROLES for key in record[role]}
                self.assertEqual(available, selected)
                self.assertEqual(record["notice"], generator.notice(self.policy, self.references, self.framework, component))
                self.assertEqual(record["id"], generator.canonical_digest({k:v for k,v in record.items() if k != "id"}))

    def test_snapshot_fingerprint_tracks_policy_bibliography_version_and_activity(self):
        def snapshot(version="3.0.0", activity="calculation"):
            return generator.snapshot(self.cff, self.policy, self.references, self.framework, "looper", version, activity)
        original = snapshot()
        self.assertNotEqual(original["id"], snapshot("3.0.1")["id"])
        self.assertNotEqual(original["id"], snapshot(activity="analysis")["id"])
        reordered = dict(reversed(list(self.policy.items())))
        self.assertEqual(original["catalog_sha256"], generator.canonical_digest({"cff": self.cff, "policy": reordered}))
        self.cff["preferred-citation"]["title"] += " corrected metadata"
        self.assertNotEqual(original["catalog_sha256"], snapshot()["catalog_sha256"])
        self.assertNotEqual(original["id"], snapshot()["id"])

    @unittest.skipUnless(shutil.which("c++"), "C++ compiler unavailable")
    def test_generated_cpp_round_trip_escaping(self):
        self.references[self.framework]["title"] = 'Quotes " backslash \\ semicolon ; Unicode ö —'
        expected = generator.notice(self.policy, self.references, self.framework, "framework")
        source = self.root / "escaping.cpp"
        source.write_text(
            '#include <iostream>\n#include <string>\n'
            'struct citation_entry { const char* component; const char* text; };\n'
            + generator.cpp_data(self.policy, self.references, self.framework)
            + '\nint main() { for (auto e : citation_entries) if (std::string(e.component) == "framework") std::cout << e.text; }\n',
            encoding="utf-8")
        binary = self.root / "escaping"
        subprocess.run(["c++", "-std=c++17", str(source), "-o", str(binary)], check=True, capture_output=True)
        self.assertEqual(subprocess.check_output([str(binary)], text=True), expected)

    @unittest.skipUnless(shutil.which("cmake"), "CMake unavailable")
    def test_cmake_regeneration_and_all_install_layouts(self):
        # Use the real generation/install module without configuring ALPS's
        # numerical dependencies. Only the small citation inputs are copied.
        shutil.copytree(ROOT / "script/citations", self.root / "script/citations")
        shutil.copyfile(ROOT / "script/generate_citations.py", self.root / "script/generate_citations.py")
        project = self.root / "CMakeLists.txt"
        project.write_text(
            'cmake_minimum_required(VERSION 3.22)\nproject(citation_test LANGUAGES NONE)\n'
            f'include("{ROOT.as_posix()}/cmake/ALPSCitations.cmake")\n', encoding="utf-8")
        build = self.root / "build"
        for mode, destination in (("native", "share/alps"), ("libraries", "share/alps"), ("wheel", "pyalps/share/alps")):
            with self.subTest(mode=mode):
                subprocess.run(["cmake", "-S", str(self.root), "-B", str(build),
                                f"-DALPS_CITATION_PYTHON={sys.executable}",
                                "-DALPS_BUILD_PYTHON=OFF",
                                f"-DALPS_BUILD_LIBS_ONLY={'ON' if mode == 'libraries' else 'OFF'}",
                                f"-DALPS_PYTHON_WHEEL={'ON' if mode == 'wheel' else 'OFF'}"],
                               check=True, capture_output=True)
                prefix = self.root / ("install-" + mode)
                subprocess.run(["cmake", "--install", str(build), "--prefix", str(prefix), "--component", "libraries"],
                               check=True, capture_output=True)
                for filename in ("CITATION.cff", "CITATIONS.yaml", "CITATION.md"):
                    self.assertTrue((prefix / destination / filename).is_file())
        data = build / "src/alps/utility/citations_data.inc"
        before = data.read_text(encoding="utf-8")
        self.cff["preferred-citation"]["title"] = "Changed using only the catalog"
        self.save()
        subprocess.run(["cmake", "--build", str(build)], check=True, capture_output=True)
        after = data.read_text(encoding="utf-8")
        self.assertNotEqual(before, after)
        self.assertIn("Changed using only the catalog", after)
        self.policy["components"]["dmrg"]["algorithm"] = ["missing_reference"]
        self.save()
        invalid = subprocess.run(["cmake", "--build", str(build)], capture_output=True, text=True)
        self.assertNotEqual(invalid.returncode, 0)
        self.assertIn("unknown reference missing_reference", " ".join((invalid.stdout + invalid.stderr).split()))


if __name__ == "__main__":
    unittest.main()
