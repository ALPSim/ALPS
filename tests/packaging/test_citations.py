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
SPEC = importlib.util.spec_from_file_location("generate_citations", ROOT / ".github/scripts/generate_citations.py")
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
        for component in ("dmft", "interaction", "hybridization", "hirschfye"):
            implementations = ["gull2011cpc"]
            if component == "hybridization":
                implementations.append("hafermann2013")
            self.assertEqual(generator.select_references(self.policy, self.framework, component)["implementation"],
                             implementations)

    def test_framework_deduplicated_but_both_roles_retained(self):
        text = generator.detailed_notice(self.policy, self.references, self.framework, "qwl")
        self.assertEqual(" ".join(text.split()).count(self.references[self.framework]["title"]), 1)
        self.assertIn("Implementation: [1], [5]", text)
        self.assertIn("Framework: [1]", text)
        self.assertNotIn("P05001", text)

    def test_citation_expected_output(self):
        cases = [
            ({"title": "Cluster algorithm", "authors": [
                {"given-names": "R. H.", "family-names": "Swendsen"},
                {"given-names": "J.-S.", "family-names": "Wang"}],
              "journal": "Physical Review Letters", "volume": "58",
              "start": "86", "end": "88", "year": 1987, "doi": "10.1103/PhysRevLett.58.86"},
             'R. H. Swendsen, J.-S. Wang, "Cluster algorithm." Physical Review Letters 58, 86–88 (1987). https://doi.org/10.1103/PhysRevLett.58.86',
             "R. H. Swendsen, J.-S. Wang, Physical Review Letters 58, 86 (1987)."),
            ({"title": "Numerical methods", "authors": [{"given-names": "A. E.", "family-names": "Feiguin"}],
              "collection-title": "Strongly Correlated Systems: Numerical Methods",
              "volume": "176", "start": "31", "end": "64", "year": 2013},
             'A. E. Feiguin, "Numerical methods." Strongly Correlated Systems: Numerical Methods 176, 31–64 (2013).',
             "A. E. Feiguin, Strongly Correlated Systems: Numerical Methods 176, 31 (2013)."),
            ({"title": "A methods book", "authors": [{"given-names": "A.", "name-particle": "van", "family-names": "Dijk"}],
              "year": 2020},
             'A. van Dijk, "A methods book." (2020).',
             "A. van Dijk, A methods book (2020)."),
            ({"title": "Release manuscript", "authors": [
                {"given-names": "F.", "family-names": "Alet"},
                {"family-names": "Chen"}, {"family-names": "Feiguin"}, {"family-names": "Wolf"}],
              "status": "in-preparation"},
             'F. Alet et al., "Release manuscript." In preparation.',
             "F. Alet et al., Release manuscript; in preparation."),
            ({"title": "Release paper", "authors": [
                {"given-names": "F.", "family-names": "Alet"},
                {"family-names": "Chen"}, {"family-names": "Feiguin"}, {"family-names": "Wolf"}],
              "year": 2026, "status": "preprint", "doi": "10.48550/arXiv.2610.03884"},
             'F. Alet et al., "Release paper." arXiv:2610.03884 (2026). Preprint. https://doi.org/10.48550/arXiv.2610.03884',
             "F. Alet et al., Release paper, arXiv:2610.03884 (2026); preprint."),
            ({"title": "Undated paper", "authors": [{"name": "ALPS collaboration"}],
              "journal": "Physics Journal", "volume": "5", "start": "7"},
             'ALPS collaboration, "Undated paper." Physics Journal 5, 7.',
             "ALPS collaboration, Physics Journal 5, 7."),
            ({"title": "Three authors", "authors": [
                {"family-names": "Troyer"}, {"family-names": "Ammon"}, {"family-names": "Heeb"}]},
             'Troyer, Ammon, Heeb, "Three authors."',
             "Troyer, Ammon, Heeb, Three authors."),
        ]
        for reference, detailed, compact in cases:
            with self.subTest(title=reference["title"]):
                self.assertEqual(generator.citation(reference), detailed)
                self.assertEqual(generator.compact_citation(reference), compact)

    def test_wolf_initials(self):
        authors = self.cff["preferred-citation"]["authors"]
        wolf, = [author for author in authors if author.get("family-names") == "Wolf"]
        self.assertEqual(wolf["given-names"], "T. M. R.")
        self.assertEqual(generator.person_name(wolf), "T. M. R. Wolf")

    def test_startup_box_is_compact_and_framework_first(self):
        for component in self.policy["components"]:
            with self.subTest(component=component):
                text = generator.notice(self.policy, self.references, self.framework, component)
                lines = text.rstrip().splitlines()
                self.assertEqual(lines[0], "*" * 80)
                self.assertEqual(lines[-1], lines[0])
                self.assertTrue(all(len(line) == 80 for line in lines))
                # DMRG's seven entries need 17 lines with recognizable authors.
                self.assertLessEqual(len(lines), 17 if component == "dmrg" else 15)
                self.assertIn("Framework: [1]", text)
                content = " ".join(line[2:-2].strip() for line in lines[1:-1])
                self.assertIn("[1] F. Alet et al., " + self.references[self.framework]["title"], content)
                self.assertNotIn("https://doi.org/", text)
                self.assertIn("Full references: --citations", text)
                details = generator.detailed_notice(self.policy, self.references, self.framework, component)
                self.assertIn("Framework: [1]", details)
                self.assertIn("[1] F. Alet et al.", details)
                if component != "framework":
                    self.assertIn("https://doi.org/", details)

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
        preferred.pop("doi", None)  # a new release paper has its own DOI
        self.cff["references"].append(old)
        self.save()
        _, policy, references, framework = generator.load_catalog(self.root)
        self.assertEqual(framework, "future_release")
        self.assertIn("A future release", generator.notice(policy, references, framework, "looper"))

    def test_checked_in_document_is_current(self):
        expected = generator.markdown(self.policy, self.references, self.framework)
        self.assertEqual((ROOT / "CITATION.md").read_text(encoding="utf-8"), expected)

    def test_checked_in_native_data_is_current(self):
        for filename, content in generator.native_outputs(ROOT, self.policy, self.references, self.framework).items():
            with self.subTest(filename=filename):
                self.assertEqual((ROOT / filename).read_text(encoding="utf-8"), content)

    @unittest.skipUnless(shutil.which("c++"), "C++ compiler unavailable")
    def test_generated_cpp_round_trip_escaping(self):
        self.references[self.framework]["title"] = 'Quotes " backslash \\ semicolon ; Unicode ö —'
        expected = generator.notice(self.policy, self.references, self.framework, "framework")
        source = self.root / "escaping.cpp"
        source.write_text(
            '#include <iostream>\n#include <string>\n'
            'struct citation_entry { const char* component; const char* text; const char* details; };\n'
            + generator.cpp_data(self.policy, self.references, self.framework)
            + '\nint main() { for (auto e : citation_entries) if (std::string(e.component) == "framework") std::cout << e.text; }\n',
            encoding="utf-8")
        binary = self.root / "escaping"
        subprocess.run(["c++", "-std=c++17", str(source), "-o", str(binary)], check=True, capture_output=True)
        self.assertEqual(subprocess.check_output([str(binary)], text=True), expected)

    @unittest.skipUnless(shutil.which("cmake"), "CMake unavailable")
    def test_native_build_needs_no_python(self):
        # Exercise the real native module with an unusable Python path. The
        # generator runs only for the explicit maintainer regeneration below.
        shutil.copytree(ROOT / ".github/scripts/citations", self.root / ".github/scripts/citations")
        shutil.copyfile(ROOT / ".github/scripts/generate_citations.py", self.root / ".github/scripts/generate_citations.py")
        shutil.copyfile(ROOT / "CITATION.md", self.root / "CITATION.md")
        project = self.root / "CMakeLists.txt"
        project.write_text(
            'cmake_minimum_required(VERSION 3.22)\nproject(citation_test LANGUAGES NONE)\n'
            f'include("{ROOT.as_posix()}/cmake/ALPSCitations.cmake")\n', encoding="utf-8")
        build = self.root / "build"
        for mode in ("native", "libraries"):
            with self.subTest(mode=mode):
                subprocess.run(["cmake", "-S", str(self.root), "-B", str(build),
                                "-DALPS_CITATION_PYTHON=/nonexistent/python",
                                "-DCMAKE_DISABLE_FIND_PACKAGE_Python3=ON",
                                f"-DALPS_BUILD_LIBS_ONLY={'ON' if mode == 'libraries' else 'OFF'}"],
                               check=True, capture_output=True)
                expected = generator.cpp_data(self.policy, self.references, self.framework)
                self.assertEqual((build / "src/alps/utility/citations_data.inc").read_text(encoding="utf-8"), expected)
                prefix = self.root / ("install-" + mode)
                subprocess.run(["cmake", "--install", str(build), "--prefix", str(prefix), "--component", "libraries"],
                               check=True, capture_output=True)
                for filename in ("CITATION.cff", "CITATIONS.yaml", "CITATION.md"):
                    self.assertTrue((prefix / "share/alps" / filename).is_file())
        data = build / "src/alps/utility/citations_data.inc"
        before = data.read_text(encoding="utf-8")
        self.cff["preferred-citation"]["title"] = "Changed using only the catalog @DO_NOT_REPLACE@"
        self.save()
        self.cff, self.policy, self.references, self.framework = generator.load_catalog(self.root)
        stale = subprocess.run(["cmake", "--build", str(build)], capture_output=True, text=True)
        self.assertNotEqual(stale.returncode, 0)
        self.assertIn("Generated citation data is stale", stale.stdout + stale.stderr)
        subprocess.run([sys.executable, str(self.root / ".github/scripts/generate_citations.py"), "--regenerate"],
                       check=True, capture_output=True)
        subprocess.run(["cmake", "--build", str(build)], check=True, capture_output=True)
        after = data.read_text(encoding="utf-8")
        self.assertNotEqual(before, after)
        self.assertIn("Changed using only the catalog", after)
        self.assertIn("@DO_NOT_REPLACE@", after)
        self.assertEqual(after, generator.cpp_data(self.policy, self.references, self.framework))
        derived = self.root / ".github/scripts/citations/generated/citations_data.inc"
        with derived.open("a", encoding="utf-8") as output:
            output.write("\n// manually edited generated data\n")
        stale = subprocess.run(["cmake", "--build", str(build)], capture_output=True, text=True)
        self.assertNotEqual(stale.returncode, 0)
        self.assertIn("Generated citation data is stale", stale.stdout + stale.stderr)
        self.policy["components"]["dmrg"]["algorithm"] = ["missing_reference"]
        self.save()
        invalid = subprocess.run([sys.executable, str(self.root / ".github/scripts/generate_citations.py"), "--regenerate"],
                                 capture_output=True, text=True)
        self.assertNotEqual(invalid.returncode, 0)
        self.assertIn("unknown reference missing_reference", " ".join((invalid.stdout + invalid.stderr).split()))


if __name__ == "__main__":
    unittest.main()
