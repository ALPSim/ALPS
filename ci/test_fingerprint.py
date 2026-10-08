"""Unit tests for ci/fingerprint.py. Run: python3 -m unittest ci/test_fingerprint.py"""

import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import fingerprint as fp  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def same(path, a, b, strip=True):
    return fp.normalize(path, a.encode(), strip) == fp.normalize(path, b.encode(), strip)


class CFamily(unittest.TestCase):
    def test_comment_only_edits_are_equal(self):
        self.assertTrue(same("a.cpp", "int x; // one\n", "int x; // two\n"))
        self.assertTrue(same("a.cpp", "/* one */ int x;\n", "/* two\n lines */ int x;\n"))
        self.assertTrue(same("a.ts", "let x = 1; // a\n", "let x = 1;\n"))
        self.assertTrue(same("a.C", "// header\nint f();\n", "int f();\n"))

    def test_whitespace_reflow_is_equal(self):
        self.assertTrue(same("a.cpp", "int  f(int a,  int b) {\n  return a+b;\n}\n",
                             "int f(int a, int b) {\n\n        return a+b;   \n}"))

    def test_wrapping_code_in_block_comment_is_different(self):
        before = "int a;\nint b;\n"
        after = "/*\nint a;\n*/\nint b;\n"
        self.assertFalse(same("a.c", before, after))

    def test_commenting_out_a_line_is_different(self):
        self.assertFalse(same("a.cpp", "x();\ny();\n", "// x();\ny();\n"))

    def test_string_contents_are_preserved(self):
        self.assertTrue(same("a.cpp", 's = "http://x";\n', 's = "http://x"; // c\n'))
        self.assertFalse(same("a.cpp", 's = "http://x";\n', 's = "http://y";\n'))
        self.assertFalse(same("a.js", "s = '/* not a comment */';\n", "s = '/* changed */';\n"))
        self.assertFalse(same("a.cpp", 's = "a  b";\n', 's = "a b";\n'))

    def test_pointer_dereference_is_code(self):
        self.assertFalse(same("a.c", "*ptr = 5;\n", "*ptr = 6;\n"))
        self.assertIn(b"*ptr = 5;", fp.normalize("a.c", b"*ptr = 5; /* c */\n", True))

    def test_directive_comments_are_code(self):
        self.assertFalse(same("a.ts", "// eslint-disable-next-line no-x\nf();\n",
                              "// eslint-disable-next-line no-y\nf();\n"))
        self.assertFalse(same("a.go", "//go:build linux\npackage a\n",
                              "//go:build darwin\npackage a\n"))
        self.assertFalse(same("a.cpp", "f(); // NOLINT\n", "f();\n"))

    def test_unterminated_input_falls_back_to_raw(self):
        for text in ('s = "abc\n', "/* open\nint x;\n", "s = `tmpl\n"):
            self.assertEqual(fp.normalize("a.js", text.encode(), True), text.encode())

    def test_backslash_continues_cpp_line_comment(self):
        self.assertFalse(same("a.cpp", "// c\nint x;\n", "// c \\\nint x;\n"))

    def test_preprocessor_lines_keep_their_newlines(self):
        self.assertFalse(same("a.h", "#define A 1\nint b;\n", "#define A 1 int b;\n"))

    def test_cpp_raw_strings_and_digit_separators(self):
        self.assertFalse(same("a.cpp", 's = R"x(// a)x";\n', 's = R"x(// b)x";\n'))
        self.assertTrue(same("a.cpp", "n = 1'000; // a\n", "n = 1'000;\n"))

    def test_js_regex_and_templates(self):
        self.assertFalse(same("a.js", "r = /a\\/\\/b/;\n", "r = /a\\/\\/c/;\n"))
        self.assertTrue(same("a.js", "r = a / b; // x\n", "r = a / b;\n"))
        self.assertFalse(same("a.ts", "t = `// ${x} a`;\n", "t = `// ${x} b`;\n"))
        nested = "t = `${`x`}`;\n"
        self.assertEqual(fp.normalize("a.ts", nested.encode(), True), nested.encode())

    def test_rust_and_go_raw_strings(self):
        self.assertFalse(same("a.rs", 's = r#"// a"#;\n', 's = r#"// b"#;\n'))
        self.assertTrue(same("a.rs", "fn f<'a>(x: &'a str) {} // c\n", "fn f<'a>(x: &'a str) {}\n"))
        self.assertFalse(same("a.go", "s := `// a`\n", "s := `// b`\n"))

    def test_sql_line_comments(self):
        self.assertTrue(same("a.sql", "select 1; -- a\n", "select 1;\n"))
        self.assertFalse(same("a.sql", "select '--a';\n", "select '--b';\n"))


class Python(unittest.TestCase):
    def test_comment_edit_is_equal(self):
        self.assertTrue(same("a.py", "x = 1  # one\n\n# block\n", "x = 1  # two\n"))

    def test_docstring_edit_is_different(self):
        self.assertFalse(same("a.py", 'def f():\n    """One."""\n', 'def f():\n    """Two."""\n'))

    def test_indentation_change_is_different(self):
        self.assertFalse(same("a.py", "if x:\n    a()\nb()\n", "if x:\n    a()\n    b()\n"))

    def test_directive_comments_are_code(self):
        self.assertFalse(same("a.py", "import os  # noqa\n", "import os\n"))

    def test_invalid_python_falls_back_to_raw(self):
        self.assertEqual(fp.normalize("a.py", b"if x:\n  a(\n", True), b"if x:\n  a(\n")


class HashComments(unittest.TestCase):
    def test_yaml(self):
        self.assertTrue(same("a.yml", "a: 1  # x\n\n# y\nb: 2\n", "a: 1\nb: 2\n"))
        self.assertFalse(same("a.yml", "a:\n  b: 1\n", "a:\n    b: 1\n"))
        self.assertFalse(same("a.yml", "a: 'x # y'\n", "a: 'x # z'\n"))
        self.assertFalse(same("a.yml", "a: don't\nb: 'x # y'\n", "a: don't\nb: 'x # z'\n"))
        self.assertFalse(same("a.yml", "url: a#b\n", "url: a#c\n"))

    def test_yaml_block_scalars_are_verbatim(self):
        self.assertFalse(same("a.yml", "run: |\n  # a\n  echo\nnext: 1\n",
                              "run: |\n  # b\n  echo\nnext: 1\n"))
        self.assertTrue(same("a.yml", "run: |\n  echo\n# c\nnext: 1\n", "run: |\n  echo\nnext: 1\n"))

    def test_toml(self):
        self.assertTrue(same("a.toml", 'a = "x"  # c\n', 'a = "x"\n'))
        self.assertFalse(same("a.toml", 'a = "x # y"\n', 'a = "x # z"\n'))
        self.assertFalse(same("a.toml", 'a = """\n# y\n"""\n', 'a = """\n# z\n"""\n'))


class Settings(unittest.TestCase):
    def test_strip_comments_false_compares_raw(self):
        self.assertFalse(same("a.cpp", "int x; // a\n", "int x; // b\n", strip=False))

    def test_unknown_extensions_are_raw(self):
        self.assertFalse(same("a.sh", "echo # a\n", "echo # b\n"))


class Globs(unittest.TestCase):
    def test_double_star(self):
        m = fp.Matcher(["services/api/**"])
        self.assertTrue(m("services/api/a/b.ts"))
        self.assertTrue(m("services/api/x.ts"))
        self.assertFalse(m("services/apix/a.ts"))

    def test_single_star_and_question_mark(self):
        self.assertTrue(fp.Matcher(["*.md"])("README.md"))
        self.assertFalse(fp.Matcher(["*.md"])("docs/README.md"))
        self.assertTrue(fp.Matcher(["**/*.md"])("README.md"))
        self.assertTrue(fp.Matcher(["**/*.md"])("a/b/c.md"))
        self.assertTrue(fp.Matcher(["a?.txt"])("ab.txt"))
        self.assertFalse(fp.Matcher(["a?.txt"])("a/.txt"))
        self.assertFalse(fp.Matcher(["a.txt"])("abtxt"))


CONFIG = {
    "global": ["ci/**"],
    "ignore": ["docs/**", "**/*.md"],
    "areas": {
        "api": {"paths": ["services/api/**", "libs/shared/**"], "strip_comments": True},
        "web": {"paths": ["apps/web/**", "libs/shared/**"], "strip_comments": True},
        "lint": {"paths": ["**/*.ts"], "strip_comments": False},
    },
}


class Repository(unittest.TestCase):
    """Change detection and keys against a throwaway git repository."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        cwd = os.getcwd()
        os.chdir(self.tmp.name)
        self.addCleanup(os.chdir, cwd)
        self.git("init", "-q")
        self.write("services/api/a.ts", "export const a = 1;\n")
        self.write("apps/web/b.ts", "export const b = 2;\n")
        self.write("libs/shared/c.ts", "export const c = 3;\n")
        self.write("docs/guide.txt", "hello\n")
        self.config = fp.Config(json.dumps(CONFIG).encode())
        self.base = self.commit()

    def git(self, *args):
        env = dict(os.environ, GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t",
                   GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@t")
        return subprocess.run(("git",) + args, check=True, env=env, stdout=subprocess.PIPE,
                              universal_newlines=True).stdout.strip()

    def write(self, path, text):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w") as handle:
            handle.write(text)

    def commit(self):
        self.git("add", "-A")
        self.git("commit", "-q", "--allow-empty", "-m", "c")
        return self.git("rev-parse", "HEAD")

    def detect(self, base=None):
        blobs = fp.Blobs()
        try:
            return {a: v[0] for a, v in
                    fp.detect_changes(self.config, base or self.base, "HEAD", blobs).items()}
        finally:
            blobs.close()

    def keys(self):
        blobs = fp.Blobs()
        try:
            return fp.area_keys(self.config, "HEAD", blobs, [], "")
        finally:
            blobs.close()

    def test_comment_only_edit_affects_only_raw_areas(self):
        self.write("services/api/a.ts", "// note\nexport const a = 1;\n")
        self.commit()
        self.assertEqual(self.detect(), {"api": False, "web": False, "lint": True})

    def test_code_edit_affects_its_area(self):
        self.write("services/api/a.ts", "export const a = 2;\n")
        self.commit()
        self.assertEqual(self.detect(), {"api": True, "web": False, "lint": True})

    def test_shared_library_affects_dependants(self):
        self.write("libs/shared/c.ts", "export const c = 4;\n")
        self.commit()
        self.assertEqual(self.detect(), {"api": True, "web": True, "lint": True})

    def test_docs_affect_nothing_and_unknown_paths_affect_everything(self):
        self.write("docs/guide.txt", "changed\n")
        self.commit()
        self.assertEqual(self.detect(), {"api": False, "web": False, "lint": False})
        self.write("tools/new.sh", "echo\n")
        self.commit()
        self.assertEqual(self.detect(), {"api": True, "web": True, "lint": True})

    def test_adding_or_removing_ignored_files_affects_everything(self):
        self.write("docs/new.txt", "x\n")
        self.commit()
        self.assertTrue(all(self.detect().values()))
        before = self.keys()
        os.remove("docs/new.txt")
        self.commit()
        self.assertNotEqual(before["api"], self.keys()["api"])

    def test_global_and_missing_base_affect_everything(self):
        self.assertTrue(all(fp.detect_changes(self.config, "0" * 40, "HEAD", None)[a][0]
                            for a in CONFIG["areas"]))
        self.assertTrue(all(self.detect(base="f" * 40).values()))
        self.write("ci/x.json", "{}\n")
        self.commit()
        self.assertTrue(all(self.detect().values()))

    def test_comment_only_commit_keeps_strip_area_keys(self):
        before = self.keys()
        self.write("services/api/a.ts", "/* note */ export const a = 1;\n")
        self.commit()
        after = self.keys()
        self.assertEqual(before["api"], after["api"])
        self.assertNotEqual(before["lint"], after["lint"])

    def test_key_changes_with_content_and_reverts_back(self):
        before = self.keys()
        self.write("services/api/a.ts", "export const a = 9;\n")
        self.commit()
        self.assertNotEqual(before["api"], self.keys()["api"])
        self.assertEqual(before["web"], self.keys()["web"])
        self.git("revert", "--no-edit", "HEAD")
        self.assertEqual(before, self.keys())

    def test_key_is_independent_of_listing_order(self):
        entries = fp.tree_entries("HEAD")
        expected = self.keys()
        with mock.patch.object(fp, "tree_entries", return_value=list(reversed(entries))):
            self.assertEqual(expected, self.keys())

    def test_citation_inputs_invalidate_every_area(self):
        # Native configure checks raw checksums, including generated Markdown
        # and comments. Cached passes must not hide stale citation snapshots.
        self.config = fp.Config.load(os.path.join(ROOT, "ci", "areas.json"))
        for path in ("CITATION.cff", "CITATIONS.yaml", "CITATION.md",
                     ".github/scripts/generate_citations.py",
                     ".github/scripts/citations/policy.schema.json",
                     ".github/scripts/citations/generated/citations_data.inc",
                     ".github/scripts/citations/generated/snapshots.cmake"):
            with self.subTest(path=path):
                self.write(path, "# original\n")
                base = self.commit()
                before = self.keys()
                self.write(path, "# changed\n")
                self.commit()
                self.assertTrue(all(self.detect(base).values()))
                after = self.keys()
                for area in self.config.areas:
                    self.assertNotEqual(before[area], after[area], area)


class Workflow(unittest.TestCase):
    def test_ci_workflow_matches_areas(self):
        config = fp.Config.load(os.path.join(ROOT, "ci", "areas.json"))
        with open(os.path.join(ROOT, ".github", "workflows", "ci.yml")) as handle:
            self.assertEqual(fp.check_workflow(config, handle.read()), [])

    def test_areas_cover_the_inputs_their_jobs_read(self):
        # tests/ci/test_numeric_bindings.py (static) scans these trees, and
        # run_wheel_tests.py (packaging) runs tests/cmake against the wheel.
        config = fp.Config.load(os.path.join(ROOT, "ci", "areas.json"))
        for area, path in (("static", "src/alps/ietl/include/ietl/jd.h"),
                           ("static", "python/pyalps/cpp/pyalea.cpp"),
                           ("static", "tests/integration/CMakeLists.txt"),
                           ("static", "third_party/boost_numeric_bindings/README.md"),
                           ("packaging", "tests/cmake/test_sdk_contracts.py")):
            with self.subTest(area=area, path=path):
                self.assertTrue(config.in_footprint(area, path))

    def test_decide_requires_a_real_key_for_a_marker_hit(self):
        env = {"AFFECTED_API": "true", "KEY_API": "", "HIT_API": "true",
               "AFFECTED_WEB": "true", "KEY_WEB": "a" * 64, "HIT_WEB": "true",
               "AFFECTED_LINT": "false", "GITHUB_STEP_SUMMARY": ""}
        path = os.path.join(tempfile.mkdtemp(), "areas.json")
        with open(path, "w") as handle:
            json.dump(CONFIG, handle)
        out = subprocess.run([sys.executable, fp.__file__, "--decide", "--config", path],
                             env=dict(os.environ, **env), stdout=subprocess.PIPE,
                             universal_newlines=True, check=True).stdout.split()
        self.assertEqual(out, ["api_run=true", "web_run=false", "lint_run=false"])


if __name__ == "__main__":
    unittest.main()
