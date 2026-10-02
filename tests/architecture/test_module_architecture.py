"""Exercise module boundaries independently of ALPS's current directory layout."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[2] / ".github/scripts/check_module_architecture.py"
SPEC = importlib.util.spec_from_file_location("module_architecture", SCRIPT)
architecture = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(architecture)


class SourceTree:
    def __init__(self, root):
        self.root = root
        self.manifest = {"schema_version": 1, "source_root": str(root), "modules": []}

    def write(self, path, content=""):
        destination = self.root / path
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(content)
        return path

    def module(self, name, *, depends=(), header=None, source=None):
        root = f"src/alps/{name}"
        module = {
            "name": name, "root": root, "public_include_roots": [f"{root}/include"],
            "public_headers": [], "private_files": [], "test_roots": [f"{root}/tests"],
            "depends": list(depends), "generated_headers": [],
        }
        self.manifest["modules"].append(module)
        if header is not None:
            module["public_headers"].append(self.write(f"{root}/include/alps/{name}.hpp", header))
        if source is not None:
            module["private_files"].append(self.write(f"{root}/src/{name}.cpp", source))
        return module

    def audit(self):
        return architecture.audit_manifest(self.manifest)


def test_declared_public_and_private_dependencies_and_generated_headers(tmp_path):
    tree = SourceTree(tmp_path)
    foundation = tree.module("foundation", header="#include <vector>\n")
    foundation["generated_headers"] = ["alps/config.h"]
    solver = tree.module("solver", depends=["foundation"],
                         header="#include <alps/foundation.hpp>\n#include <alps/config.h>\n",
                         source='#include <alps/solver.hpp>\n#include "private.hpp"\n')
    solver["private_files"].append(tree.write("src/alps/solver/src/private.hpp", "#include <boost/variant.hpp>\n"))
    tree.write("src/alps/solver/tests/helper.hpp", "#include <alps/test_only.hpp>\n")
    tree.write("src/alps/solver/CMakeLists.txt", "module setup")
    tree.write("src/alps/foundation/config.h.in", "template")
    result = tree.audit()
    assert result["errors"] == []
    assert result["public_header_owners"]["alps/config.h"] == "foundation"
    assert result["modules"]["solver"]["observed_dependencies"] == ["foundation"]
    assert result["observed_cycles"] == []


def test_duplicate_public_spelling_is_rejected(tmp_path):
    tree = SourceTree(tmp_path)
    tree.module("one", header="")
    second = tree.module("two", source="")
    second["public_headers"] = [tree.write("src/alps/two/include/alps/one.hpp")]
    assert any("Duplicate public include alps/one.hpp" in error for error in tree.audit()["errors"])


def test_each_production_file_has_exactly_one_owner(tmp_path):
    tree = SourceTree(tmp_path)
    one = tree.module("one", header="", source="")
    one["private_files"].append(one["private_files"][0])
    tree.write("src/alps/one/src/forgotten.cpp")
    tree.write("src/alps/one/include/alps/unregistered.hpp")
    tree.write("src/alps/stranded.cpp")
    errors = tree.audit()["errors"]
    assert any("File owned more than once" in error for error in errors)
    assert sum("Unowned production file" in error for error in errors) == 3


def test_dependency_in_disabled_conditional_branch_is_still_required(tmp_path):
    tree = SourceTree(tmp_path)
    tree.module("one", header="#if defined(OTHER_PLATFORM)\n#include <alps/two.hpp>\n#endif\n")
    tree.module("two", header="")
    errors = tree.audit()["errors"]
    assert any("one -> two" in error and ":2:" in error for error in errors)


def test_comments_and_raw_string_contents_do_not_create_dependencies(tmp_path):
    tree = SourceTree(tmp_path)
    tree.module("one", header='''// #include <alps/missing.hpp>
/*
#include <alps/also_missing.hpp>
*/
const char* sample = R"sample(
#include <alps/not_a_directive.hpp>
)sample";
#include <vector> // external dependency
''')
    assert tree.audit()["errors"] == []


def test_physical_include_into_another_module_is_rejected_even_when_declared(tmp_path):
    tree = SourceTree(tmp_path)
    tree.module("one", source='#include "../../two/include/alps/two.hpp"\n', depends=["two"])
    tree.module("two", header="")
    errors = tree.audit()["errors"]
    assert any("physical include crosses from one into two" in error for error in errors)
    assert not any("undeclared module dependency" in error for error in errors)


def test_public_header_cannot_depend_on_own_private_header(tmp_path):
    tree = SourceTree(tmp_path)
    module = tree.module("one", header='#include "../../src/private.hpp"\n')
    module["private_files"].append(tree.write("src/alps/one/src/private.hpp"))
    assert any("public header includes a private file" in error for error in tree.audit()["errors"])


@pytest.mark.parametrize("spelling", ["alps/deprecated.hpp", "ietl/deprecated.h"])
def test_exact_documented_unresolved_exemption_is_reported_as_debt(tmp_path, spelling):
    tree = SourceTree(tmp_path)
    module = tree.module("one", header=f"#include <{spelling}>\n")
    module["allowed_unresolved"] = [{"file": module["public_headers"][0], "include": spelling,
                                     "reason": "Shipped historical adapter; upstream backend was removed."}]
    result = tree.audit()
    assert result["errors"] == []
    assert result["counts"]["known_unresolved"] == 1
    assert result["known_debt"][0]["include"] == spelling


@pytest.mark.parametrize("change", [
    {"file": "src/alps/one/include/alps/*.hpp"},
    {"include": "alps/*"},
    {"reason": "  "},
    {"include": "boost/unknown.hpp"},
])
def test_exemptions_cannot_be_broad_or_unjustified(tmp_path, change):
    tree = SourceTree(tmp_path)
    module = tree.module("one", header="#include <alps/missing.hpp>\n")
    exemption = {"file": module["public_headers"][0], "include": "alps/missing.hpp", "reason": "Legacy debt."}
    module["allowed_unresolved"] = [exemption | change]
    errors = tree.audit()["errors"]
    assert any("exemption needs exact file/include and a reason" in error for error in errors)
    assert any("unresolved first-party include" in error for error in errors)


def test_exemption_does_not_hide_another_file_or_become_permanent(tmp_path):
    tree = SourceTree(tmp_path)
    module = tree.module("one", header="#include <alps/missing.hpp>\n",
                         source="#include <alps/missing.hpp>\n")
    module["allowed_unresolved"] = [{"file": module["public_headers"][0], "include": "alps/missing.hpp",
                                     "reason": "Legacy debt."}]
    result = tree.audit()
    assert len(result["errors"]) == 1
    assert "/src/one.cpp" in result["errors"][0]
    tree.module("missing", header="")
    module["depends"] = ["missing"]
    assert any("Unused unresolved-include exemption" in error for error in tree.audit()["errors"])


def test_cycles_are_recorded_and_unknown_module_dependencies_rejected(tmp_path):
    tree = SourceTree(tmp_path)
    one = tree.module("one", depends=["two"], header="#include <alps/two.hpp>\n")
    tree.module("two", depends=["one"], header="#include <alps/one.hpp>\n")
    result = tree.audit()
    assert result["errors"] == []
    assert result["observed_cycles"] == [["one", "two"]]
    assert result["declared_cycles"] == [["one", "two"]]
    one["depends"].append("typo")
    assert any("unknown declared dependency: typo" in error for error in tree.audit()["errors"])


def test_generated_header_dependencies_have_an_owner(tmp_path):
    tree = SourceTree(tmp_path)
    one = tree.module("one", header="")
    one["generated_headers"] = ["alps/config.h"]
    tree.module("two", header="#include <alps/config.h>\n")
    assert any("two -> one" in error for error in tree.audit()["errors"])


def test_explicit_bundled_dependency_sources_and_headers_can_belong_to_module(tmp_path):
    tree = SourceTree(tmp_path)
    module = tree.module("osiris", header="#include <alps/osiris/xdrcore.h>\n")
    module["public_include_roots"].append("third_party/xdr/include")
    module["public_headers"].append(tree.write("third_party/xdr/include/alps/osiris/xdrcore.h"))
    module["private_files"].append(tree.write("third_party/xdr/src/xdr.c", "#include <alps/osiris/xdrcore.h>\n"))
    assert tree.audit()["errors"] == []
    tree.write("third_party/xdr/include/alps/osiris/stranded.h")
    assert any("Unowned production file" in error for error in tree.audit()["errors"])


def test_nonliteral_include_is_visible_for_manual_review(tmp_path):
    tree = SourceTree(tmp_path)
    tree.module("one", header="#include ALPS_BACKEND_HEADER\n")
    result = tree.audit()
    assert result["errors"] == []
    assert result["nonliteral_includes"][0]["expression"] == "ALPS_BACKEND_HEADER"


def test_cli_writes_report_and_fails_for_stranded_source(tmp_path):
    tree = SourceTree(tmp_path)
    tree.module("one", header="")
    tree.write("src/alps/unowned.cpp")
    manifest = tmp_path / "modules.json"
    manifest.write_text(json.dumps(tree.manifest))
    report = tmp_path / "report.json"
    result = subprocess.run([sys.executable, str(SCRIPT), "--manifest", str(manifest),
                             "--write-report", str(report)], text=True, capture_output=True)
    assert result.returncode == 1
    assert "Unowned production file" in result.stderr
    assert json.loads(report.read_text())["counts"]["modules"] == 1
