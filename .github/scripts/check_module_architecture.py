#!/usr/bin/env python3
"""Audit CMake-declared ALPS module ownership and literal include dependencies.

The manifest comes from the configured modules' actual file sets and source lists.
This is a source-boundary check, not a replacement for compilation: all conditional
branches are inspected, while nonliteral preprocessor includes are reported for
manual review. Declared cycles are recorded rather than silently forbidden.
"""

import argparse
from collections import defaultdict
import json
from pathlib import Path
import re
import sys


CODE_SUFFIXES = {".h", ".hh", ".hpp", ".hxx", ".ipp", ".inl", ".tpp", ".tcc",
                 ".c", ".cc", ".cpp", ".cxx", ".f", ".f90", ".f95"}
FIRST_PARTY_PREFIXES = ("alps/", "ietl/")
INCLUDE = re.compile(
    r'^[ \t]*#[ \t]*include(?:[ \t]|\\\r?\n)*(?:[<"]([^>"\n]+)[>"]|([^\n]+))', re.M)
COMMENTS_AND_LITERALS = re.compile(
    r'(?P<raw>(?:u8|u|U|L)?R"(?P<delimiter>[^\s\\()]*)\([\s\S]*?\)(?P=delimiter)")'
    r'|(?P<quoted>"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\')'
    r'|(?P<comment>//[^\n]*|/\*[\s\S]*?\*/)')


def includes(path):
    """Yield (line, literal spelling, macro expression), ignoring commented code."""
    def strip(match):
        if match.lastgroup in {"comment", "raw"}:
            return re.sub(r"[^\n]", " ", match.group())
        return match.group()

    source = COMMENTS_AND_LITERALS.sub(strip, path.read_text(errors="replace"))
    for match in INCLUDE.finditer(source):
        yield source.count("\n", 0, match.start()) + 1, match.group(1), match.group(2)


def cycles(graph):
    """Return deterministic strongly connected components containing a cycle."""
    index = 0
    indices, low, stack, active, components = {}, {}, [], set(), []

    def visit(node):
        nonlocal index
        indices[node] = low[node] = index
        index += 1
        stack.append(node)
        active.add(node)
        for other in sorted(graph[node]):
            if other not in indices:
                visit(other)
                low[node] = min(low[node], low[other])
            elif other in active:
                low[node] = min(low[node], indices[other])
        if low[node] == indices[node]:
            component = []
            while True:
                other = stack.pop()
                active.remove(other)
                component.append(other)
                if other == node:
                    break
            if len(component) > 1 or node in graph[node]:
                components.append(sorted(component))

    for node in sorted(graph):
        if node not in indices:
            visit(node)
    return sorted(components)


def audit_manifest(manifest):
    """Return a JSON-serializable report; every contract violation is in errors."""
    if manifest.get("schema_version") != 1:
        raise ValueError("Expected module manifest schema_version 1")
    root = Path(manifest["source_root"]).resolve()
    errors, debt, macros = [], [], []
    modules, ownership, public, file_kinds = {}, {}, {}, {}
    roots, test_roots, all_include_roots = {}, {}, set()
    exemptions, used_exemptions = {}, set()

    def error(message):
        errors.append(message)

    def source_path(value, context):
        if not isinstance(value, str) or not value or Path(value).is_absolute():
            error(f"{context}: expected a nonempty repository-relative path")
            return None
        path = (root / value).resolve()
        if not path.is_relative_to(root):
            error(f"{context}: path escapes the source root: {value}")
            return None
        return path

    def relative(path):
        return path.relative_to(root).as_posix()

    def add_public(spelling, name, path=None):
        if spelling in public:
            error(f"Duplicate public include {spelling}: {public[spelling][0]} and {name}")
        else:
            public[spelling] = (name, path)

    for module in manifest["modules"]:
        name = module["name"]
        if name in modules:
            error(f"Duplicate module: {name}")
            continue
        modules[name] = module
        module_root = source_path(module["root"], f"{name} root")
        if module_root is None:
            continue
        roots[name] = module_root
        if not module_root.is_dir():
            error(f"{name}: module directory does not exist: {module['root']}")
        include_roots = []
        for value in module["public_include_roots"]:
            path = source_path(value, f"{name} include root")
            if path is not None:
                include_roots.append(path)
                all_include_roots.add(path)
        test_roots[name] = []
        for value in module.get("test_roots", []):
            path = source_path(value, f"{name} test root")
            if path is not None:
                if not path.is_relative_to(module_root):
                    error(f"{name}: test root lies outside its module: {value}")
                else:
                    test_roots[name].append(path)
        for field in ("public_headers", "private_files"):
            for value in module[field]:
                path = source_path(value, f"{name} {field}")
                if path is None:
                    continue
                if path in ownership:
                    error(f"File owned more than once: {value} ({ownership[path]}, {name})")
                    continue
                ownership[path], file_kinds[path] = name, field
                if not path.is_file():
                    error(f"{name}: declared file does not exist: {value}")
                if any(path.is_relative_to(base) for base in test_roots[name]):
                    error(f"{name}: production file is also inside a test root: {value}")
                if field == "public_headers":
                    bases = [base for base in include_roots if path.is_relative_to(base)]
                    if len(bases) != 1:
                        error(f"{name}: public header needs exactly one include root: {value}")
                    else:
                        add_public(path.relative_to(bases[0]).as_posix(), name, path)
                elif any(path.is_relative_to(base) for base in include_roots):
                    error(f"{name}: private file is inside a public include root: {value}")
        for spelling in module.get("generated_headers", []):
            if (not isinstance(spelling, str) or not spelling
                    or Path(spelling).is_absolute() or ".." in Path(spelling).parts):
                error(f"{name}: invalid generated header spelling: {spelling!r}")
            else:
                add_public(spelling, name)
        for entry in module.get("allowed_unresolved", []):
            filename, spelling, reason = entry.get("file"), entry.get("include"), entry.get("reason")
            if (not isinstance(filename, str) or not isinstance(spelling, str)
                    or not isinstance(reason, str) or not reason.strip()
                    or any(c in filename + spelling for c in "*?[]")
                    or not spelling.startswith(FIRST_PARTY_PREFIXES)):
                error(f"{name}: unresolved-include exemption needs exact file/include and a reason: {entry}")
                continue
            path = source_path(filename, f"{name} unresolved-include exemption")
            if path is None:
                continue
            key = (path, spelling)
            if ownership.get(path) != name:
                error(f"{name}: unresolved-include exemption is not for its owned file: {filename}")
            elif key in exemptions:
                error(f"{name}: duplicate unresolved-include exemption: {filename}: {spelling}")
            else:
                exemptions[key] = reason.strip()

    for name, base in roots.items():
        for other, other_base in roots.items():
            if name != other and base.is_relative_to(other_base):
                error(f"Overlapping module roots: {name} and {other}")
    declared = {name: set(module["depends"]) for name, module in modules.items()}
    for name, dependencies in declared.items():
        for dependency in sorted(dependencies - modules.keys()):
            error(f"{name}: unknown declared dependency: {dependency}")
    # Include the module parent area so a file stranded between modules is visible.
    scan_roots = {root / "src/alps", *roots.values(), *all_include_roots}
    all_test_roots = [path for paths in test_roots.values() for path in paths]
    for path in sorted({path.resolve() for base in scan_roots if base.is_dir()
                        for path in base.rglob("*") if path.is_file()
                        and path.suffix.lower() in CODE_SUFFIXES}):
        if any(path.is_relative_to(base) for base in all_test_roots):
            continue
        if path not in ownership:
            error(f"Unowned production file: {relative(path)}")

    observed = {name: set() for name in modules}
    evidence = defaultdict(list)
    for path, owner in sorted(ownership.items()):
        if not path.is_file() or path.suffix.lower() not in CODE_SUFFIXES:
            continue
        for line, spelling, macro in includes(path):
            location = f"{relative(path)}:{line}"
            if macro:
                macros.append({"file": relative(path), "line": line, "expression": macro.strip()})
                continue
            destination = public.get(spelling)
            local = (path.parent / spelling).resolve()
            if destination is None and local in ownership:
                destination = (ownership[local], local)
            if destination is None:
                if spelling.startswith(FIRST_PARTY_PREFIXES):
                    key = (path, spelling)
                    if key in exemptions:
                        used_exemptions.add(key)
                        debt.append({"file": relative(path), "line": line, "include": spelling,
                                     "reason": exemptions[key]})
                    else:
                        error(f"{location}: unresolved first-party include <{spelling}>")
                elif local.is_file() and local.is_relative_to(root):
                    error(f"{location}: include resolves to an unowned source file: {spelling}")
                elif ".." in Path(spelling).parts or spelling.startswith(("src/", "/")):
                    error(f"{location}: unresolved physical source-path include: {spelling}")
                continue
            dependency, target = destination
            if owner != dependency:
                observed[owner].add(dependency)
                evidence[(owner, dependency)].append(
                    {"file": relative(path), "line": line, "include": spelling})
                if spelling not in public:
                    error(f"{location}: physical include crosses from {owner} into {dependency}: {spelling}")
                if dependency not in declared[owner]:
                    error(f"{location}: undeclared module dependency {owner} -> {dependency}: {spelling}")
            if (file_kinds[path] == "public_headers" and target is not None
                    and file_kinds[target] == "private_files"):
                error(f"{location}: public header includes a private file: {spelling}")
    for path, spelling in sorted(exemptions.keys() - used_exemptions):
        error(f"Unused unresolved-include exemption: {relative(path)}: {spelling}")
    safe_declared = {name: dependencies & modules.keys() for name, dependencies in declared.items()}
    return {
        "schema_version": 1,
        "counts": {"modules": len(modules), "public_headers": len(public),
                   "production_files": len(ownership), "known_unresolved": len(used_exemptions)},
        "modules": {name: {"declared_dependencies": sorted(declared[name]),
                           "observed_dependencies": sorted(observed[name])} for name in sorted(modules)},
        "public_header_owners": {spelling: owner for spelling, (owner, _) in sorted(public.items())},
        "edges": [{"from": a, "to": b, "includes": records}
                  for (a, b), records in sorted(evidence.items())],
        "observed_cycles": cycles(observed),
        "declared_cycles": cycles(safe_declared),
        "known_debt": debt,
        "nonliteral_includes": macros,
        "errors": sorted(set(errors)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--write-report", type=Path)
    args = parser.parse_args()
    report = audit_manifest(json.loads(args.manifest.read_text()))
    if args.write_report:
        args.write_report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    for message in report["errors"]:
        print(message, file=sys.stderr)
    print("Module architecture: " + ", ".join(f"{key}={value}" for key, value in report["counts"].items())
          + f", errors={len(report['errors'])}")
    return int(bool(report["errors"]))


if __name__ == "__main__":
    raise SystemExit(main())
