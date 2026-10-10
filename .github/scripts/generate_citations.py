#!/usr/bin/env python3
"""Validate the two citation authorities and generate CLI data/documentation.

Requires .github/scripts/citations/requirements.txt. No network access is performed.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import textwrap

try:
    import jsonschema
    import yaml
except ImportError as exc:
    raise SystemExit(
        "Citation generation requires PyYAML and jsonschema. Install with: "
        "python -m pip install -r .github/scripts/citations/requirements.txt"
    ) from exc


ROOT = Path(__file__).resolve().parents[2]
SCHEMAS = Path(__file__).resolve().parent / "citations"
ROLES = ("framework", "algorithm", "implementation")
KEY_DESCRIPTION = "ALPS reference key"


class UniqueLoader(yaml.SafeLoader):
    """Reject duplicate YAML keys instead of silently losing policy entries."""


def unique_mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if not isinstance(key, str):
            raise ValueError(f"Expected a string mapping key at {key_node.start_mark}")
        if key in result:
            raise ValueError(f"Duplicate YAML key {key!r} at {key_node.start_mark}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueLoader.add_constructor("tag:yaml.org,2002:map", unique_mapping)
# CFF dates are strings for JSON Schema validation, including unquoted YAML dates.
UniqueLoader.add_constructor("tag:yaml.org,2002:timestamp", lambda loader, node: loader.construct_scalar(node))


def read_yaml(path):
    return yaml.load(path.read_text(encoding="utf-8"), Loader=UniqueLoader)


def validate_schema(data, filename):
    schema = json.loads((SCHEMAS / filename).read_text(encoding="utf-8"))
    jsonschema.Draft7Validator(schema, format_checker=jsonschema.FormatChecker()).validate(data)


def reference_key(reference):
    keys = [item["value"] for item in reference.get("identifiers", [])
            if item["type"] == "other" and item.get("description") == KEY_DESCRIPTION]
    if len(keys) != 1 or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", keys[0]):
        raise ValueError(f"Reference {reference['title']!r} must have one {KEY_DESCRIPTION}")
    return keys[0]


def load_catalog(root):
    cff = read_yaml(root / "CITATION.cff")
    policy = read_yaml(root / "CITATIONS.yaml")
    validate_schema(cff, "cff-1.2.0.schema.json")
    validate_schema(policy, "policy.schema.json")
    if "preferred-citation" not in cff:
        raise ValueError("CITATION.cff must specify the framework preferred-citation")
    preferred = cff["preferred-citation"]
    references = {}
    dois = set()
    for reference in [preferred] + cff.get("references", []):
        key = reference_key(reference)
        if key in references:
            raise ValueError(f"Duplicate reference key: {key}")
        doi = reference.get("doi", "").lower()
        if doi and doi in dois:
            raise ValueError(f"Duplicate DOI: {doi}; use one bibliographic record")
        if doi:
            dois.add(doi)
        references[key] = reference
    components = policy["components"]
    for component, entry in components.items():
        for role in ("algorithm", "implementation"):
            for key in entry[role]:
                if key not in references:
                    raise ValueError(f"{component}.{role}: unknown reference {key}")
        for dependency in entry.get("uses", []):
            if dependency not in components:
                raise ValueError(f"{component}.uses: unknown component {dependency}")
    # Check all components, including those not currently used by a binary.
    for component in components:
        component_order(components, component)
    return cff, policy, references, reference_key(preferred)


def component_order(components, component, active=()):
    if component in active:
        raise ValueError("Cyclic component dependency: " + " -> ".join(active + (component,)))
    result = [component]
    for dependency in components[component].get("uses", []):
        for item in component_order(components, dependency, active + (component,)):
            if item not in result:
                result.append(item)
    return result


def select_references(policy, framework, component):
    roles = {role: [] for role in ROLES}
    for item in component_order(policy["components"], component):
        for role in ("algorithm", "implementation"):
            for key in policy["components"][item][role]:
                if key not in roles[role]:
                    roles[role].append(key)
    roles["framework"] = [framework]
    return roles


def person_name(person):
    if "name" in person:
        return person["name"]
    return " ".join(person[field] for field in ("given-names", "name-particle", "family-names", "name-suffix") if person.get(field))


def citation_authors(reference):
    authors = [person_name(author) for author in reference["authors"]]
    return authors[0] + " et al." if len(authors) > 3 else ", ".join(authors)


def arxiv_id(reference):
    """arXiv identifier of a preprint recorded with its arXiv DOI, else None."""
    doi = reference.get("doi", "")
    prefix = "10.48550/arxiv."
    return doi[len(prefix):] if doi.lower().startswith(prefix) else None


def citation(reference):
    text = citation_authors(reference) + ', "' + reference["title"] + '."'
    venue = reference.get("journal") or reference.get("collection-title")
    if not venue and arxiv_id(reference):
        venue = "arXiv:" + arxiv_id(reference)
    if venue:
        text += " " + venue
        if "volume" in reference:
            text += " " + str(reference["volume"])
        if "start" in reference:
            text += ", " + str(reference["start"])
            if "end" in reference:
                text += "–" + str(reference["end"])
    if "year" in reference:
        text += " (" + str(reference["year"]) + ")"
    if venue or "year" in reference:
        text += "."
    if "status" in reference:
        text += " " + reference["status"].replace("-", " ").capitalize() + "."
    if "notes" in reference:
        text += " " + reference["notes"]
    if "doi" in reference:
        text += " https://doi.org/" + reference["doi"]
    elif "url" in reference:
        text += " " + reference["url"]
    return text


def detailed_notice(policy, references, framework, component):
    entry = policy["components"][component]
    selected = select_references(policy, framework, component)
    keys = list(dict.fromkeys(key for role in ROLES for key in selected[role]))
    numbers = {key: index + 1 for index, key in enumerate(keys)}
    lines = ["Recommended citations for " + entry["name"] + ":",
             textwrap.fill(policy["policy"]["request"], width=88)]
    for role in ROLES:
        if selected[role]:
            labels = ", ".join(f"[{numbers[key]}]" for key in selected[role])
            lines.append("  " + role.capitalize() + ": " + labels)
    for key in keys:
        lines.append(textwrap.fill(citation(references[key]), width=88,
                                   initial_indent=f"  [{numbers[key]}] ", subsequent_indent="      "))
    if entry.get("note"):
        lines.append(textwrap.fill(entry["note"], width=88))
    lines.append(textwrap.fill(policy["policy"]["license_note"], width=88))
    return "\n".join(lines) + "\n\n"


def compact_citation(reference):
    """Short authors and publication coordinates, or unpublished title/status.

    Preprints keep their title and add the arXiv identifier.
    """
    venue = reference.get("journal") or reference.get("collection-title")
    text = citation_authors(reference) + ", " + (venue or reference["title"])
    if venue:
        if "volume" in reference:
            text += " " + str(reference["volume"])
        if "start" in reference:
            text += ", " + str(reference["start"])
    elif arxiv_id(reference):
        text += ", arXiv:" + arxiv_id(reference)
    if "year" in reference:
        text += " (" + str(reference["year"]) + ")"
    if "status" in reference:
        text += "; " + reference["status"].replace("-", " ")
    return text + "."


def notice(policy, references, framework, component):
    """Short startup box; complete bibliography is available through --citations."""
    selected = select_references(policy, framework, component)
    keys = list(dict.fromkeys(key for role in ROLES for key in selected[role]))
    numbers = {key: index + 1 for index, key in enumerate(keys)}
    entry = policy["components"][component]
    width = 76
    lines = textwrap.wrap("Recommended citations for " + entry["name"] + ":", width)
    labels = [role.capitalize() + ": " + ", ".join(f"[{numbers[key]}]" for key in selected[role])
              for role in ROLES if selected[role]]
    lines += textwrap.wrap("; ".join(labels), width)
    for key in keys:
        lines += textwrap.wrap(compact_citation(references[key]), width,
                               initial_indent=f"[{numbers[key]}] ", subsequent_indent="    ")
    if entry.get("note"):
        lines += textwrap.wrap(entry["note"], width)
    lines += textwrap.wrap(policy["policy"]["license_note"], width)
    lines.append("Full references: --citations")
    border = "*" * (width + 4)
    return "\n".join([border] + ["* " + line.ljust(width) + " *" for line in lines] + [border]) + "\n\n"


def cpp_data(policy, references, framework):
    lines = ["// Generated from CITATION.cff and CITATIONS.yaml. Do not edit.",
             "static const citation_entry citation_entries[] = {"]
    for component in sorted(policy["components"]):
        # JSON escaping is also valid for these C++ UTF-8 string literals.
        lines.append("  {" + json.dumps(component) + ",")
        for render in (notice, detailed_notice):
            value = render(policy, references, framework, component)
            pieces = value.splitlines(keepends=True)
            lines.extend("   " + json.dumps(line, ensure_ascii=False) for line in pieces[:-1])
            lines.append("   " + json.dumps(pieces[-1], ensure_ascii=False) + ",")
        lines.append("  },")
    lines.extend(["};", ""])
    return "\n".join(lines)


def native_outputs(root, policy, references, framework):
    """Checked-in documentation and C++ data for builds without Python."""
    generated = ".github/scripts/citations/generated/"
    outputs = {"CITATION.md": markdown(policy, references, framework),
               generated + "citations_data.inc": cpp_data(policy, references, framework)}
    inputs = ("CITATION.cff", "CITATIONS.yaml", ".github/scripts/generate_citations.py",
              ".github/scripts/citations/cff-1.2.0.schema.json", ".github/scripts/citations/policy.schema.json")
    hashes = {path: hashlib.sha256((root / path).read_bytes()).hexdigest() for path in inputs}
    hashes.update({path: hashlib.sha256(content.encode("utf-8")).hexdigest() for path, content in outputs.items()})
    lines = ["# Generated by .github/scripts/generate_citations.py --regenerate. Do not edit.",
             "set(_alps_citation_generated_format 1)",
             "set(_alps_citation_files " + " ".join(hashes) + ")",
             "set(_alps_citation_hashes " + " ".join(hashes.values()) + ")"]
    outputs[generated + "snapshots.cmake"] = "\n".join(lines) + "\n"
    return outputs


def markdown(policy, references, framework):
    def links(keys):
        return ", ".join(f"[{key}](#{key.lower()})" for key in keys) or "N/A"

    lines = ["<!-- Generated by .github/scripts/generate_citations.py; do not edit. -->",
             "# Citing ALPS", "",
             "[CITATIONS.yaml](CITATIONS.yaml) is the citation-policy authority; "
             "[CITATION.cff](CITATION.cff) is the bibliographic authority. "
             "CLI notices and this page are generated from these files.", "",
             policy["policy"]["request"], "", policy["policy"]["license_note"], ""]
    for role in ROLES:
        lines.extend(["- **" + role.capitalize() + ":** " + policy["policy"][role]])
    lines.extend(["", "## Framework paper", "", citation(references[framework]), "",
                  "## Application citation matrix", "",
                  "The framework paper applies to every component. A paper serving several roles "
                  "is listed once in CLI output. Components in the Uses column contribute their "
                  "references as well.", "",
                  "| Component | Algorithm | Implementation | Uses |",
                  "| --- | --- | --- | --- |"])
    for component, entry in policy["components"].items():
        if component == "framework":
            continue
        uses = ", ".join(entry.get("uses", [])) or "—"
        lines.append(f"| `{component}` | {links(entry['algorithm'])} | {links(entry['implementation'])} | {uses} |")
    for entry in policy["components"].values():
        if entry.get("note"):
            lines.extend(["", entry["note"]])
    lines.extend(["", "## References", ""])
    for key, reference in references.items():
        lines.extend([f"### {key}", "", citation(reference), ""])
    source = policy["source"]
    lines.extend(["## Policy provenance", "",
                  f"Initially transcribed from {source['document']} at revision "
                  f"`{source['revision']}`: {source['sections']}.", "",
                  "Maintenance instructions: [.github/scripts/citations/README.md](.github/scripts/citations/README.md).", ""])
    return "\n".join(lines)


def write_if_changed(path, content):
    if not path.exists() or path.read_bytes() != content.encode("utf-8"):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="\n") as output:
            output.write(content)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--cpp", type=Path, help="Write generated C++ data to this build-tree path")
    parser.add_argument("--markdown", type=Path, help="Write generated readable guidance")
    parser.add_argument("--regenerate", action="store_true", help="Refresh checked-in native data and CITATION.md")
    parser.add_argument("--check", action="store_true", help="Require all checked-in generated files to be current")
    args = parser.parse_args()
    try:
        _, policy, references, framework = load_catalog(args.root)
        doc = markdown(policy, references, framework)
        if args.check or args.regenerate:
            outputs = native_outputs(args.root, policy, references, framework)
            for filename, content in outputs.items():
                path = args.root / filename
                if args.check and (not path.exists() or path.read_bytes() != content.encode("utf-8")):
                    raise ValueError(filename + " is stale; run python .github/scripts/generate_citations.py --regenerate")
                if args.regenerate:
                    write_if_changed(path, content)
        if args.cpp:
            write_if_changed(args.cpp, cpp_data(policy, references, framework))
        if args.markdown:
            write_if_changed(args.markdown, doc)
    except (ValueError, OSError, yaml.YAMLError, jsonschema.ValidationError) as exc:
        print(f"Citation validation failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
