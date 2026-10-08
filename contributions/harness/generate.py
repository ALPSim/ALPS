#!/usr/bin/env python3
"""Generate the harness glue for a submission from its manifest.yml.

Usage: generate.py <submission_dir> <out_dir>

Writes <out_dir>/algorithm.hpp (the Solver wrapper for the manifest's
contract, or a copy of the submission's own for `contract: custom`),
<out_dir>/selection.hpp (which catalog cases the manifest declares) and
<out_dir>/name. Every manifest value that reaches C++ is checked against a
fixed list or pattern first, so nothing from a submission is pasted in raw.

Standard library only: the manifest is a flat subset of YAML, `key: value`
and `key: [a, b]` (or a `- item` block list).
"""
import re
import shutil
import sys
from pathlib import Path

# --- What the catalog offers; keep in step with problems.cpp ----------------

MODELS = ["tight-binding"]
GEOMETRIES = ["1d", "2d"]
LATTICES = ["open chain lattice", "chain lattice", "dimer",
            "open square lattice", "square lattice"]
QUANTITIES = ["energy", "correlation"]
# Model parameters a `model_parameters` function can take, with their C++ type.
PARAMETERS = {"L": "int", "N": "int", "t": "double", "V": "double"}

# --- Contracts ----------------------------------------------------------------
# Each contract fixes the contributor's function signature and what the
# generated wrapper can derive from it.

CONTRACTS = {
    "dense_eigenvalues": {
        "quantities": ["energy"],
        "signature": "std::vector<double> {entry}(std::vector<double> h, int n)",
    },
    "dense_eigensystem": {
        "quantities": ["energy", "correlation"],
        "signature": "void {entry}(std::vector<double> h, int n, "
                     "std::vector<double>& values, std::vector<double>& vectors)",
    },
    "model_parameters": {
        "quantities": ["energy"],
        "signature": "double {entry}({inputs})",
    },
    "custom": {
        "quantities": QUANTITIES,
        "signature": None,
    },
}

NAME_RE = re.compile(r"^[A-Za-z0-9_-]+$")
IDENT_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
RESERVED = {"main", "onboard", "std"}

KEYS = {"name", "description", "entry", "contract", "models", "geometry",
        "lattices", "quantities", "inputs", "representation"}


class ManifestError(Exception):
    pass


def parse_manifest(path):
    data, last_key = {}, None
    for lineno, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        item = re.match(r"^\s+-\s*(.+)$", line)
        if item:
            if last_key is None or not isinstance(data.get(last_key), list):
                raise ManifestError(f"line {lineno}: list item without a list key")
            data[last_key].append(unquote(item.group(1)))
            continue
        kv = re.match(r"^([A-Za-z_]+):\s*(.*)$", line)
        if not kv:
            raise ManifestError(f"line {lineno}: expected 'key: value', got {raw!r}")
        key, value = kv.group(1), kv.group(2).strip()
        if key not in KEYS:
            raise ManifestError(f"line {lineno}: unknown key '{key}' (known: {', '.join(sorted(KEYS))})")
        if key in data:
            raise ManifestError(f"line {lineno}: '{key}' given twice")
        if value == "":
            data[key] = []
        elif value.startswith("["):
            if not value.endswith("]"):
                raise ManifestError(f"line {lineno}: unterminated list")
            inner = value[1:-1].strip()
            data[key] = [unquote(v) for v in inner.split(",")] if inner else []
        else:
            data[key] = unquote(value)
        last_key = key
    return data


def unquote(s):
    s = s.strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in "\"'":
        return s[1:-1]
    return s


def as_list(data, key):
    v = data.get(key)
    if v is None:
        return None
    return v if isinstance(v, list) else [v]


def check_subset(key, values, allowed):
    bad = [v for v in values if v not in allowed]
    if bad:
        raise ManifestError(f"{key}: unknown {bad}; choose from {allowed}")


def validate(data, submission):
    name = data.get("name")
    if not isinstance(name, str) or not NAME_RE.match(name):
        raise ManifestError("name: required, letters, digits, '_' or '-'")

    contract = data.get("contract")
    if contract not in CONTRACTS:
        raise ManifestError(f"contract: required, one of {list(CONTRACTS)}")

    models = as_list(data, "models")
    if not models:
        raise ManifestError(f"models: required, from {MODELS}")
    check_subset("models", models, MODELS)

    geometry = as_list(data, "geometry") or GEOMETRIES
    check_subset("geometry", geometry, GEOMETRIES)

    lattices = as_list(data, "lattices")
    if lattices is not None:
        check_subset("lattices", lattices, LATTICES)

    provided = CONTRACTS[contract]["quantities"]
    quantities = as_list(data, "quantities")
    if quantities is None:
        quantities = provided if contract != "custom" else ["energy"]
    check_subset("quantities", quantities, provided)

    entry, inputs = data.get("entry"), as_list(data, "inputs")
    has_own_wrapper = (submission / "algorithm.hpp").exists()
    if contract == "custom":
        if not has_own_wrapper:
            raise ManifestError("contract: custom needs your own algorithm.hpp")
    else:
        if has_own_wrapper:
            raise ManifestError(f"contract: {contract} generates algorithm.hpp; "
                                "remove yours, or use contract: custom")
        if not isinstance(entry, str) or not IDENT_RE.match(entry) or entry in RESERVED:
            raise ManifestError("entry: required, the name of your C++ function")

    if contract == "model_parameters":
        if not inputs:
            raise ManifestError(f"inputs: required for model_parameters, from {list(PARAMETERS)}")
        check_subset("inputs", inputs, list(PARAMETERS))
        if len(set(inputs)) != len(inputs):
            raise ManifestError("inputs: a parameter is listed twice")
    elif inputs:
        raise ManifestError("inputs: only used with contract: model_parameters")

    return {"name": name, "contract": contract, "entry": entry, "models": models,
            "geometry": geometry, "lattices": lattices, "quantities": quantities,
            "inputs": inputs or []}


# --- Code generation ----------------------------------------------------------

HEADER = """\
#pragma once
// GENERATED by contributions/harness/generate.py from manifest.yml. Do not edit.
"""

WRAPPERS = {
    "dense_eigenvalues": """\
#include "problem.hpp"
#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <vector>

// Contract dense_eigenvalues: all eigenvalues of the symmetric n x n row-major
// matrix h, in any order.
{signature};

namespace onboard {{

class Algorithm final : public Solver {{
public:
    explicit Algorithm(const MatrixProblem& p) : p_(p) {{}}

    std::optional<Estimate> groundStateEnergy() const override {{
        const int n = p_.numSites();
        std::vector<double> levels = {entry}(p_.singleParticleMatrix(), n);
        if (static_cast<int>(levels.size()) != n) return std::nullopt;
        std::sort(levels.begin(), levels.end());

        double e = 0.0;
        for (int k = 0; k < p_.numParticles(); ++k) e += levels[k];
        if (!std::isfinite(e)) return std::nullopt;
        return Estimate{{e}};
    }}

    const char* name() const override {{ return "{name}"; }}

private:
    const MatrixProblem& p_;
}};

inline std::unique_ptr<const Solver> makeAlgorithm(const MatrixProblem& p) {{
    return std::make_unique<Algorithm>(p);
}}

}} // namespace onboard
""",
    "dense_eigensystem": """\
#include "problem.hpp"
#include "solver.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <optional>
#include <vector>

// Contract dense_eigensystem: for the symmetric n x n row-major matrix h, fill
// `values` with all n eigenvalues and `vectors` with the matching normalized
// eigenvectors, eigenvector k at vectors[k*n .. k*n + n).
{signature};

namespace onboard {{

class Algorithm final : public Solver {{
public:
    explicit Algorithm(const MatrixProblem& p) : p_(p) {{
        const int n = p_.numSites();
        std::vector<double> values, vectors;
        {entry}(p_.singleParticleMatrix(), n, values, vectors);
        if (static_cast<int>(values.size()) != n || static_cast<int>(vectors.size()) != n * n)
            return;

        // The numParticles() lowest states, in order.
        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {{ return values[a] < values[b]; }});

        double e = 0.0;
        std::vector<double> g(static_cast<std::size_t>(n) * n, 0.0);
        for (int k = 0; k < p_.numParticles(); ++k) {{
            const int s = order[k];
            e += values[s];
            const double* v = &vectors[static_cast<std::size_t>(s) * n];
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j) g[i * n + j] += v[i] * v[j];
        }}
        if (std::isfinite(e)) energy_ = Estimate{{e}};
        correlation_ = MatrixEstimate{{std::move(g)}};
    }}

    bool provides(Quantity) const override {{ return true; }}
    std::optional<Estimate> groundStateEnergy() const override {{ return energy_; }}
    std::optional<MatrixEstimate> correlation() const override {{ return correlation_; }}
    const char* name() const override {{ return "{name}"; }}

private:
    const MatrixProblem& p_;
    std::optional<Estimate> energy_;
    std::optional<MatrixEstimate> correlation_;
}};

inline std::unique_ptr<const Solver> makeAlgorithm(const MatrixProblem& p) {{
    return std::make_unique<Algorithm>(p);
}}

}} // namespace onboard
""",
    "model_parameters": """\
#include "problem.hpp"
#include "solver.hpp"

#include <cmath>
#include <memory>
#include <optional>

// Contract model_parameters: the ground-state energy from named model
// parameters, passed in the order the manifest's `inputs` lists them.
{signature};

namespace onboard {{

class Algorithm final : public Solver {{
public:
    explicit Algorithm(const Problem& p) : p_(p) {{}}

    std::optional<Estimate> groundStateEnergy() const override {{
{fetch}
        const double e = {entry}({args});
        if (!std::isfinite(e)) return std::nullopt;
        return Estimate{{e}};
    }}

    const char* name() const override {{ return "{name}"; }}

private:
    const Problem& p_;
}};

inline std::unique_ptr<const Solver> makeAlgorithm(const Problem& p) {{
    return std::make_unique<Algorithm>(p);
}}

}} // namespace onboard
""",
}


def cpp_list(values):
    return "{" + ", ".join(f'"{v}"' for v in values) + "}"


def selection_hpp(m):
    lines = [HEADER, '#include "problem.hpp"', "", "#include <cstring>",
             "#include <initializer_list>", "", "namespace onboard {", "",
             "// nullptr if the manifest declares this case, otherwise why not.",
             "inline const char* notSelected(const Problem& p, Quantity q) {",
             "    auto in = [](const char* s, std::initializer_list<const char*> xs) {",
             "        for (const char* x : xs) if (std::strcmp(s, x) == 0) return true;",
             "        return false;",
             "    };",
             f"    if (!in(p.model(), {cpp_list(m['models'])})) return \"model not declared\";",
             f"    if (!in(p.geometry(), {cpp_list(m['geometry'])})) return \"geometry not declared\";"]
    if m["lattices"] is not None:
        lines.append(f"    if (!in(p.lattice(), {cpp_list(m['lattices'])})) return \"lattice not declared\";")
    for key in m["inputs"]:
        lines.append(f'    if (!p.parameter("{key}")) return "no parameter {key}";')
    lines.append(f"    if (!in(quantityName(q), {cpp_list(m['quantities'])})) return \"quantity not declared\";")
    lines += ["    return nullptr;", "}", "", "} // namespace onboard", ""]
    return "\n".join(lines)


def algorithm_hpp(m):
    contract, entry, inputs = m["contract"], m["entry"], m["inputs"]
    params = ", ".join(f"{PARAMETERS[k]} {k}" for k in inputs)
    signature = CONTRACTS[contract]["signature"].format(entry=entry, inputs=params)
    fetch = "\n".join(
        f'        const auto {k}_ = p_.parameter("{k}");\n'
        f"        if (!{k}_) return std::nullopt;" for k in inputs)
    args = ", ".join(f"static_cast<{PARAMETERS[k]}>(*{k}_)" for k in inputs)
    body = WRAPPERS[contract].format(signature=signature, entry=entry, name=m["name"],
                                     fetch=fetch, args=args)
    return HEADER + "\n" + body


def main():
    if len(sys.argv) != 3:
        sys.exit("usage: generate.py <submission_dir> <out_dir>")
    submission, out = Path(sys.argv[1]), Path(sys.argv[2])
    manifest = submission / "manifest.yml"
    try:
        if not manifest.exists():
            raise ManifestError("submission is missing manifest.yml")
        m = validate(parse_manifest(manifest), submission)
    except ManifestError as e:
        print(f"::error::manifest.yml: {e}")
        sys.exit(1)

    out.mkdir(parents=True, exist_ok=True)
    if m["contract"] == "custom":
        shutil.copyfile(submission / "algorithm.hpp", out / "algorithm.hpp")
    else:
        (out / "algorithm.hpp").write_text(algorithm_hpp(m))
    (out / "selection.hpp").write_text(selection_hpp(m))
    (out / "name").write_text(m["name"] + "\n")
    signature = CONTRACTS[m["contract"]]["signature"]
    if signature:
        params = ", ".join(f"{PARAMETERS[k]} {k}" for k in m["inputs"])
        (out / "signature").write_text(signature.format(entry=m["entry"], inputs=params) + "\n")

    print(f"Submission '{m['name']}': contract {m['contract']}, models {m['models']}, "
          f"geometry {m['geometry']}, quantities {m['quantities']}"
          + (f", lattices {m['lattices']}" if m["lattices"] is not None else "")
          + (f", inputs {m['inputs']}" if m["inputs"] else ""))


if __name__ == "__main__":
    main()
