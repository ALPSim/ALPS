# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
# Adapted from skilledwolf/alpsx, examples/demo_support/__init__.py, source worktree
# 98ed/ALPS at 1adba871e (including its local edits), SHA-256 c6ef62796dade877a42425bb64f231e45d64d31c52c9ed2ebb1f58e2d3dc5cdb.
# Port: pyalps model/classic-DMRG APIs; final-state RDMs and capability-gated QDK export.
"""Small independent references and artifact helpers for the demonstrations.

Dense/occupation-bit references here are deliberately limited to tiny examples;
they are checks of the production solvers, not alternative solver backends.
"""

import argparse
import datetime
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

PAULI = {"I": np.eye(2, dtype=complex), "X": np.array([[0, 1], [1, 0]], dtype=complex),
         "Y": np.array([[0, -1j], [1j, 0]]), "Z": np.diag([1, -1]).astype(complex)}


def parser(description):
    result = argparse.ArgumentParser(description=description)
    result.add_argument("--output", type=Path, help="save validated data, provenance and a PNG figure")
    return result


def site_matrix(n, factors):
    """Independent Kronecker reference, site 0 most significant."""
    result = np.ones((1, 1), dtype=complex)
    for site in range(n):
        result = np.kron(result, PAULI[factors.get(site, "I")])
    return result


def pauli_matrix(operator):
    """Read the QDK container directly, without the ALPS spin importer."""
    data = operator.to_dict()
    coefficients = np.asarray(data["coefficients"]["real"]) + 1j * np.asarray(data["coefficients"]["imag"])
    result = np.zeros((2**operator.num_qubits,) * 2, dtype=complex)
    for word, coefficient in zip(data["pauli_strings"], coefficients, strict=True):
        term = np.ones((1, 1), dtype=complex)
        for p in word:
            term = np.kron(term, PAULI[p])
        result += coefficient * term
    return result


def bit_reversal(n):
    return np.array([int(f"{i:0{n}b}"[::-1], 2) for i in range(2**n)])


def fock_matrix(states, products):
    """Apply ordered CAR monomials to explicitly listed occupation integers.

    Mode 0 is the least significant occupation bit. Each product is
    (coefficient, [(mode, creation), ...]) in left-to-right operator order.
    Only these small reference examples enumerate the Hilbert space.
    """
    lookup = {state: i for i, state in enumerate(states)}
    matrix = np.zeros((len(states), len(states)), dtype=complex)
    for coefficient, factors in products:
        for column, original in enumerate(states):
            state, value = original, complex(coefficient)
            for mode, creation in reversed(factors):
                occupied = bool(state & (1 << mode))
                if occupied == creation:
                    value = 0.
                    break
                value *= (-1) ** (state & ((1 << mode) - 1)).bit_count()
                state ^= 1 << mode
            if value and state in lookup:
                matrix[lookup[state], column] += value
    return matrix


def finish(args, data, started, plot):
    elapsed = time.perf_counter() - started
    print(f"Validated in {elapsed:.2f} s")
    if args.output is None:
        return
    directory = args.output.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    script = Path(sys.argv[0]).resolve()
    repository = Path(__file__).resolve().parents[4]

    def git(*arguments):
        result = subprocess.run(["git", "-C", str(repository), *arguments],
                                capture_output=True, text=True, timeout=5, check=False)
        return result.stdout.strip() if result.returncode == 0 else "unavailable"

    versions = {}
    for package in ("pyalps", "qdk", "qdk-chemistry", "qdk-chemistry-alps", "numpy", "scipy"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "source checkout"
    data["provenance"] = {"revision": git("rev-parse", "HEAD"),
                          "generated_utc": datetime.datetime.now(datetime.UTC).isoformat(),
                          "platform": platform.platform(),
                          "working_tree_dirty": bool(git("status", "--porcelain")),
                          "tracked_diff_sha256": hashlib.sha256(git("diff", "HEAD").encode()).hexdigest(),
                          "script": script.name,
                          "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
                          "versions": versions, "python": sys.version,
                          "runtime_seconds": elapsed, "saved_results": True}

    def convert(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)

    (directory / f"{script.stem}.json").write_text(json.dumps(data, default=convert, indent=2) + "\n")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "alps-demo-mpl"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with plt.rc_context({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}):
        figure = plot(plt, data)
        figure.savefig(directory / f"{script.stem}.png", dpi=170, bbox_inches="tight")
        plt.close(figure)
    print(f"Saved results and figure: {directory / script.stem}")


def model_matrix(model):
    """Small ED reference in site-product order; never used by DMRG."""
    from pyalps import ed, models
    sites = [model.site_for(i) for i in range(model.num_sites)]
    problem = ed.prepare_product(models.located(sites[0], "I", 0), sites, [], model.fermion_convention)
    indices = np.zeros(len(problem.states), dtype=int)
    for site, dim in enumerate(problem.site_dims):
        indices = indices * int(dim) + np.asarray(problem.states)[:, site]
    matrix = problem.operator_matrix(model.operator_decomposition().op).to_scipy().toarray()
    result = np.empty_like(matrix)
    result[np.ix_(indices, indices)] = matrix
    return result


def export_operator(model, observable=None, encoding=None, tolerance=0.):
    """Use the native Pauli compiler for a Hamiltonian or an observable."""
    from pyalps import models
    if observable is None:
        return models.export_pauli(model, encoding=encoding, tolerance=tolerance)
    basis = models.physical_basis(model)
    op = models.compile_operator(observable, basis)
    return models.export_pauli(op, basis, encoding=encoding, tolerance=tolerance,
                              require_hermitian=False, require_even=False)
