# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/alps/src/alps/models/validate.py. Changes: the targets are the
# solvers pyalps binds (None, "ed", "dmrg"); alpsx's QMC targets (sse,
# looper, worm, qwl, classical_mc) and their checks are dropped, as no
# pyalps QMC solver takes a kit Hamiltonian; the checks call the bound
# compile_hamiltonian and validate_hermiticity, the ed and dmrg targets add a
# Hermiticity check, unbound ParametricHamiltonian symbols are an error of
# their own, and the suggestions name pyalps calls.
"""Target-aware checks of a ``Hamiltonian`` before it reaches a solver.

>>> from pyalps import models
>>> report = models.validate(models.heisenberg(models.lattice.chain(64)), "ed")
>>> report.ok, [d.code for d in report]
(True, ['ed.large_space'])
>>> report.raise_for_errors()                    # raises ValueError on errors only

``validate`` takes a ``Hamiltonian`` or ``ParametricHamiltonian`` (a
``CatalogModel`` is checked by the ED and DMRG entry points as they compile
it). Diagnostics carry a level (error, warning, info), a stable code, a
message and an optional suggestion.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Diagnostic:
    """One finding: ``level`` (error, warning, info), a stable ``code``, ``message`` and ``suggestion``."""

    level: str
    code: str
    message: str
    suggestion: str | None = None

    def __str__(self) -> str:
        out = f"[{self.level.upper()}] {self.code}: {self.message}"
        if self.suggestion:
            out += f"\n  -> {self.suggestion}"
        return out


@dataclass
class ValidationReport:
    """The diagnostics of :func:`validate`; true when no error was found."""

    target: str | None
    diagnostics: list[Diagnostic] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    @property
    def errors(self) -> list[Diagnostic]:
        return [d for d in self.diagnostics if d.level == "error"]

    @property
    def warnings(self) -> list[Diagnostic]:
        return [d for d in self.diagnostics if d.level == "warning"]

    def raise_for_errors(self) -> None:
        if self.ok:
            return
        body = "\n".join(str(d) for d in self.errors)
        raise ValueError(f"validate(target={self.target!r}) reported {len(self.errors)} error(s):\n{body}")

    def _add(self, *args) -> None:
        self.diagnostics.append(Diagnostic(*args))

    def __bool__(self) -> bool:
        return self.ok

    def __iter__(self):
        return iter(self.diagnostics)

    def __len__(self) -> int:
        return len(self.diagnostics)

    def __repr__(self) -> str:
        return f"ValidationReport(target={self.target!r}, errors={len(self.errors)}, warnings={len(self.warnings)})"


KNOWN_TARGETS = {None, "ed", "dmrg"}

# The unrestricted product dimension above which ED deserves a sector (alpsx's threshold).
LARGE_SPACE = 2 ** 14


def validate(model, target: str | None = None) -> ValidationReport:
    """Structural checks of ``model`` and, for ``target`` "ed" or "dmrg", that solver's."""
    from . import Hamiltonian
    from .parametric import ParametricHamiltonian

    if not isinstance(model, (Hamiltonian, ParametricHamiltonian)):
        raise TypeError(f"validate takes a Hamiltonian or ParametricHamiltonian, not {type(model).__name__}")
    report = ValidationReport(target=target)
    if target not in KNOWN_TARGETS:
        report._add("warning", "unknown_target", f"target {target!r} is not a solver pyalps binds",
                    "known targets: 'ed', 'dmrg' or None; running the structural checks only")
    if isinstance(model, ParametricHamiltonian):
        free = sorted(s.name for s in model.free_symbols())
        if free:
            report._add("error", "unbound_parameters", f"the coefficients depend on unbound symbols {free}",
                        "call model.bind(...) before validating or solving")
            return report
        model = model.bind()
    if not _generic_checks(model, report):
        return report
    if target == "ed":
        _hermiticity_check(model, report)
        _ed_checks(model, report)
    elif target == "dmrg":
        _hermiticity_check(model, report)
        _dmrg_checks(model, report)
    return report


def _generic_checks(model, report: ValidationReport) -> bool:
    """Terms and the declared conservation laws; false when nothing further can run."""
    from . import compile_hamiltonian

    if model.expand().empty:
        report._add("warning", "empty", "the Hamiltonian has no terms")
    try:
        compile_hamiltonian(model)
    except ValueError as exc:
        if "conserve" in str(exc):
            report._add("error", "conservation_broken", str(exc),
                        "drop the charge from conserved_quantum_numbers or remove the term that breaks it")
        else:
            report._add("error", "compile_failed", str(exc))
        return False
    return True


def _hermiticity_check(model, report: ValidationReport) -> None:
    from . import validate_hermiticity

    try:
        validate_hermiticity(model)
    except ValueError as exc:
        if "must be Hermitian" in str(exc):
            report._add("error", "not_hermitian", str(exc),
                        "add the Hermitian conjugates (add_hc=True, or TermSum.add_hc())")
        else:
            report._add("warning", "hermiticity_unchecked", str(exc))


def _ed_checks(model, report: ValidationReport) -> None:
    # The unrestricted product dimension, not the size of any sector solved; in
    # logarithms, as a long chain's dimension has more digits than str() prints.
    digits = sum(math.log10(d) for d in model.site_dims)
    if digits >= math.log10(LARGE_SPACE):
        size = str(math.prod(model.site_dims)) if digits < 15 else f"about 10^{digits:.0f}"
        report._add("warning", "ed.large_space",
                    f"the unrestricted product basis has {size} states; the cost depends on the sector",
                    "fix charges with ed.prepare_product(model, sector={...}) and ask eigenstates(count=...) "
                    "for the low-lying levels only")


def _dmrg_checks(model, report: ValidationReport) -> None:
    lattice = model.lattice
    if lattice.dimension >= 2:
        report._add("warning", "dmrg.dimension",
                    f"the lattice is {lattice.dimension}-dimensional; DMRG cost grows exponentially with the width",
                    "use a chain or a narrow ladder, or ed for small clusters")
    elif lattice.periodic_directions():
        report._add("warning", "dmrg.periodic", "the chain is periodic; DMRG converges far faster on open chains",
                    "build the lattice with bc='open' if the boundary allows")
    if model.conserved_quantum_numbers:
        report._add("info", "dmrg.sector_required",
                    f"the model conserves {list(model.conserved_quantum_numbers)}; DMRG solves one sector at a time",
                    "pass sector={...} to classic_dmrg.solve")


__all__ = ["Diagnostic", "ValidationReport", "validate", "KNOWN_TARGETS"]
