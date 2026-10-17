# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
#
# The parametric side of alpsx's Python Hamiltonian (skilledwolf/alpsx,
# private repository, branch agent/v3-execution-and-api-audit at d37e51227,
# python/alps/src/alps/models/hamiltonian.py: symbolic coefficients,
# Hamiltonian.bind and free_symbols), written anew over the bound
# Hamiltonian and TermSum, whose coefficients are numbers.
"""Hamiltonians with symbolic coefficients, bound to numbers on demand.

>>> from pyalps import models
>>> J, h = models.symbols("J", "h")
>>> spin = models.spin()
>>> model = models.ParametricHamiltonian(models.lattice.chain(8), spin, conserved_quantum_numbers=["Sz"])
>>> _ = model.add_bond(J, lambda i, j: spin.Sz(i) * spin.Sz(j) + 0.5 * (spin.Sp(i) * spin.Sm(j)).add_hc())
>>> _ = model.add_site(-h, lambda i: spin.Sz(i))
>>> sorted(s.name for s in model.free_symbols())
['J', 'h']
>>> numeric = model.bind(J=1.0, h=0.1)              # a models.Hamiltonian

``bind`` with some symbols left returns another ``ParametricHamiltonian``;
the solvers take the bound ``Hamiltonian`` only.
"""

from __future__ import annotations

from .symbols import Expr, Symbol, _is_scalar, _symbol_values


def _terms_type():
    from . import TermSum

    return TermSum


class SymbolicTermSum:
    """Located terms with symbolic coefficients: a sum of ``coefficient * TermSum`` parts.

    It arises from ``Symbol * TermSum`` (either order) and keeps the operator
    algebra: sums with TermSum and SymbolicTermSum values, products with
    numbers, symbols and operators, ``dag``/``add_hc``. ``bind`` returns a
    TermSum when every symbol is bound.
    """

    __slots__ = ("parts",)

    def __init__(self, parts=()):
        self.parts = [(Expr.of(c), t) for c, t in parts]

    @classmethod
    def of(cls, value) -> "SymbolicTermSum":
        if isinstance(value, SymbolicTermSum):
            return value
        if isinstance(value, _terms_type()):
            return cls([(1.0, value)])
        raise TypeError(f"not an operator: {type(value).__name__}")

    def free_symbols(self) -> set[Symbol]:
        return {s for c, _ in self.parts for s in c.free_symbols()}

    def bind(self, values=None, **kwargs):
        """The parts with ``values`` substituted: a TermSum when nothing stays symbolic."""
        bound = [(c.bind(values, **kwargs), t) for c, t in self.parts]
        if any(isinstance(c, Expr) for c, _ in bound):
            return SymbolicTermSum(bound)
        total = _terms_type()()
        for c, t in bound:
            total += c * t
        return total

    def dag(self) -> "SymbolicTermSum":
        return SymbolicTermSum([(c.conjugate(), t.dag()) for c, t in self.parts])

    def add_hc(self) -> "SymbolicTermSum":
        return self + self.dag()

    # -- algebra --------------------------------------------------------------

    def _coerce(self, other):
        if isinstance(other, (SymbolicTermSum, _terms_type())):
            return SymbolicTermSum.of(other)
        return None

    def __add__(self, other):
        other = self._coerce(other)
        return NotImplemented if other is None else SymbolicTermSum(self.parts + other.parts)

    def __radd__(self, other):
        other = self._coerce(other)
        return NotImplemented if other is None else SymbolicTermSum(other.parts + self.parts)

    def __sub__(self, other):
        other = self._coerce(other)
        return NotImplemented if other is None else self + (-other)

    def __rsub__(self, other):
        other = self._coerce(other)
        return NotImplemented if other is None else other + (-self)

    def __neg__(self):
        return SymbolicTermSum([(-c, t) for c, t in self.parts])

    def __mul__(self, other):
        if _is_scalar(other) or isinstance(other, (Symbol, Expr)):
            factor = Expr.of(other)
            return SymbolicTermSum([(c * factor, t) for c, t in self.parts])
        operator = self._coerce(other)
        if operator is None:
            return NotImplemented
        return SymbolicTermSum([(ca * cb, ta * tb) for ca, ta in self.parts for cb, tb in operator.parts])

    def __rmul__(self, other):
        if _is_scalar(other) or isinstance(other, (Symbol, Expr)):
            return self.__mul__(other)
        operator = self._coerce(other)
        if operator is None:
            return NotImplemented
        return operator.__mul__(self)

    def __truediv__(self, other):
        if not _is_scalar(other):
            return NotImplemented
        return self * (1.0 / other)

    def __len__(self):
        return len(self.parts)

    def __repr__(self):
        return "SymbolicTermSum(" + " + ".join(f"({c})*{t!r}" for c, t in self.parts) + ")"


class ParametricHamiltonian:
    """A ``Hamiltonian`` whose coefficients may be symbols, bound to numbers by ``bind``.

    The constructor takes the arguments of ``models.Hamiltonian``;
    ``add_bond``, ``add_site`` and ``add`` take those of the Hamiltonian
    methods, with a number, Symbol or Expr as coefficient and specifications
    that may return SymbolicTermSum values. Specifications are called when
    added, as for ``Hamiltonian``.
    """

    def __init__(self, lattice, sites, fermion_convention="strict", conserved_quantum_numbers=()):
        from . import Hamiltonian

        self._arguments = (lattice, sites, fermion_convention, list(conserved_quantum_numbers))
        self._skeleton = Hamiltonian(*self._arguments)  # validates the arguments; never solved
        self._pieces = []  # (SymbolicTermSum, add_hc), in the order added

    @property
    def lattice(self):
        return self._skeleton.lattice

    @property
    def num_sites(self):
        return self._skeleton.num_sites

    @property
    def site_dims(self):
        return self._skeleton.site_dims

    @property
    def conserved_quantum_numbers(self):
        return self._skeleton.conserved_quantum_numbers

    def site_for(self, site):
        return self._skeleton.site_for(site)

    def add_bond(self, coeff, spec, bond_type=None, add_hc=False):
        for bond in self.lattice.bonds:
            if bond_type is None or bond.bond_type == bond_type:
                self._pieces.append((coeff * SymbolicTermSum.of(spec(bond.source, bond.target)), add_hc))
        return self

    def add_site(self, coeff, spec, site_type=None, add_hc=False):
        for site in self.lattice.sites:
            if site_type is None or site.site_type == site_type:
                self._pieces.append((coeff * SymbolicTermSum.of(spec(site.index)), add_hc))
        return self

    def add(self, term, add_hc=False):
        self._pieces.append((SymbolicTermSum.of(term), add_hc))
        return self

    def free_symbols(self) -> set[Symbol]:
        return {s for piece, _ in self._pieces for s in piece.free_symbols()}

    def bind(self, values=None, **kwargs):
        """A ``Hamiltonian`` with ``values`` substituted, or a ``ParametricHamiltonian`` if symbols remain."""
        _symbol_values(values, kwargs)  # rejects non-numeric values up front
        bound = [(piece.bind(values, **kwargs), add_hc) for piece, add_hc in self._pieces]
        if any(isinstance(piece, SymbolicTermSum) for piece, _ in bound):
            partial = ParametricHamiltonian(*self._arguments)
            partial._pieces = [(SymbolicTermSum.of(piece), add_hc) for piece, add_hc in bound]
            return partial
        from . import Hamiltonian

        model = Hamiltonian(*self._arguments)
        for piece, add_hc in bound:
            model.add(piece, add_hc)
        return model

    def __repr__(self):
        free = ", ".join(sorted(s.name for s in self.free_symbols())) or "none"
        return f"ParametricHamiltonian({self.num_sites} sites, {len(self._pieces)} pieces, free symbols: {free})"


def require_bound(model):
    """Raise a clear error for a ParametricHamiltonian handed to a solver."""
    if isinstance(model, ParametricHamiltonian):
        free = ", ".join(sorted(s.name for s in model.free_symbols())) or "none"
        raise TypeError(f"a ParametricHamiltonian (free symbols: {free}) needs bind(...) before it is solved")


__all__ = ["SymbolicTermSum", "ParametricHamiltonian"]
