# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/alps/src/alps/models/symbols.py. Changes: numbers.Number scalars
# (NumPy scalars included) instead of (int, float, complex); a Symbol or Expr
# times a pyalps TermSum is a SymbolicTermSum (parametric.py), as alpsx's
# coefficients were stored in its Python Hamiltonian; no __float__,
# __complex__ or __index__, so an unbound parameter never converts silently.
"""Symbolic polynomial parameters for Hamiltonian coefficients.

>>> from pyalps import models
>>> J, h = models.symbols("J", "h")
>>> spin = models.spin()
>>> bond = J * spin.Sz(0) * spin.Sz(1) - h * spin.Sz(0)    # a SymbolicTermSum
>>> sorted(s.name for s in bond.free_symbols())
['J', 'h']
>>> numeric = bond.bind(J=1.0, h=0.5)                        # a TermSum

``Symbol`` and ``Expr`` form a commutative algebra over the complex numbers:
sums, products (polynomials with positive integer exponents), negation,
division by a number and complex conjugation. ``Expr.bind`` substitutes
values and returns a number once every symbol is bound, else an ``Expr``.
"""

from __future__ import annotations

import numbers
from dataclasses import dataclass, field
from typing import Iterable

Monomial = frozenset  # frozenset[(Symbol, exponent)]


def _is_scalar(value) -> bool:
    return isinstance(value, numbers.Number) and not isinstance(value, bool)


def _is_operator(value) -> bool:
    from . import TermSum
    from .parametric import SymbolicTermSum

    return isinstance(value, (TermSum, SymbolicTermSum))


def _scale(operator, coefficient: "Expr"):
    from .parametric import SymbolicTermSum

    return SymbolicTermSum.of(operator) * coefficient


@dataclass(frozen=True)
class Symbol:
    """A named scalar parameter."""

    name: str

    def __post_init__(self):
        if not self.name:
            raise ValueError("Symbol must have a non-empty name")

    def to_expr(self) -> "Expr":
        return Expr.from_symbol(self)

    def __add__(self, other):
        return self.to_expr() + other

    def __radd__(self, other):
        return self.to_expr().__radd__(other)

    def __sub__(self, other):
        return self.to_expr() - other

    def __rsub__(self, other):
        return self.to_expr().__rsub__(other)

    def __mul__(self, other):
        return self.to_expr() * other

    def __rmul__(self, other):
        return self.to_expr().__rmul__(other)

    def __truediv__(self, other):
        return self.to_expr() / other

    def __neg__(self):
        return -self.to_expr()

    def __pow__(self, n):
        return self.to_expr() ** n

    def conjugate(self):
        return self.to_expr().conjugate()

    def free_symbols(self) -> set["Symbol"]:
        return {self}

    def bind(self, values=None, **kwargs):
        return self.to_expr().bind(values, **kwargs)

    def __repr__(self) -> str:
        return self.name


def _normalise_monomial(items: Iterable[tuple[Symbol, int]]) -> Monomial:
    """Combine repeated symbols and drop zero exponents."""
    acc: dict[Symbol, int] = {}
    for sym, exp in items:
        if exp == 0:
            continue
        acc[sym] = acc.get(sym, 0) + exp
    return frozenset((s, e) for s, e in acc.items() if e != 0)


def _mul_monomials(a: Monomial, b: Monomial) -> Monomial:
    return _normalise_monomial(list(a) + list(b))


@dataclass(frozen=True)
class Expr:
    """A polynomial in :class:`Symbol` parameters with complex coefficients."""

    monomials: dict = field(default_factory=dict)  # Monomial -> number

    @classmethod
    def from_constant(cls, c) -> "Expr":
        if c == 0:
            return cls({})
        return cls({frozenset(): _as_scalar(c)})

    @classmethod
    def from_symbol(cls, s: Symbol) -> "Expr":
        return cls({frozenset({(s, 1)}): 1.0})

    @classmethod
    def of(cls, value) -> "Expr":
        """A number, Symbol or Expr as an Expr."""
        if isinstance(value, Expr):
            return value
        if isinstance(value, Symbol):
            return value.to_expr()
        return cls.from_constant(value)

    def is_zero(self) -> bool:
        return not self.monomials or all(c == 0 for c in self.monomials.values())

    def is_constant(self) -> bool:
        return all(m == frozenset() for m in self.monomials)

    def constant_value(self):
        return self.monomials.get(frozenset(), 0)

    def free_symbols(self) -> set[Symbol]:
        return {s for mon in self.monomials for s, _ in mon}

    # -- arithmetic ----------------------------------------------------------

    def __neg__(self) -> "Expr":
        return Expr({m: -c for m, c in self.monomials.items()})

    def __add__(self, other):
        if isinstance(other, Symbol):
            other = other.to_expr()
        if isinstance(other, Expr):
            out = dict(self.monomials)
            for m, c in other.monomials.items():
                out[m] = out.get(m, 0) + c
            return Expr(_strip_zeros(out))
        if _is_scalar(other):
            return self + Expr.from_constant(other)
        return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, (Expr, Symbol)) or _is_scalar(other):
            return self + (-other)
        return NotImplemented

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if isinstance(other, Symbol):
            other = other.to_expr()
        if isinstance(other, Expr):
            out: dict = {}
            for ma, ca in self.monomials.items():
                for mb, cb in other.monomials.items():
                    m = _mul_monomials(ma, mb)
                    out[m] = out.get(m, 0) + ca * cb
            return Expr(_strip_zeros(out))
        if _is_scalar(other):
            if other == 0:
                return Expr.from_constant(0)
            return Expr({m: c * other for m, c in self.monomials.items()})
        if _is_operator(other):
            return _scale(other, self)
        return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)  # coefficients commute with operators and numbers

    def __truediv__(self, other):
        if _is_scalar(other):
            if other == 0:
                raise ZeroDivisionError
            return Expr({m: c / other for m, c in self.monomials.items()})
        return NotImplemented

    def __pow__(self, n: int):
        if not isinstance(n, numbers.Integral) or n < 0:
            raise ValueError("Expr ** n requires non-negative integer n")
        out = Expr.from_constant(1)
        for _ in range(int(n)):
            out = out * self
        return out

    def conjugate(self) -> "Expr":
        return Expr(_strip_zeros({m: c.conjugate() for m, c in self.monomials.items()}))

    def bind(self, values=None, **kwargs):
        """Substitute symbols with numbers; a number when all are bound, else an Expr.

        Keys may be ``Symbol`` instances or names; keyword arguments are names.
        """
        sym_values = _symbol_values(values, kwargs)
        result_terms: dict = {}
        for mon, coeff in self.monomials.items():
            num_val = coeff
            remaining: list[tuple[Symbol, int]] = []
            for sym, exp in mon:
                if sym in sym_values:
                    num_val *= sym_values[sym] ** exp
                else:
                    remaining.append((sym, exp))
            key = _normalise_monomial(remaining)
            result_terms[key] = result_terms.get(key, 0) + num_val
        stripped = _strip_zeros(result_terms)
        if not stripped:
            return 0.0
        if list(stripped) == [frozenset()]:
            v = stripped[frozenset()]
            if isinstance(v, complex) and v.imag == 0:
                return v.real
            return v
        return Expr(stripped)

    def __repr__(self) -> str:
        if not self.monomials:
            return "0"
        parts = []
        for mon, c in sorted(self.monomials.items(), key=lambda kv: sorted((s.name, e) for s, e in kv[0])):
            cs = _format_scalar(c)
            if not mon:
                parts.append(cs)
                continue
            mon_str = "*".join(s.name if e == 1 else f"{s.name}**{e}" for s, e in sorted(mon, key=lambda x: x[0].name))
            parts.append(mon_str if cs == "1" else f"-{mon_str}" if cs == "-1" else f"{cs}*{mon_str}")
        return " + ".join(parts).replace("+ -", "- ")


def _symbol_values(values, kwargs) -> dict:
    out: dict = {}
    for key, value in {**(values or {}), **kwargs}.items():
        if not _is_scalar(value):
            raise TypeError(f"the value of {key} must be a number, not {type(value).__name__}")
        out[key if isinstance(key, Symbol) else Symbol(str(key))] = value
    return out


def _as_scalar(c):
    if _is_scalar(c):
        return c
    raise TypeError(f"Expected a scalar, got {type(c).__name__}")


def _strip_zeros(d):
    return {m: c for m, c in d.items() if c != 0}


def _format_scalar(c) -> str:
    if isinstance(c, complex):
        if c.imag == 0:
            return _format_scalar(c.real)
        return f"({c.real:+g}{c.imag:+g}j)"
    if isinstance(c, float):
        return f"{c:g}"
    return str(c)


def symbols(*names: str) -> tuple[Symbol, ...]:
    """One :class:`Symbol` per name."""
    return tuple(Symbol(n) for n in names)


__all__ = ["Symbol", "Expr", "symbols"]
