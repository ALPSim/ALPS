# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Programmatic lattices and models (the C++ ``ALPS::models`` kit).

``pyalps.models.lattice`` holds the unit cells, the named lattices and the
legacy XML bridge; this module holds the local bases, the operator algebra,
``Hamiltonian`` with the built-in model factories, compiled operators, physical
bases, the qubit/Pauli export and the XML catalog importer ``CatalogModel``;
symbolic coefficients (``Symbol``, ``ParametricHamiltonian``), ``validate``
and the ``has_sign_problem`` check.

>>> from pyalps import models
>>> lat = models.lattice
>>> model = models.heisenberg(lat.chain(8), J=1.0)
>>> compiled = models.compile_hamiltonian(model)
>>> compiled.op.site_dims, len(compiled.op.terms)
([2, 2, 2, 2, 2, 2, 2, 2], 24)

The bindings are optional: they are present when pyalps was built with
``PYALPS_BUILD_MODELS`` (the default) against an SDK with the ``models``
component.
"""

import sys as _sys

try:
    from .._ext import models_c as _models_c
except ImportError as error:  # pragma: no cover - build without PYALPS_BUILD_MODELS
    raise ImportError(
        "pyalps was built without the models bindings (PYALPS_BUILD_MODELS=OFF "
        "or an SDK without ALPS::models)"
    ) from error

from .._ext.models_c import *  # noqa: F401,F403

lattice = _models_c.lattice
# ``import pyalps.models.lattice`` and ``from pyalps.models.lattice import chain``
# need a sys.modules entry; the attribute alone would not serve them.
_sys.modules[__name__ + ".lattice"] = lattice

_CatalogModel = getattr(_models_c, "CatalogModel", None)  # absent without the legacy model library


# -- operator sugar -------------------------------------------------------------------
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227: Site.__getattr__ and
# OpFactory of python/alps/src/alps/models/{sites,operators}.py, here over
# the bound SiteBasis and TermSum.located; site["name"] added.

class LocatedOperator:
    """``site.Sz`` or ``site["Sz"]``: called with a site index, ``located(site, "Sz", i)``."""

    __slots__ = ("basis", "name")

    def __init__(self, basis, name):
        self.basis, self.name = basis, name

    def __call__(self, site):
        return TermSum.located(self.basis, self.name, site)

    def __repr__(self):
        return f"LocatedOperator({self.name!r}, {self.basis!r})"


def _site_getattr(self, name):
    # Reached only when normal attribute lookup misses.
    if name.startswith("_") or not self.has_op(name):
        raise AttributeError(f"SiteBasis has no attribute or operator {name!r}")
    return LocatedOperator(self, name)


def _site_getitem(self, name):
    if not isinstance(name, str) or not self.has_op(name):
        raise KeyError(name)
    return LocatedOperator(self, name)


SiteBasis.__getattr__ = _site_getattr
SiteBasis.__getitem__ = _site_getitem


def compile_observable(model, observable):
    """``observable`` as a ``CompiledOperator`` on the sites of ``model``.

    ``model`` is a ``Hamiltonian``, a ``CatalogModel`` or a ``PhysicalBasis``.
    ``observable`` is a ``CompiledOperator`` (returned unchanged), a
    ``TermSum`` of located terms (each factor's site basis must match the
    model's at its site, so terms built on another local state order are
    rejected), or an operator name summed over every site that defines it
    (on a catalog model, ``CatalogModel.compile_operator(name)``).
    """
    require_bound(model)
    if isinstance(observable, CompiledOperator):
        return observable
    if isinstance(observable, TermSum):
        basis = model if isinstance(model, PhysicalBasis) else model.physical_basis()
        return compile_operator(observable, basis)
    if isinstance(observable, str):
        if _CatalogModel is not None and isinstance(model, _CatalogModel):
            return model.compile_operator(observable)
        if isinstance(model, Hamiltonian):
            total = TermSum()
            for site in range(model.num_sites):
                basis = model.site_for(site)
                if basis.has_op(observable):
                    total += TermSum.located(basis, observable, site)
            if total.empty:
                raise ValueError(f"no site of the model defines an operator {observable!r}")
            return compile_operator(total, model.physical_basis())
        raise TypeError("an operator name needs a Hamiltonian or CatalogModel to resolve it")
    raise TypeError(f"cannot compile an observable of type {type(observable).__name__}")


# -- matrices -------------------------------------------------------------------------

def _to_sparse_matrix(self, sector=None):
    """The Hamiltonian as a scipy CSR matrix on its product basis or one charge sector.

    Rows and columns follow ``pyalps.ed.prepare_product(model, sector).states``;
    ``sector=None`` is the whole space for a ``Hamiltonian`` and the
    ``*_total`` sector for a ``CatalogModel``.
    """
    from .. import ed

    problem = ed.prepare_product(self) if sector is None else ed.prepare_product(self, sector=sector)
    return problem.hamiltonian.to_scipy()


def _to_dense_matrix(self, sector=None):
    """The same as a dense NumPy array (see ``to_sparse_matrix``)."""
    return _to_sparse_matrix(self, sector).toarray()


for _model_type in (Hamiltonian, _CatalogModel):
    if _model_type is not None:
        _model_type.to_sparse_matrix = _to_sparse_matrix
        _model_type.to_dense_matrix = _to_dense_matrix
del _model_type

# -- symbolic coefficients and validation ---------------------------------------------

from .symbols import Expr, Symbol, symbols  # noqa: E402
from .parametric import ParametricHamiltonian, SymbolicTermSum, require_bound  # noqa: E402
from .validate import Diagnostic, ValidationReport, validate  # noqa: E402

Hamiltonian.validate = lambda self, target=None: validate(self, target)
ParametricHamiltonian.validate = lambda self, target=None: validate(self, target)

__all__ = [name for name in dir(_models_c) if not name.startswith("_")] + [
    "LocatedOperator", "compile_observable", "Symbol", "Expr", "symbols", "SymbolicTermSum",
    "ParametricHamiltonian", "Diagnostic", "ValidationReport", "validate"]
