# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Exact diagonalization of compiled models (the C++ ``ALPS::ed`` component).

>>> from pyalps import models, ed
>>> model = models.heisenberg(models.lattice.chain(12))
>>> problem = ed.prepare_product(model, sector={"Sz": 0.0})     # dimension 924
>>> result = problem.eigenstates(count=3)
>>> result.eigenvalues[0]                                        # doctest: +SKIP
-5.387390917445...
>>> spectrum = ed.model_spectrum(model, all_momenta=True, count=2)  # Sz and k resolved

``prepare_product``, ``prepare_projected``, ``prepare_finite``,
``model_spectrum`` and ``finite_spectrum`` also take an XML catalog model
(``models.CatalogModel``) in place of a ``Hamiltonian``: its compiled
Hamiltonian runs on its physical basis, in the sector of its ``*_total``
parameters unless ``sector=`` says otherwise (``{}`` is the whole space).
A catalog carries no translations, so momenta need ``translations=``.

Eigenvectors are returned as ``(dimension, count)`` arrays whose column ``j``
belongs to eigenvalue ``j``; sparse matrices carry ``to_scipy()``. Every solve
releases the GIL. ``lanczos(matrix, ...)`` runs the deflated Lanczos driver on
a scipy CSR matrix.
"""

from . import models as _models  # registers the model types first  # noqa: F401

try:
    from ._ext import ed_c as _ed_c
except ImportError as error:  # pragma: no cover - build without PYALPS_BUILD_MODELS
    raise ImportError(
        "pyalps was built without the exact diagonalization bindings (PYALPS_BUILD_MODELS=OFF "
        "or an SDK without ALPS::ed)"
    ) from error

from ._ext.ed_c import *  # noqa: F401,F403

_CatalogModel = getattr(_models, "CatalogModel", None)  # absent without the legacy model library


def _is_catalog(model):
    return _CatalogModel is not None and isinstance(model, _CatalogModel)


def _catalog_inputs(catalog, sector):
    """A catalog model's compiled Hamiltonian, physical basis and sector (its own unless given)."""
    physical = catalog.physical_basis()
    return catalog.compile_hamiltonian(), physical, catalog.sector() if sector is None else dict(sector)


def _check_momentum(momentum, translations, all_momenta=False):
    if (len(momentum) or all_momenta) and not len(translations):
        raise ValueError("a CatalogModel carries no lattice translations: pass translations= "
                         "(SpectrumTranslation values) to resolve momenta")


def _catalog_product(catalog, sector=None, max_full_dimension=10_000_000):
    op, physical, sector = _catalog_inputs(catalog, sector)
    return _ed_c.prepare_product(op, physical, sector, max_full_dimension)


def _catalog_projected(catalog, sector=None, momentum=(), translations=()):
    _check_momentum(momentum, translations)
    op, physical, sector = _catalog_inputs(catalog, sector)
    return _ed_c.prepare_projected(op, projected_spaces(physical), list(translations), sector, list(momentum))


def _catalog_finite(catalog, irrep, sector=None):
    op, physical, sector = _catalog_inputs(catalog, sector)
    return _ed_c.prepare_finite(op, projected_spaces(physical), irrep, sector)


def _catalog_model_spectrum(catalog, sector=None, momentum=(), all_momenta=False, count=None, which="SA",
                            tolerance=None, max_iterations=None, translations=()):
    _check_momentum(momentum, translations, all_momenta)
    op, physical, sector = _catalog_inputs(catalog, sector)
    return _ed_c.model_spectrum(op, projected_spaces(physical), list(translations), sector, list(momentum),
                                all_momenta, count, which, tolerance, max_iterations)


def _catalog_finite_spectrum(catalog, irreps, sector=None, count=None, which="SA", tolerance=None,
                             max_iterations=None):
    op, physical, sector = _catalog_inputs(catalog, sector)
    return _ed_c.finite_spectrum(op, projected_spaces(physical), irreps, sector, count, which, tolerance,
                                 max_iterations)


def _accepting_catalogs(name, catalog_form, signature):
    """``name`` from ed_c, with a CatalogModel first argument routed to ``catalog_form``."""
    native = getattr(_ed_c, name)

    def dispatch(model, *args, **kwargs):
        if _is_catalog(model):
            return catalog_form(model, *args, **kwargs)
        return native(model, *args, **kwargs)

    dispatch.__name__ = dispatch.__qualname__ = name
    dispatch.__doc__ = (f"{name}{signature}\n\nA CatalogModel runs its compiled Hamiltonian on its physical "
                        f"basis; ``sector`` defaults to the catalog's own.\n\nThe other forms:\n\n{native.__doc__}")
    return dispatch


prepare_product = _accepting_catalogs(
    "prepare_product", _catalog_product, "(catalog, sector=None, max_full_dimension=10000000)")
prepare_projected = _accepting_catalogs(
    "prepare_projected", _catalog_projected, "(catalog, sector=None, momentum=(), translations=())")
prepare_finite = _accepting_catalogs("prepare_finite", _catalog_finite, "(catalog, irrep, sector=None)")
model_spectrum = _accepting_catalogs(
    "model_spectrum", _catalog_model_spectrum,
    "(catalog, sector=None, momentum=(), all_momenta=False, count=None, which='SA', tolerance=None, "
    "max_iterations=None, translations=())")
finite_spectrum = _accepting_catalogs(
    "finite_spectrum", _catalog_finite_spectrum,
    "(catalog, irreps, sector=None, count=None, which='SA', tolerance=None, max_iterations=None)")


def lanczos(matrix, count=1, which="SA", tolerance=0.0, max_iterations=0, want_vectors=True):
    """Deflated Lanczos on a scipy sparse matrix (real symmetric or complex Hermitian).

    The matrix is converted to CSR; ``which`` is ``SA``, ``LA`` or ``LM`` and
    ``tolerance``/``max_iterations`` at or below zero select the driver defaults.
    Degenerate eigenvalues come back with multiplicity.
    """
    import scipy.sparse

    csr = scipy.sparse.csr_matrix(matrix)
    if csr.shape[0] != csr.shape[1]:
        raise ValueError("lanczos needs a square matrix")
    return _ed_c.lanczos_csr(csr.indptr, csr.indices, csr.data, csr.shape[0], count, which,
                             tolerance, max_iterations, want_vectors)


__all__ = [name for name in dir(_ed_c) if not name.startswith("_")] + ["lanczos"]
