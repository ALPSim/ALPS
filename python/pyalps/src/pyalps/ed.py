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
