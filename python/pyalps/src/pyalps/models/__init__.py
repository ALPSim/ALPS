# Copyright (C) 2026 by the ALPS collaboration
# SPDX-License-Identifier: MIT

"""Programmatic lattices and models (the C++ ``ALPS::models`` kit).

``pyalps.models.lattice`` holds the unit cells, the named lattices and the
legacy XML bridge; this module holds the local bases, the operator algebra,
``Hamiltonian`` with the built-in model factories, compiled operators, physical
bases, the qubit/Pauli export and the XML catalog importer ``CatalogModel``.

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

__all__ = [name for name in dir(_models_c) if not name.startswith("_")]
