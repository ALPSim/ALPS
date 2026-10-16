# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/src/qdk_chemistry_alps/plugin.py, with the
# alps_classic_dmrg registration of branch WIP/classic-dmrg-chemistry at
# d4384a7f8. Changes: the requirement check imports pyalps.models and
# pyalps.ed instead of alpsx's MPS State API; alps_dmrg (alpsx's MPS DMRG,
# not ported) is not registered; alps_classic_dmrg is registered only when
# pyalps has its classic_dmrg bindings.
"""qdk-chemistry plugin registration for the ALPS backends.

Two registration paths, matching two generations of qdk-chemistry:

* **Released (<= 2.1.0)**: no entry-point discovery exists; plugins register
  by being imported (the bundled-PySCF pattern). Call :func:`load` — or just
  ``import qdk_chemistry_alps``, whose ``__init__`` calls it.
* **Later builds**: the ``qdk_chemistry.plugins`` entry-point group loads
  :class:`AlpsPlugin` automatically at ``import qdk_chemistry``.

Both paths are idempotent against each other, in either import order:
each skips algorithm names the other has already registered, and
``AlpsPlugin`` is materialized lazily (module ``__getattr__``) so this
module never imports ``qdk_chemistry`` at its own import time — doing so
would fire qdk-chemistry's plugin discovery back into this module while
it is still initializing.
"""

from __future__ import annotations

# AlpsPlugin is provided by the module __getattr__ below (lazy import).
__all__ = ["AlpsPlugin", "load"]  # noqa: F822

_loaded = False


def _require_alps() -> None:
    try:
        from pyalps import ed, models  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "qdk-chemistry-alps requires pyalps with its models and ed "
            "bindings (pyalps.models, pyalps.ed). Install the ALPS SDK, then "
            "build python/pyalps of the ALPS repository against it; see "
            "README.md in python/qdk-chemistry-alps."
        ) from exc


def _has_classic_dmrg() -> bool:
    try:
        from pyalps import classic_dmrg  # noqa: F401
    except ImportError:
        return False
    return True


def _calculators():
    """The calculator factories this pyalps build supports, by name."""
    from .ed_calculator import AlpsEdCalculator

    factories = {"alps_ed": AlpsEdCalculator}
    if _has_classic_dmrg():
        from .classic_dmrg_calculator import AlpsClassicDmrgCalculator

        factories["alps_classic_dmrg"] = AlpsClassicDmrgCalculator
    return factories


def _register(register) -> None:
    from qdk_chemistry.algorithms import available

    registered = available("multi_configuration_calculator")
    for name, factory in _calculators().items():
        if name not in registered:
            register(factory)


def load() -> None:
    """Register the ALPS algorithms with qdk-chemistry (idempotent)."""
    global _loaded
    if _loaded:
        return
    _loaded = True

    _require_alps()

    from qdk_chemistry.algorithms import register

    _register(register)


# --- entry-point path (qdk-chemistry builds with plugin discovery) --------

def _build_entry_point_plugin():
    from qdk_chemistry.plugins import PluginRegistrar, QdkChemistryPlugin

    class AlpsPlugin(QdkChemistryPlugin):
        """Entry-point plugin for qdk-chemistry's automatic discovery."""

        api_version = 1

        def register(self, registrar: PluginRegistrar) -> None:
            _require_alps()
            # Entry-point discovery imports this module, whose package
            # __init__ has already run the import-time registration path —
            # skip names that are present so the two paths stay idempotent.
            _register(registrar.register_algorithm)
            global _loaded
            _loaded = True

    return AlpsPlugin


_plugin_cls = None


def __getattr__(name: str):
    global _plugin_cls
    if name == "AlpsPlugin":
        if _plugin_cls is None:
            _plugin_cls = _build_entry_point_plugin()
        return _plugin_cls
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
