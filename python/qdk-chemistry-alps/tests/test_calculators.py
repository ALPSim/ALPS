# Copyright (c) 2026 Tobias Wolf
# SPDX-License-Identifier: MIT
#
# Adapted from skilledwolf/alpsx (private repository), branch
# agent/v3-execution-and-api-audit at d37e51227,
# python/qdk-chemistry-alps/tests/test_calculators.py: the exact
# diagonalization cases. The MPS DMRG cases are not ported (alps_dmrg is not
# part of this package); test_classic_calculator.py covers the classic
# backend. Changes: the "engine" setting is gone; a determinant-sign check
# against MACIS is added.
"""The exact diagonalization calculator as qdk-chemistry sees it.

The energy it returns, the wavefunction container it packages it in and the
settings that steer it.
"""

from __future__ import annotations

import numpy as np

import qdk_chemistry_alps  # noqa: F401  (registers the backends)


def _calculator(name: str, **settings):
    from qdk_chemistry.algorithms import create

    calc = create("multi_configuration_calculator", name)
    for key, value in settings.items():
        calc.settings().set(key, value)
    return calc


def _exact_energy(qdk_h, n_alpha: int, n_beta: int) -> float:
    """Reference energy from qdk-chemistry's own CAS solver."""
    return float(_calculator("macis_cas").run(qdk_h, n_alpha, n_beta)[0])


def _coefficients(wavefunction) -> dict[str, float]:
    return {
        det.to_string(): coefficient
        for det, coefficient in zip(
            wavefunction.get_active_determinants(),
            wavefunction.get_coefficients(), strict=True,
        )
    }


def test_ed_matches_macis_on_a_hubbard_chain(hubbard_qdk):
    qdk_h = hubbard_qdk(6, u=4.0)
    energy, wavefunction = _calculator("alps_ed").run(qdk_h, 3, 3)

    assert abs(energy - _exact_energy(qdk_h, 3, 3)) < 1e-8
    assert wavefunction.get_container_type() == "state_vector"


def test_ed_ci_vector_matches_macis_up_to_global_phase(hubbard_qdk):
    """Determinant labels and signs agree with qdk-chemistry's own solver."""
    qdk_h = hubbard_qdk(4, u=4.0)
    ours = _coefficients(_calculator("alps_ed").run(qdk_h, 2, 2)[1])
    reference = _coefficients(_calculator("macis_cas").run(qdk_h, 2, 2)[1])

    significant = {det for det, c in reference.items() if abs(c) > 1e-8}
    assert significant <= ours.keys()
    overlap = sum(np.conj(reference[det]) * ours[det] for det in significant)
    assert abs(abs(overlap) - 1.0) < 1e-9


def test_ed_settings_are_declared():
    settings = _calculator("alps_ed").settings()
    for key in ("integral_cutoff", "ci_coefficient_cutoff"):
        assert settings.get(key) is not None


def test_ci_cutoff_shortens_the_expansion(hubbard_qdk):
    qdk_h = hubbard_qdk(4, u=4.0)
    full = _calculator("alps_ed").run(qdk_h, 2, 2)[1]
    trimmed = _calculator("alps_ed", ci_coefficient_cutoff=0.05).run(qdk_h, 2, 2)[1]
    assert len(trimmed.get_coefficients()) < len(full.get_coefficients())
    assert abs(np.linalg.norm(trimmed.get_coefficients()) - 1.0) < 1e-12


def test_spin_imbalanced_sector(hubbard_qdk):
    qdk_h = hubbard_qdk(3, u=4.0)
    energy, _ = _calculator("alps_ed").run(qdk_h, 2, 1)
    assert abs(energy - _exact_energy(qdk_h, 2, 1)) < 1e-8
