# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
"""Parameters written for pyalps.cthyb, translated for CT-HYB-SEGMENT."""
import pytest

from pyalps.cthyb_segment_compat import from_hybridization


def hybridization_parms(**changes):
    parms = {
        "SWEEPS": 1000, "THERMALIZATION": 10, "N_MEAS": 50, "MAX_TIME": 60,
        "SEED": 0, "N_ORBITALS": 2, "N_TAU": 1000, "N_MATSUBARA": 500,
        "BETA": 45.0, "U": 4.0, "MU": 2.0, "DELTA": "delta.dat",
        "TEXT_OUTPUT": 1, "MEASURE_nnw": 1, "BASENAME": "run",
    }
    parms.update(changes)
    return parms


def test_names_and_prefixes():
    result = from_hybridization(hybridization_parms())
    assert (result["N"], result["NMATSUBARA"], result["FLAVORS"]) == (1000, 500, 2)
    for key in ("N_TAU", "N_MATSUBARA", "N_ORBITALS", "N_MEAS", "DELTA", "TEXT_OUTPUT"):
        assert key not in result
    assert result["cthyb.N_MEAS"] == 50
    assert result["cthyb.DELTA"] == "delta.dat"
    assert result["cthyb.MEASURE_nnw"] == 1
    # Shared and unknown parameters pass through unchanged.
    assert (result["BETA"], result["MAX_TIME"], result["SEED"], result["BASENAME"]) == (45.0, 60, 0, "run")


def test_sweeps_count_single_updates():
    result = from_hybridization(hybridization_parms())
    assert result["cthyb.SWEEPS"] == 1000 * 50
    assert result["cthyb.THERMALIZATION"] == 10 * 50


def test_chemical_potential_relative_to_half_filling():
    # CT-HYB-SEGMENT adds U/2 to MU; the hybridization solver does not.
    assert from_hybridization(hybridization_parms())["MU"] == 0.0
    assert from_hybridization(hybridization_parms(MU=3.0, U=2.0))["MU"] == 2.0


def test_interorbital_interaction_defaults_like_the_hybridization_solver():
    assert from_hybridization(hybridization_parms(J=0.5))["Uprime"] == 4.0 - 2 * 0.5
    assert from_hybridization(hybridization_parms(J=0.5))["cthyb.J"] == 0.5
    result = from_hybridization(hybridization_parms(**{"U'": 3.0}))
    assert result["Uprime"] == 3.0 and "U'" not in result


def test_frequency_measurement_off_unless_requested():
    assert from_hybridization(hybridization_parms())["cthyb.MEASURE_freq"] is False
    assert from_hybridization(hybridization_parms(MEASURE_freq=1))["cthyb.MEASURE_freq"] == 1


def test_mu_vector_and_u_matrix():
    result = from_hybridization(hybridization_parms(MU_VECTOR="mu.dat"))
    # Both solvers take MU_VECTOR as is and ignore MU next to it.
    assert result["MU_VECTOR"] == "mu.dat" and "MU" not in result
    with pytest.raises(ValueError, match="MU_VECTOR"):
        from_hybridization(hybridization_parms(U_MATRIX="u.dat"))
    result = from_hybridization(hybridization_parms(U_MATRIX="u.dat", MU_VECTOR="mu.dat"))
    assert "Uprime" not in result and result["U_MATRIX"] == "u.dat"
