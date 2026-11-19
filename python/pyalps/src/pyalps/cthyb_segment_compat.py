# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
"""Run parameters written for ``pyalps.cthyb`` on CT-HYB-SEGMENT.

The segment solver (``pyalps.cthyb_segment``, the ct-hyb-segment program)
descends from the same code as the hybridization solver (``pyalps.cthyb``)
but reads its parameters differently:

* its own options carry a ``cthyb.`` prefix, and N_TAU, N_MATSUBARA and
  N_ORBITALS are called N, NMATSUBARA and FLAVORS;
* a sweep is one update, not N_MEAS updates, so SWEEPS and THERMALIZATION
  are multiplied by N_MEAS;
* it adds U/2 to MU, which the hybridization solver takes as is;
* the inter-orbital interaction is Uprime (default 0), not U' (default U-2J);
* it measures G(i omega_n) unless told otherwise (MEASURE_freq), the
  hybridization solver only when asked to.

``from_hybridization`` translates a parameter dictionary accordingly, and
``solve`` runs the translated parameters, so that a script written for
``pyalps.cthyb`` can use either solver.
"""

_RENAMED = {"N_TAU": "N", "N_MATSUBARA": "NMATSUBARA", "N_ORBITALS": "FLAVORS"}
_SWEEP_COUNTS = ("SWEEPS", "THERMALIZATION")
_PREFIXED = (
    "N_MEAS", "N_HISTOGRAM_ORDERS", "N_LEGENDRE", "N_W", "N_w2", "N_nn",
    "MEASURE_freq", "MEASURE_time", "MEASURE_legendre", "MEASURE_g2w",
    "MEASURE_h2w", "MEASURE_nn", "MEASURE_nnt", "MEASURE_nnw",
    "MEASURE_sector_statistics", "COMPUTE_VERTEX", "SPINFLIP", "GLOBALFLIP",
    "ACCURATE_COVARIANCE", "TEXT_OUTPUT", "DELTA", "DELTA_IN_HDF5",
    "RET_INT_K", "K_IN_HDF5", "BASEPATH", "DMFT_FRAMEWORK", "J")


def from_hybridization(parms):
    """Return the CT-HYB-SEGMENT parameters equivalent to ``pyalps.cthyb`` ones."""
    result = {}
    for key, value in parms.items():
        if key in _RENAMED:
            result[_RENAMED[key]] = value
        elif key in _SWEEP_COUNTS:
            result["cthyb." + key] = int(value) * int(parms.get("N_MEAS", 1))
        elif key in _PREFIXED:
            result["cthyb." + key] = value
        elif key not in ("MU", "U'"):
            result[key] = value

    result.setdefault("cthyb.MEASURE_freq", False)
    if "U_MATRIX" not in parms and "U" in parms:
        J = parms.get("J", 0.0)
        result["Uprime"] = parms.get("U'", parms["U"] - 2 * J)
    if "MU" in parms and "MU_VECTOR" not in parms:
        if "U" not in parms or "U_MATRIX" in parms:
            raise ValueError("CT-HYB-SEGMENT reads MU relative to U/2; with U_MATRIX, "
                             "give the orbital energies as MU_VECTOR instead")
        result["MU"] = parms["MU"] - parms["U"] / 2.0
    return result


def solve(parms):
    """Solve the impurity problem of ``pyalps.cthyb`` parameters with CT-HYB-SEGMENT."""
    from pyalps import cthyb_segment
    cthyb_segment.solve(from_hybridization(parms))
