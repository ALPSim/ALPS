"""Molecules: QDK CASCI wave function vs the MPS of ALPS dmrg, full (pq|rs).

    ALPS_BIN=<dmrg with SAVE_MPS> ALPS_XML_PATH=<alps>/src/alps/resources \
    python compare_casci_mps.py        # needs qdk-chemistry, h5py

dmrg writes its ground state as an MPS with SAVE_MPS=1 (ALPSim/ALPS#216), so
ALPS_BIN must point to a dmrg built with that change and the n-site terms.
The input is the one compare_casci_dmrg.py uses: orbitals are sites, every
(pq|rs) a NSITETERM. The MPS is contracted to the amplitude of every
determinant and compared with QDK's CI vector:

- the determinant order of QDK (all alpha modes, then all beta) and ALPS
  (site-major, up before down on a site) differ by blocked_sign, as for the
  ED eigenvector in qdk_alps_ed.py;
- |<CASCI|MPS>| must be 1 (to 1e-8) and each CI coefficient must agree (to
  1e-4, set by DMRG convergence; a wrong sign gives errors of ~0.1).

H3+ is left out: dmrg needs at least 4 sites.
"""

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "alps_checks"))

import compare_casci_nsite as ref                     # noqa: E402
import mpsref                                         # noqa: E402
from compare_casci_dmrg import alps_dmrg              # noqa: E402
from qdk_alps_ed import blocked_sign                  # noqa: E402
from qdk_chemistry.algorithms import create           # noqa: E402

TOL = 1e-8
_LOCAL = {"0": (0, 0), "u": (1, 0), "d": (0, 1), "2": (1, 1)}


def amplitudes(mps, dets):
    """MPS amplitude of each determinant (QDK occupation strings), QDK's sign convention."""
    psi = mpsref.contract(mps["tensors"])
    psi = psi / np.linalg.norm(psi)
    basis = mps["site_basis"][mps["site_type"][0]]
    index = {(int(u), int(d)): k for k, (u, d) in enumerate(zip(basis["Nup"], basis["Ndown"]))}
    d = len(basis["Nup"])
    n = len(mps["tensors"])
    out = []
    for occ in dets:
        # QDK pads the occupation string past the active orbitals
        assert set(occ[n:]) <= {"0"}, occ
        occ = occ[:n]
        flat = 0
        for ch in occ:
            flat = flat * d + index[_LOCAL[ch]]
        out.append(blocked_sign(occ) * psi[flat])
    covered = float(np.sum(np.square(out)))
    return np.array(out), covered, n


def main():
    casci = create("multi_configuration_calculator")
    print(f"{'case':30}{'sector':>8}{'dets':>6}{'|<CASCI|MPS>|':>16}{'max|dc|':>10}"
          f"{'MPS weight':>12}{'E_dmrg-E_CASCI':>16}{'t':>7}")
    print("-" * 105)
    bad = 0
    for label, ham, sectors in ref.cases():
        if label.startswith("H3+"):
            continue
        for na, nb in sectors:
            e_qdk, wf = casci.run(ham, na, nb)
            dets = [str(c) for c in wf.get_active_determinants()]
            c = np.asarray(wf.get_coefficients(), float)
            c /= np.linalg.norm(c)
            t0 = time.time()
            mps = {}
            # 8 sweeps: coefficients converge as the square root of the energy
            e_dmrg = alps_dmrg(ham, na, nb, sweeps=8, mps=mps)
            dt = time.time() - t0
            a, covered, _ = amplitudes(mps, dets)
            ov = float(np.dot(c, a))
            dc = float(np.max(np.abs(np.sign(ov) * a - c)))
            # Coefficients converge as the square root of the DMRG energy error,
            # which varies from run to run (~1e-13 to 1e-10 here), so allow 1e-4;
            # a wrong sign anywhere moves them by ~0.1.
            ok = abs(1 - abs(ov)) < TOL and dc < 1e-4 and abs(1 - covered) < TOL and abs(e_dmrg - e_qdk) < 1e-6
            bad += not ok
            print(f"{label:30}{f'({na},{nb})':>8}{len(dets):>6}{abs(ov):16.12f}{dc:10.1e}"
                  f"{covered:12.9f}{e_dmrg - e_qdk:16.1e}{dt:6.1f}s" + ("" if ok else "  <-- MISMATCH"),
                  flush=True)
    print("-" * 105)
    print("ALL AGREE" if not bad else f"{bad} MISMATCH(ES)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
