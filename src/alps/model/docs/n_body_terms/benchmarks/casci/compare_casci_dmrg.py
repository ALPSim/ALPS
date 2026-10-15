"""Molecular Hamiltonians, full (pq|rs), QDK CASCI vs ALPS dmrg with n-site terms.

    ALPS_BIN=<alps>/_build/default/bin ALPS_XML_PATH=<alps>/src/alps/resources \
    python compare_casci_dmrg.py        # needs qdk-chemistry

The same molecules and ALPS input as compare_casci_nsite.py, solved with dmrg
instead of fulldiag. The orbitals are the DMRG sites, in QDK's order; every
(pq|rs) is a NSITETERM on 2, 3 or 4 of them. With enough states DMRG is
exact for these sizes, so the ground state must match CASCI. H3+ (3 orbitals)
is left out: ALPS dmrg needs at least 4 sites.
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import compare_casci_nsite as ref                     # noqa: E402
from qdk2alps import parm_text, to_alps               # noqa: E402
from qdk_alps_ed import _tool                         # noqa: E402
from qdk_chemistry.algorithms import create           # noqa: E402

TOL = 1e-6
_ENERGY = re.compile(r'<SCALAR_AVERAGE name="Energy"><MEAN>([^<]+)</MEAN>')


def alps_dmrg(ham, na, nb, states=200, sweeps=4, mps=None):
    """Ground-state energy. With mps={} also SAVE_MPS=1 (ALPSim/ALPS#216), and the MPS read into it."""
    tmp = tempfile.mkdtemp(prefix="alps_dmrg_")
    try:
        out = to_alps(ham, tmp, name="qdk-dmrg")
        with open(out["parm"], "a") as fh:
            p = {"CONSERVED_QUANTUMNUMBERS": "Nup,Ndown", "Nup_total": na, "Ndown_total": nb,
                 "SWEEPS": sweeps, "MAXSTATES": states, "LANCZOS_TOLERANCE": 1e-12}
            if mps is not None:
                p["SAVE_MPS"] = 1
            fh.write(parm_text(p))
            fh.write("{ }\n")
        base = os.path.basename(out["parm"])
        subprocess.run([_tool("parameter2xml"), base], cwd=tmp, check=True, capture_output=True)
        subprocess.run([_tool("dmrg"), "--write-xml", base + ".in.xml"], cwd=tmp,
                       check=True, capture_output=True)
        task = next(f for f in os.listdir(tmp) if f.endswith(".out.xml") and ".task" in f)
        with open(os.path.join(tmp, task)) as fh:
            levels = sorted(float(x) for x in _ENERGY.findall(fh.read()))
        if mps is not None:
            sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "alps_checks"))
            from alpsrun import _read_mps
            mps.update(_read_mps(os.path.join(tmp, task.replace(".out.xml", ".out.h5"))))
        return levels[0] + float(ham.get_core_energy())
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    casci = create("multi_configuration_calculator")
    print(f"{'case':30}{'sector':>8}{'QDK CASCI':>16}{'ALPS dmrg':>16}{'|dE|':>9}{'t_dmrg':>8}")
    print("-" * 87)
    bad = 0
    for label, ham, sectors in ref.cases():
        if label.startswith("H3+"):
            continue
        for na, nb in sectors:
            e_qdk, _ = casci.run(ham, na, nb)
            t0 = time.time()
            e_dmrg = alps_dmrg(ham, na, nb)
            dt = time.time() - t0
            ok = abs(e_qdk - e_dmrg) < TOL
            bad += not ok
            print(f"{label:30}{f'({na},{nb})':>8}{e_qdk:16.10f}{e_dmrg:16.10f}"
                  f"{abs(e_qdk - e_dmrg):9.1e}{dt:7.1f}s" + ("" if ok else "  <-- MISMATCH"))
    print("-" * 87)
    print("ALL AGREE" if not bad else f"{bad} MISMATCH(ES)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
