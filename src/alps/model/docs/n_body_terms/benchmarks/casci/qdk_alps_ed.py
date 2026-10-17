"""ALPS exact diagonalisation as a QDK MultiConfigurationCalculator.

Importing this module registers `alps_ed`, after which it is used exactly like
QDK's own CI solver -- same factory, same call, same return type:

    import qdk_alps_ed                                   # registers the backend
    from qdk_chemistry.algorithms import create

    calc  = create("multi_configuration_calculator", "alps_ed")
    E, wf = calc.run(hamiltonian, n_alpha, n_beta)

    wf.get_active_determinants()    # Slater determinants
    wf.get_coefficients()           # CI coefficients

Swapping `"alps_ed"` for the default is the only change a caller makes.

Under the hood: qdk2alps writes an ALPS custom lattice (plus a Hamiltonian
former when ALPS has no matching model), fulldiag diagonalises the sector, and
the eigenvector is translated back into QDK's determinant basis.

Two conversions happen on the way back, both covered in checks/:

* ALPS prints |Nup Ndown> per site; QDK uses an occupation string over orbitals
  ('0', 'u', 'd', '2'). Site index maps to orbital index directly.
* QDK orders modes blocked (all alpha, then all beta), ALPS interleaved
  (site-major), so each determinant picks up (-1)**(alpha-above-beta pairs).

ALPS binaries are found via $ALPS_BIN, default ~/alps-install/bin.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile

import numpy as np

from qdk_chemistry.algorithms import MultiConfigurationCalculator, available, register
from qdk_chemistry.data import Configuration, StateVectorContainer, Wavefunction

from qdk2alps import parm_text, to_alps

__all__ = ["AlpsEdCalculator", "load", "blocked_sign", "ALPS_BIN"]

ALPS_BIN = os.environ.get("ALPS_BIN", os.path.expanduser("~/alps-install/bin"))

_ENERGY = re.compile(r'<SCALAR_AVERAGE name="Energy"><MEAN>([^<]+)</MEAN>')
_KET = re.compile(r"\|\s*(\d+)\s+(\d+)\s*>")
_BASIS = re.compile(r"Basis:\s*\{(.*?)\}", re.S)
_VEC0 = re.compile(r"Eigenvector# 0:\s*\[(.*?)\]\s*\n", re.S)
_COMPLEX = re.compile(r"\(([-\d.eE+]+),")

_OCC = {(0, 0): "0", (1, 0): "u", (0, 1): "d", (1, 1): "2"}


class AlpsNotFound(RuntimeError):
    pass


def _tool(name: str) -> str:
    path = os.path.join(ALPS_BIN, name)
    if os.path.isfile(path):
        return path
    found = shutil.which(name)
    if found:
        return found
    raise AlpsNotFound(
        f"{name} not found in {ALPS_BIN!r} or on PATH. "
        "Set ALPS_BIN to your ALPS install's bin directory."
    )


def blocked_sign(occ: str) -> int:
    """Sign relating QDK's blocked mode order to ALPS's interleaved one.

    QDK orders all alpha modes then all beta; ALPS orders site-major. Moving a
    beta operator past every alpha operator on a higher orbital costs one
    transposition each.
    """
    A = [i for i, c in enumerate(occ) if c in "u2"]
    B = [i for i, c in enumerate(occ) if c in "d2"]
    return (-1) ** sum(sum(1 for a in A if a > b) for b in B)


def _parse(xml_text: str, stdout: str):
    """(sorted levels, [(occupation string, coefficient)] or None)."""
    levels = sorted(float(x) for x in _ENERGY.findall(xml_text))

    basis = _BASIS.search(stdout)
    vec = _VEC0.search(stdout)
    if not (basis and vec):
        return levels, None

    dets = []
    for line in basis.group(1).strip().splitlines():
        kets = _KET.findall(line)
        if kets:
            dets.append("".join(_OCC[(int(a), int(b))] for a, b in kets))
    coeffs = [float(x) for x in _COMPLEX.findall(vec.group(1))]
    if len(dets) != len(coeffs):
        return levels, None
    return levels, list(zip(dets, coeffs))


class AlpsEdCalculator(MultiConfigurationCalculator):
    """ALPS fulldiag behind QDK's MultiConfigurationCalculator interface."""

    def __init__(self, keep_workdir: bool = False):
        super().__init__()
        self.keep_workdir = keep_workdir
        self.last_levels: list[float] = []
        self.last_workdir: str | None = None

    def name(self) -> str:
        return "alps_ed"

    def run(self, hamiltonian, n_active_alpha_electrons, n_active_beta_electrons):
        na, nb = int(n_active_alpha_electrons), int(n_active_beta_electrons)
        tmp = tempfile.mkdtemp(prefix="alps_ed_")
        self.last_workdir = tmp
        try:
            out = to_alps(hamiltonian, tmp, name="qdk-ed")
            with open(out["parm"], "a") as fh:
                fh.write(parm_text({"CONSERVED_QUANTUMNUMBERS": "Nup,Ndown",
                                    "Nup_total": na,
                                    "Ndown_total": nb,
                                    "PRINT_EIGENVECTORS": 1}))
                fh.write("{ }\n")

            base = os.path.basename(out["parm"])
            subprocess.run([_tool("parameter2xml"), base], cwd=tmp,
                           check=True, capture_output=True)
            proc = subprocess.run([_tool("fulldiag"), "--write-xml", base + ".in.xml"],
                                  cwd=tmp, check=True, capture_output=True, text=True)

            results = [f for f in os.listdir(tmp)
                       if f.endswith(".out.xml") and ".task" in f]
            if not results:
                raise RuntimeError(f"fulldiag produced no task output in {tmp}")
            with open(os.path.join(tmp, results[0])) as fh:
                xml_text = fh.read()

            levels, state = _parse(xml_text, proc.stdout)
            if not levels:
                raise RuntimeError("no energies parsed from fulldiag output")

            # ALPS's model language has no constant term, so E_core never
            # reached it. It shifts every level uniformly.
            core = float(hamiltonian.get_core_energy())
            levels = [e + core for e in levels]
            self.last_levels = levels

            return levels[0], self._wavefunction(hamiltonian, state)
        finally:
            if not self.keep_workdir:
                shutil.rmtree(tmp, ignore_errors=True)
                self.last_workdir = None

    @staticmethod
    def _wavefunction(hamiltonian, state):
        if state is None:
            raise RuntimeError(
                "fulldiag did not print an eigenvector; cannot build a "
                "Wavefunction. This needs PRINT_EIGENVECTORS, which run() sets."
            )
        dets, coeffs = [], []
        for occ, c in state:
            dets.append(Configuration.from_spin_half_string(occ))
            coeffs.append(blocked_sign(occ) * c)

        vec = np.asarray(coeffs, dtype=float)
        norm = np.linalg.norm(vec)
        if norm:
            vec = vec / norm
        # Fix the free global phase the same way every time, so repeated runs
        # and comparisons against QDK are reproducible.
        if vec[int(np.argmax(np.abs(vec)))] < 0:
            vec = -vec

        return Wavefunction(
            StateVectorContainer(vec.reshape(-1, 1), dets, hamiltonian.get_orbitals()))


def load() -> None:
    """Register `alps_ed`. Idempotent; safe to call more than once."""
    if "alps_ed" not in available("multi_configuration_calculator"):
        register(lambda: AlpsEdCalculator())


load()
