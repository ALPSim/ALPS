# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]
                       / 'tutorials' / 'optical-lattice-01-bandstructure'))
import bandstructure  # noqa: E402


def test_reference_values():
    # Values printed by the removed pyalps.dwa.bandstructure for the same lattice.
    t, U = bandstructure.hubbard_parameters([8., 8., 8.], [843., 843., 843.], 114.8, 86.99, 200)
    assert np.allclose(t, 4.77051, atol=1e-4)
    assert abs(U - 38.7018) < 1e-3
    assert np.allclose(U / t, 8.11272, atol=1e-4)


def test_wannier_normalized():
    for V0 in (4., 8., 16.):
        _, w, dx = bandstructure.band_1d(V0, 200, 20)
        assert abs(bandstructure.trapezoid(w**2, dx) - 1) < 1e-9


def test_deeper_lattice():
    # A deeper lattice suppresses hopping and increases the onsite interaction.
    t, U = bandstructure.hubbard_parameters([4., 8., 8.], [843., 843., 843.], 114.8, 86.99, 200)
    t8, U8 = bandstructure.hubbard_parameters([8., 8., 8.], [843., 843., 843.], 114.8, 86.99, 200)
    assert t[0] > t[1] == t[2]
    assert U < U8


if __name__ == '__main__':
    test_reference_values()
    test_wannier_normalized()
    test_deeper_lattice()
