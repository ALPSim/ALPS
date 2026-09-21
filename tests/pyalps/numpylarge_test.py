# ****************************************************************************
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# ALPS Libraries
#
# Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>
#               2026       by the ALPS collaboration
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
#
# ****************************************************************************

# Formerly numpylarge.py: it doubled an array from 2**10 up to 2**29 doubles
# (4 GiB for the final write alone), printed the sizes, asserted nothing, and
# left its foo%d.h5 files in the working directory. pytest never collected it,
# so it also never cleaned up after itself.
#
# What it was actually exercising is worth keeping: the split-file archive
# ("foo%d.h5" plus mode "l"), which rolls over into foo0.h5, foo1.h5, ... once
# a member file passes the driver's threshold, and large contiguous writes
# through the numpy save path. Both are covered below at a size that belongs in
# a test suite, with assertions and in a temporary directory.

import glob
import os

import numpy as np
import pytest

import pyalps.hdf5 as hdf5


# 2**20 doubles = 8 MiB, enough to exercise a chunked write without making the
# suite slow. The original's 2**29 ceiling tested the operating system, not ALPS.
MAX_EXPONENT = 20


@pytest.mark.parametrize("mode", ["a", "al"])
def test_growing_arrays_roundtrip(tmp_path, mode):
    """Doubling writes must read back byte-for-byte, plain or split-file."""
    pattern = str(tmp_path / "foo%d.h5") if "l" in mode else str(tmp_path / "foo.h5")

    written = {}
    with hdf5.archive(pattern, mode) as archive:
        size = 2 ** 10
        while size <= 2 ** MAX_EXPONENT:
            values = np.linspace(0.0, 1.0, size)
            archive[str(size)] = values
            written[str(size)] = values
            size *= 2

    with hdf5.archive(pattern, "r" + ("l" if "l" in mode else "")) as archive:
        for key, values in written.items():
            restored = archive[key]
            assert restored.shape == values.shape, key
            assert restored.dtype == values.dtype, key
            np.testing.assert_array_equal(restored, values)


def test_split_archive_creates_numbered_member_files(tmp_path):
    """Mode "l" must lay the data out across foo0.h5, foo1.h5, ... ."""
    pattern = str(tmp_path / "foo%d.h5")

    with hdf5.archive(pattern, "al") as archive:
        archive["/block"] = np.zeros(2 ** MAX_EXPONENT)

    members = sorted(glob.glob(str(tmp_path / "foo*.h5")))
    assert members, "the split-file driver wrote no member files"
    assert os.path.basename(members[0]) == "foo0.h5"

    with hdf5.archive(pattern, "rl") as archive:
        np.testing.assert_array_equal(archive["/block"], np.zeros(2 ** MAX_EXPONENT))


def test_empty_and_single_element_arrays(tmp_path):
    """The size-0 and size-1 edges of the same write path."""
    path = str(tmp_path / "edges.h5")

    with hdf5.archive(path, "a") as archive:
        archive["/empty"] = np.empty(0)
        archive["/single"] = np.array([np.pi])

    with hdf5.archive(path, "r") as archive:
        assert archive["/empty"].shape == (0,)
        np.testing.assert_array_equal(archive["/single"], np.array([np.pi]))
