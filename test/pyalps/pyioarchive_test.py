 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 #                                                                                 #
 # ALPS Project: Algorithms and Libraries for Physics Simulations                  #
 #                                                                                 #
 # ALPS Libraries                                                                  #
 #                                                                                 #
 # Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   #
 #                                                                                 #
 # ALPS Project: https://alps.comp-phys.org/                                       #
 # SPDX-License-Identifier: MIT                                                    #
 #                                                                                 #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Archive handle lifecycle: alps::hdf5::archive reference-counts one
# archivecontext per (filename, flags), so re-opening a file that already has a
# live handle shares that context. This exercises the overlap sequence -- open
# for write, keep a reader alive across two further write handles, then drop the
# reader last -- which is the only coverage in the suite for that path.
#
# Renamed from pyioarchive.py so pytest collects it; under the retired ctest
# runner the file name did not have to match a collection pattern.

import os
import tempfile

import pyalps.hdf5 as h5


def test_overlapping_archive_handles():
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "blubb")

        o = h5.archive(path, 'w')
        o['/a'] = 0
        del o

        i = h5.archive(path, 'r')
        o = h5.archive(path, 'w')
        o['/a'] = 0
        del o

        del i
        o = h5.archive(path, 'w')
        o['/a'] = 0
        del o

        with h5.archive(path, 'r') as ar:
            assert ar['/a'] == 0


if __name__ == '__main__':
    test_overlapping_archive_handles()
    print("SUCCESS")
