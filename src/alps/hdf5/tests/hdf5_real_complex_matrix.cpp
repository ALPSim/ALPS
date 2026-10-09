/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Michele Dolfi <dolfim@phys.ethz.ch>                *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5/matrix.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/numeric/matrix.hpp>

#include <boost/filesystem.hpp>

#include <vector>
#include <complex>
#include <iostream>

TEST(Hdf5, RealComplexMatrixRejectsMissingComplexMetadata) {
    alps::testing::TemporaryDirectory temporary;
    alps::numeric::matrix<double> source(4, 4, 1.5);
    { alps::hdf5::archive archive((temporary.path() / "real_complex.h5").string(), "w"); archive["/matrix"] << source; }
    alps::hdf5::archive archive((temporary.path() / "real_complex.h5").string(), "r");
    alps::numeric::matrix<std::complex<double>> destination;
    EXPECT_THROW(archive["/matrix"] >> destination, alps::hdf5::archive_error);
}
