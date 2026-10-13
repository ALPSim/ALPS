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

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/utility/vectorio.hpp>

#include <boost/filesystem.hpp>

#include <vector>
#include <complex>
#include <iostream>

template <class T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const& v)
{
	os << "[" <<alps::write_vector(v, " ", 6) << "]";
	return os;
}

TEST(Hdf5, RealComplexVecRejectsMissingComplexMetadata) {
    alps::testing::TemporaryDirectory temporary;
    std::vector<double> source(6, 3.2);
    { alps::hdf5::archive archive((temporary.path() / "real_complex.h5").string(), "w"); archive["/vec"] << source; }
    alps::hdf5::archive archive((temporary.path() / "real_complex.h5").string(), "r");
    std::vector<std::complex<double>> destination;
    EXPECT_THROW(archive["/vec"] >> destination, alps::hdf5::archive_error);
}
