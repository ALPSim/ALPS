/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5.hpp>
#include <alps/version.h>

#include <boost/filesystem.hpp>

#include <iostream>
#include <string>

TEST(Hdf5Compatibility, HistoricalFortranString) {
    boost::filesystem::path infile(ALPS_TEST_SOURCE_DIR);
    infile /= "hdf5_fortran_string.h5";
    ASSERT_TRUE(boost::filesystem::exists(infile)) << infile;
    alps::hdf5::archive archive(infile);
    EXPECT_TRUE(archive.is_datatype<std::string>("/fortran_string"));
    std::string value;
    archive["/fortran_string"] >> value;
    // The legacy fixed-width reader includes its terminating NUL in the
    // returned std::string. Preserve this observable behavior during migration.
    EXPECT_EQ(value, std::string("N_total\0", 8));
}
