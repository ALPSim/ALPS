// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// A parameter group written with h5py may hold a Python bool, stored as an
// enum {FALSE=0, TRUE=1} over int8. It loads as the text "1" or "0", which
// define<bool>() parses.
#include <alps/params.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/testing/temporary_directory.hpp>

#include <hdf5.h>

#include <string>

#include <gtest/gtest.h>

TEST(ParamsHdf5Enum, LoadsAnH5pyBool) {
    alps::testing::TemporaryDirectory temporary;
    const std::string file_name = (temporary.path() / "h5py-parameters.h5").string();
    {
        hid_t file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hid_t group = H5Gcreate2(file, "/parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t type = H5Tenum_create(H5T_NATIVE_SCHAR);
        signed char value = 0;
        H5Tenum_insert(type, "FALSE", &value);
        value = 1;
        H5Tenum_insert(type, "TRUE", &value);
        hid_t space = H5Screate(H5S_SCALAR);
        hid_t flag = H5Dcreate2(group, "FLAG", type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(flag, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
        int sweeps = 1000;
        hid_t sweeps_set = H5Dcreate2(group, "SWEEPS", H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(sweeps_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sweeps);
        H5Dclose(sweeps_set);
        H5Dclose(flag);
        H5Sclose(space);
        H5Tclose(type);
        H5Gclose(group);
        H5Fclose(file);
    }
    alps::params parameters;
    {
        alps::hdf5::archive ar(file_name, "r");
        ar["/parameters"] >> parameters;
    }
    parameters.define<bool>("FLAG", "an h5py bool").define<int>("SWEEPS", "sweeps");
    ASSERT_TRUE(parameters.ok());
    EXPECT_TRUE(parameters["FLAG"].as<bool>());
    EXPECT_EQ(parameters["SWEEPS"].as<int>(), 1000);
}
