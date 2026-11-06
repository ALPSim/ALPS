// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// Enumerations written by other tools read as their base type. h5py stores a
// Python bool as an enum {FALSE=0, TRUE=1} over int8, which the archive used to
// refuse as an unknown datatype.
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/testing/temporary_directory.hpp>

#include <hdf5.h>

#include <string>
#include <vector>

#include <gtest/gtest.h>

TEST(Hdf5Enum, H5pyBoolsReadAsTheirBaseType) {
    alps::testing::TemporaryDirectory temporary;
    const std::string file_name = (temporary.path() / "enum.h5").string();
    {
        // Write a scalar, an array and an attribute the way h5py writes bools.
        hid_t file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hid_t type = H5Tenum_create(H5T_NATIVE_SCHAR);
        signed char value = 0;
        H5Tenum_insert(type, "FALSE", &value);
        value = 1;
        H5Tenum_insert(type, "TRUE", &value);

        hid_t scalar_space = H5Screate(H5S_SCALAR);
        hid_t scalar = H5Dcreate2(file, "/flag", type, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(scalar, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
        hid_t attribute = H5Acreate2(scalar, "attr", type, scalar_space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attribute, type, &value);
        H5Aclose(attribute);
        H5Dclose(scalar);
        H5Sclose(scalar_space);

        hsize_t extent = 3;
        signed char values[3] = {1, 0, 1};
        hid_t array_space = H5Screate_simple(1, &extent, NULL);
        hid_t array = H5Dcreate2(file, "/flags", type, array_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(array, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, values);
        H5Dclose(array);
        H5Sclose(array_space);

        H5Tclose(type);
        H5Fclose(file);
    }

    alps::hdf5::archive ar(file_name, "r");
    EXPECT_TRUE(ar.is_datatype<bool>("/flag"));
    bool flag = false;
    ar["/flag"] >> flag;
    EXPECT_TRUE(flag);
    bool attr = false;
    ar["/flag/@attr"] >> attr;
    EXPECT_TRUE(attr);
    std::vector<int> flags;
    ar["/flags"] >> flags;
    EXPECT_EQ(flags, (std::vector<int>{1, 0, 1}));
}
