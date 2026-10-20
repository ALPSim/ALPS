/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */

// Adapted from ALPSCore 7146b9e1 hdf5/test/hdf5_open_intmode.cpp. Changes:
// ALPS's temporary directory, and the WRITE|REPLACE case ALPSCore left
// commented out is enabled, since ALPS implements REPLACE.

#include <alps/hdf5/archive.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>

#include <string>

/// Test for "a","w","r" open modes
TEST(hdf5,OpenModesAWR) {
    alps::testing::TemporaryDirectory temporary;
    const std::string h5name=(temporary.path() / "hdf5_open.h5").string();
    {
        ASSERT_THROW(alps::hdf5::archive ar(h5name,alps::hdf5::archive::READ),  alps::hdf5::archive_not_found);
        ASSERT_THROW(alps::hdf5::archive ar2(h5name,alps::hdf5::archive::READ), alps::hdf5::archive_not_found);
    }
    {
        alps::hdf5::archive ar(h5name,alps::hdf5::archive::WRITE); // same as "a"
        ar["/int1"]=1;
    }
    {
        alps::hdf5::archive ar(h5name,alps::hdf5::archive::WRITE); // same as "a"
        ar["/int2"]=2;
    }
    {
        // ALPS implements REPLACE: write to a copy, renamed over the file on close.
        alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE|alps::hdf5::archive::REPLACE);
        ar["/int3"]=3;
    }
    {
        alps::hdf5::archive ar(h5name,alps::hdf5::archive::READ);
        int n=0;
        ar["/int1"] >> n;
        EXPECT_EQ(1,n);
        ar["/int2"] >> n;
        EXPECT_EQ(2,n);
        ar["/int3"] >> n;
        EXPECT_EQ(3,n);
        EXPECT_THROW(ar["/int4"] << 4, alps::hdf5::archive_error);
    }
    {
        alps::hdf5::archive ar(h5name,0);
        int n=0;
        ar["/int1"] >> n;
        EXPECT_EQ(1,n);
        ar["/int2"] >> n;
        EXPECT_EQ(2,n);
        ar["/int3"] >> n;
        EXPECT_EQ(3,n);
        EXPECT_THROW(ar["/int4"] << 4, alps::hdf5::archive_error);
    }
}

/// Test for incorrect and unusual open modes
TEST(hdf5,OpenModesWrong) {
    alps::testing::TemporaryDirectory temporary;
    const std::string h5name=(temporary.path() / "hdf5_open.h5").string();
    {
      alps::hdf5::archive ar(h5name,alps::hdf5::archive::WRITE);
      ar["int"]=1;
    }
    {
      alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE|alps::hdf5::archive::MEMORY);
      SUCCEED();
    }
    {
      EXPECT_THROW(alps::hdf5::archive ar(h5name,9999), alps::hdf5::wrong_mode);
    }
}
