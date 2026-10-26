/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 */

// Adapted from ALPSCore 7146b9e1 hdf5/test/hdf5_open_strmode.cpp, which tests
// the mode validation from ALPSCore 33e26902. Changes: ALPS's temporary
// directory; "wl" is valid in ALPS (split files, which ALPSCore removed), so
// the invalid examples are "wx" and "R"; OpenModesCombined and
// OpenModesCharacter are new.

#include <alps/hdf5/archive.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>

#include <string>

/// Test for "a","w","r" open modes
TEST(hdf5,OpenModesAWR) {
    alps::testing::TemporaryDirectory temporary;
    const std::string h5name=(temporary.path() / "hdf5_open.h5").string();
    {
        ASSERT_THROW(alps::hdf5::archive ar(h5name,"r"),  alps::hdf5::archive_not_found);
        ASSERT_THROW(alps::hdf5::archive ar2(h5name,"r"), alps::hdf5::archive_not_found);
    }
    {
        alps::hdf5::archive ar(h5name,"w"); // same as "a"
        ar["/int1"]=1;
    }
    {
        alps::hdf5::archive ar(h5name,"w"); // same as "a"
        ar["/int2"]=2;
    }
    {
        alps::hdf5::archive ar(h5name,"a");
        ar["/int3"]=3;
    }
    {
        alps::hdf5::archive ar(h5name,"r");
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
        alps::hdf5::archive ar(h5name,"");
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
        alps::hdf5::archive ar(h5name);
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

/// Test for incorrect open modes
TEST(hdf5,OpenModesWrong) {
    alps::testing::TemporaryDirectory temporary;
    const std::string h5name=(temporary.path() / "hdf5_open.h5").string();
    {
      alps::hdf5::archive ar(h5name,"w");
      ar["int"]=1;
    }
    {
      EXPECT_THROW(alps::hdf5::archive ar(h5name,"wx"), alps::hdf5::wrong_mode);
    }
    {
      EXPECT_THROW(alps::hdf5::archive ar(h5name,"R"), alps::hdf5::wrong_mode);
    }
    {
      EXPECT_THROW(alps::hdf5::archive ar(h5name,"no_such_mode"), alps::hdf5::wrong_mode);
    }
}

/// ALPS keeps the "l" (split file) and "m" (memory) letters and combinations.
TEST(hdf5,OpenModesCombined) {
    alps::testing::TemporaryDirectory temporary;
    const std::string h5name=(temporary.path() / "hdf5_open.h5").string();
    {
      alps::hdf5::archive ar(h5name,"wc");
      ar["/int"]=1;
    }
    {
      alps::hdf5::archive ar(h5name,"ra");
      int n=0;
      ar["/int"] >> n;
      EXPECT_EQ(1,n);
    }
}


/// The single-letter constructors reject unknown letters like the string modes.
TEST(hdf5,OpenModesCharacter) {
    alps::testing::TemporaryDirectory temporary;
    const std::string h5name=(temporary.path() / "hdf5_open.h5").string();
    {
      alps::hdf5::archive ar(h5name,'w');
      ar["/int"]=1;
    }
    EXPECT_THROW(alps::hdf5::archive ar(h5name,'x'), alps::hdf5::wrong_mode);
    EXPECT_THROW(alps::hdf5::archive ar(h5name,static_cast<signed char>('x')), alps::hdf5::wrong_mode);
    {
      alps::hdf5::archive ar(h5name,static_cast<signed char>('a'));
      ar["/int2"]=2;
    }
    alps::hdf5::archive ar(h5name,'r');
    int n=0;
    ar["/int2"] >> n;
    EXPECT_EQ(2,n);
    EXPECT_THROW(ar["/int3"] << 3, alps::hdf5::archive_error);
}
