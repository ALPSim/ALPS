/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

// Adapted from ALPSCore 7146b9e1 hdf5/test/hdf5_attributes.cpp, added there
// with the attribute path fix in ALPSCore c7c64ca0. Changes: the fixture uses
// ALPS's TemporaryDirectory, and the AlpsSpelling* cases check that ALPS's own
// "path/@name" spelling keeps naming the same object as ALPSCore's "path@name";
// Delete is new, because ALPSCore does not implement delete_attribute.
// See src/alps/hdf5/tests/README.md.

#include <alps/hdf5/archive.hpp>
#include <alps/testing/temporary_directory.hpp>

#include <string>

#include <gtest/gtest.h>

class TestHDF5Attributes : public ::testing::Test {
  public:
    alps::testing::TemporaryDirectory temporary_;
    std::string fname_ = (temporary_.path() / "attr.h5").string();
};

namespace ah5=alps::hdf5;

TEST_F(TestHDF5Attributes, Data) {
    {
        ah5::archive ar(fname_,"w");
        ar["/int"] << int(123);
        ar["/int@attr"]="My data attribute";
    }
    ah5::archive ar(fname_,"r");
    ASSERT_TRUE(ar.is_attribute("/int@attr"));
    std::string attr;
    ar["/int@attr"] >> attr;
    EXPECT_EQ("My data attribute", attr);
}

TEST_F(TestHDF5Attributes, Group) {
    {
        ah5::archive ar(fname_,"w");
        ar["/group/int"]=int(123);
        ar["/group@attr"]="My group attribute";
    }
    ah5::archive ar(fname_,"r");
    ASSERT_TRUE(ar.is_attribute("/group@attr"));
    std::string attr;
    ar["/group@attr"] >> attr;
    EXPECT_EQ("My group attribute", attr);
}

TEST_F(TestHDF5Attributes, Root) {
    {
        ah5::archive ar(fname_,"w");
        ar["/group/int"]=int(123);
        ar["/@attr"]="My root attribute";
    }
    ah5::archive ar(fname_,"r");
    ASSERT_TRUE(ar.is_attribute("/@attr"));
    std::string attr;
    ar["/@attr"] >> attr;
    EXPECT_EQ("My root attribute", attr);
}

TEST_F(TestHDF5Attributes, Current) {
    {
        ah5::archive ar(fname_,"w");
        ar["/group/int2"]=int(321);
        ar.set_context("/group");
        ar["@attr"]="My group context attribute";
    }
    ah5::archive ar(fname_,"r");
    ar.set_context("/group");
    ASSERT_TRUE(ar.is_attribute("@attr"));
    std::string attr;
    ar["@attr"] >> attr;
    EXPECT_EQ("My group context attribute", attr);
}

TEST_F(TestHDF5Attributes, CurrentRoot) {
    {
        ah5::archive ar(fname_,"w");
        ar["/group/int"]=int(123);
        ar["@attr"]="My root context attribute";
    }
    ah5::archive ar(fname_,"r");
    ASSERT_TRUE(ar.is_attribute("/@attr"));
    std::string attr;
    ar["@attr"] >> attr;
    EXPECT_EQ("My root context attribute", attr);
}

TEST_F(TestHDF5Attributes, AlpsSpellingData) {
    {
        ah5::archive ar(fname_,"w");
        ar["/int"] << int(123);
        ar["/int/@attr"]="My data attribute";
    }
    ah5::archive ar(fname_,"r");
    ASSERT_TRUE(ar.is_attribute("/int@attr"));
    std::string attr;
    ar["/int/@attr"] >> attr;
    EXPECT_EQ("My data attribute", attr);
}

TEST_F(TestHDF5Attributes, AlpsSpellingGroup) {
    {
        ah5::archive ar(fname_,"w");
        ar["/group/int"]=int(123);
        ar["/group@attr"]="My group attribute";
    }
    ah5::archive ar(fname_,"r");
    ASSERT_TRUE(ar.is_attribute("/group/@attr"));
    std::string attr;
    ar["/group/@attr"] >> attr;
    EXPECT_EQ("My group attribute", attr);
    ar.set_context("/group");
    EXPECT_FALSE(ar.is_attribute("int@missing"));
    EXPECT_TRUE(ar.is_attribute("@attr"));
}

// ALPSCore's delete_attribute only throws "Not implemented"; ALPS's resolves
// the owner as the other attribute operations do, in both spellings.
TEST_F(TestHDF5Attributes, Delete) {
    ah5::archive ar(fname_,"w");
    ar["/group/int"]=int(123);
    ar["/grou/int"]=int(321);
    ar["/grou@attr"]="Sibling attribute";
    ar["/group@attr"]="ALPSCore spelling";
    ar.delete_attribute("/group@attr");
    EXPECT_FALSE(ar.is_attribute("/group@attr"));
    EXPECT_TRUE(ar.is_attribute("/grou@attr"));
    ar["/group/@attr"]="ALPS spelling";
    ar.delete_attribute("/group/@attr");
    EXPECT_FALSE(ar.is_attribute("/group@attr"));
    ar["/@attr"]="Root attribute";
    ar.delete_attribute("/@attr");
    EXPECT_FALSE(ar.is_attribute("/@attr"));
    EXPECT_TRUE(ar.is_attribute("/grou@attr"));
}
