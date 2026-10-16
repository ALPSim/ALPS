/*
 * Copyright (C) 1998-2018 ALPS Collaboration. See COPYRIGHT.TXT
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

// Adapted from ALPSCore 7146b9e1 hdf5/test/hdf5_boost_optional.cpp, the test
// for the adapter from ALPSCore 7b5b7be5. Changes: ALPS's temporary directory;
// fixed values instead of ALPSCore's testing::datapoint helper;
// TYPED_TEST_SUITE for current GoogleTest; ReadIntoEmpty is new and covers
// the load fix in <alps/hdf5/boost_optional.hpp>.


#include <boost/optional/optional_io.hpp> // to make gtest happy with printing optionals

#include <alps/hdf5.hpp>
#include <alps/hdf5/boost_optional.hpp>

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <string>



namespace ah5=alps::hdf5;

template <typename T>
class Hdf5BoostOptionalTest : public ::testing::Test {
    alps::testing::TemporaryDirectory temporary_;
    std::string tmpfile_ = (temporary_.path() / "boost_optional.h5").string();

    // Two distinct values, standing in for ALPSCore's testing::datapoint.
    static T first() { return T(42); }
    static T second() { return T(-7); }

  public:
    typedef T value_type;

    void test_w_empty() {
        ah5::archive ar(tmpfile_, "w");
        boost::optional<value_type> empty;
        ar["empty"] << empty;
    }

    void test_rw_empty() {
        test_w_empty();
        ah5::archive ar(tmpfile_,"r");
        boost::optional<value_type> val=first();
        ar["empty"] >> val;
        EXPECT_EQ(boost::none, val);
    }

    void test_rw() {
        value_type v1=first();
        {
            ah5::archive ar(tmpfile_, "w");
            boost::optional<value_type> val=v1;
            ar["value"] << val;
        }
        ah5::archive ar(tmpfile_, "r");
        value_type v2=second();
        boost::optional<value_type> val=v2;
        ar["value"] >> val;
        ASSERT_NE(boost::none, val) << "Expected to be non-empty!";
        EXPECT_EQ(v1, *val);
    }

    void test_wr_type() { // reading the value of optional<T> from T
        value_type v1=first();
        {
            ah5::archive ar(tmpfile_, "w");
            ar["value"] << v1;
        }
        ah5::archive ar(tmpfile_, "r");
        value_type v2=second();
        boost::optional<value_type> val=v2;
        EXPECT_THROW(ar["value"] >> val, ah5::path_not_found);
        ASSERT_NE(boost::none, val) << "Expected to be non-empty!";
        EXPECT_EQ(v2, *val) << "Expected to preserve the value!";
    }

    void test_rd_type() { // reading the value of T from optional<T>
        value_type v1=first();
        {
            ah5::archive ar(tmpfile_, "w");
            boost::optional<value_type> val=v1;
            ar["value"] << val;
        }
        ah5::archive ar(tmpfile_, "r");
        value_type v2=second();
        ar["value"] >> v2;
        EXPECT_EQ(v1, v2);
    }

    void test_rw_into_empty() { // reading a stored value into an empty optional
        value_type v1=first();
        {
            ah5::archive ar(tmpfile_, "w");
            boost::optional<value_type> val=v1;
            ar["value"] << val;
        }
        ah5::archive ar(tmpfile_, "r");
        boost::optional<value_type> val;
        ar["value"] >> val;
        ASSERT_NE(boost::none, val) << "Expected to be non-empty!";
        EXPECT_EQ(v1, *val);
    }
};

typedef ::testing::Types<int> my_test_types;
TYPED_TEST_SUITE(Hdf5BoostOptionalTest, my_test_types);

TYPED_TEST(Hdf5BoostOptionalTest, WriteEmpty) { this->test_w_empty(); }
TYPED_TEST(Hdf5BoostOptionalTest, ReadWriteEmpty) { this->test_rw_empty(); }
TYPED_TEST(Hdf5BoostOptionalTest, ReadWrite) { this->test_rw(); }
TYPED_TEST(Hdf5BoostOptionalTest, WriteTypeReadOpt) { this->test_wr_type(); }
TYPED_TEST(Hdf5BoostOptionalTest, WriteOptReadType) { this->test_rd_type(); }
TYPED_TEST(Hdf5BoostOptionalTest, ReadIntoEmpty) { this->test_rw_into_empty(); }
