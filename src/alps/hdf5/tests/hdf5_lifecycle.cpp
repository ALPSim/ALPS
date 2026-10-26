// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// Archive handles: a default-constructed archive opened later, reopening
// after close (ALPSCore 38e49c73), and no copy assignment (ALPSCore 892bff78),
// which used to share a file context without counting the reference.
#include <alps/hdf5/archive.hpp>
#include <alps/testing/temporary_directory.hpp>

#include <string>
#include <type_traits>

#include <gtest/gtest.h>

namespace {
    class Hdf5Lifecycle : public ::testing::Test {
      protected:
        alps::testing::TemporaryDirectory temporary_;
        std::string first_ = (temporary_.path() / "first.h5").string();
        std::string second_ = (temporary_.path() / "second.h5").string();
    };
}

TEST(Hdf5Archive, CopiesButDoesNotAssign) {
    EXPECT_TRUE(std::is_copy_constructible<alps::hdf5::archive>::value);
    EXPECT_FALSE(std::is_copy_assignable<alps::hdf5::archive>::value);
    EXPECT_TRUE(std::is_default_constructible<alps::hdf5::archive>::value);
}

TEST_F(Hdf5Lifecycle, DefaultConstructedArchiveIsClosed) {
    alps::hdf5::archive ar;
    EXPECT_FALSE(ar.is_open());
    EXPECT_THROW(ar.get_filename(), alps::hdf5::archive_closed);
    EXPECT_THROW(ar.close(), alps::hdf5::archive_closed);
}

TEST_F(Hdf5Lifecycle, OpensAfterDefaultConstruction) {
    {
        alps::hdf5::archive ar;
        ar.open(first_, "w");
        ASSERT_TRUE(ar.is_open());
        EXPECT_EQ(ar.get_context(), "/");
        ar["/value"] << 7;
    }
    alps::hdf5::archive ar;
    ar.open(first_);
    int value = 0;
    ar["/value"] >> value;
    EXPECT_EQ(value, 7);
}

TEST_F(Hdf5Lifecycle, ReopensAnotherFileAfterClose) {
    alps::hdf5::archive ar(first_, "w");
    ar["/value"] << 1;
    ar.set_context("/somewhere");
    ar.close();
    ar.open(second_, "w");
    EXPECT_EQ(ar.get_filename(), second_);
    EXPECT_EQ(ar.get_context(), "/");
    ar["/value"] << 2;
    ar.close();
    alps::hdf5::archive first(first_), second(second_);
    int a = 0, b = 0;
    first["/value"] >> a;
    second["/value"] >> b;
    EXPECT_EQ(a, 1);
    EXPECT_EQ(b, 2);
}

TEST_F(Hdf5Lifecycle, RejectsOpeningAnOpenArchive) {
    alps::hdf5::archive ar(first_, "w");
    EXPECT_THROW(ar.open(second_, "w"), alps::hdf5::archive_opened);
    EXPECT_EQ(ar.get_filename(), first_);
}

TEST_F(Hdf5Lifecycle, OpenValidatesItsMode) {
    alps::hdf5::archive ar;
    EXPECT_THROW(ar.open(first_, "wx"), alps::hdf5::wrong_mode);
    EXPECT_FALSE(ar.is_open());
    EXPECT_THROW(ar.open(first_, "r"), alps::hdf5::archive_not_found);
    EXPECT_FALSE(ar.is_open());
}
