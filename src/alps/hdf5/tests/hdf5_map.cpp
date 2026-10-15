// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// std::map storage: one child per key below the map's group.
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/map.hpp>
#include <alps/testing/temporary_directory.hpp>

#include <map>
#include <string>

#include <gtest/gtest.h>

namespace {
    class Hdf5Map : public ::testing::Test {
      protected:
        alps::testing::TemporaryDirectory temporary_;
        std::string file_ = (temporary_.path() / "map.h5").string();
    };
}

TEST_F(Hdf5Map, RoundTripsKeysThatNeedEncoding) {
    std::map<std::string, double> values{{"a", 1.0}, {"b/c", 2.0}, {"d e", 3.0}}, restored;
    {
        alps::hdf5::archive ar(file_, "w");
        ar["/values"] << values;
    }
    alps::hdf5::archive ar(file_, "r");
    ar["/values"] >> restored;
    EXPECT_EQ(restored, values);
}

TEST_F(Hdf5Map, RoundTripsRelativeToTheContext) {
    std::map<int, double> values{{1, 0.5}, {2, 1.5}}, restored;
    {
        alps::hdf5::archive ar(file_, "w");
        ar["/simulation/marker"] << 1;
        ar.set_context("/simulation");
        ar["values"] << values;
    }
    alps::hdf5::archive ar(file_, "r");
    ar.set_context("/simulation");
    ar["values"] >> restored;
    EXPECT_EQ(restored, values);
}

// Loading through an empty path, as a class's load(archive&) does inside its
// own context, reads the entries below that context, not below the root.
TEST_F(Hdf5Map, LoadsThroughAnEmptyPathInsideTheContext) {
    std::map<std::string, double> values{{"a", 1.0}, {"b", 2.0}}, restored;
    {
        alps::hdf5::archive ar(file_, "w");
        ar["/simulation/values"] << values;
    }
    alps::hdf5::archive ar(file_, "r");
    ar.set_context("/simulation/values");
    ar[""] >> restored;
    EXPECT_EQ(restored, values);
}

// An empty map writes its group (ALPSCore 99e8fa59), so it reads back as
// empty instead of a missing path.
TEST_F(Hdf5Map, EmptyMapRoundTrips) {
    std::map<std::string, double> empty, restored{{"stale", 1.0}};
    {
        alps::hdf5::archive ar(file_, "w");
        ar["/empty"] << empty;
    }
    alps::hdf5::archive ar(file_, "r");
    ASSERT_TRUE(ar.is_group("/empty"));
    EXPECT_TRUE(ar.list_children("/empty").empty());
    restored.clear();
    ar["/empty"] >> restored;
    EXPECT_TRUE(restored.empty());
}

TEST_F(Hdf5Map, EmptyMapReplacesAStoredMap) {
    std::map<std::string, double> values{{"a", 1.0}}, empty, restored;
    {
        alps::hdf5::archive ar(file_, "w");
        ar["/values"] << values;
        ar["/values"] << empty;
    }
    alps::hdf5::archive ar(file_, "r");
    ar["/values"] >> restored;
    EXPECT_TRUE(restored.empty());
}
