// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// ALPS reads a file written by ALPSCore 7146b9e1's archive. ALPSCore stores
// bool and signed char without ALPS's __alps_type__ tag, spells attributes
// "path@name", and writes boost::optional with its own flag attribute. The
// fixture and its writer are in alpscore/; see alpscore/README.md.
#include "alpscore/items.hpp"

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/boost_optional.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/map.hpp>
#include <alps/hdf5/multi_array.hpp>
#include <alps/hdf5/pair.hpp>
#include <alps/hdf5/valarray.hpp>
#include <alps/hdf5/vector.hpp>

#include <boost/multi_array.hpp>
#include <boost/optional.hpp>

#include <string>

#include <gtest/gtest.h>

namespace {
    std::string const fixture = ALPS_TEST_SOURCE_DIR "/alpscore/archive_7146b9e1.h5";

    template <class T> void expect_stored(char const * name, T const & expected) {
        alps::hdf5::archive ar(fixture, "r");
        T value;
        ar[std::string("/data/") + name] >> value;
        EXPECT_TRUE(alps_fixture::same(value, expected)) << name;
    }
}

#define ALPS_FIXTURE_TEST(name, type, value) \
    TEST(Hdf5AlpscoreFixture, name) { typedef ALPS_FIXTURE_ARG type T; expect_stored<T>(#name, T(value)); }
ALPS_FIXTURE_ITEMS(ALPS_FIXTURE_TEST)
#undef ALPS_FIXTURE_TEST

TEST(Hdf5AlpscoreFixture, multi_array) {
    alps::hdf5::archive ar(fixture, "r");
    boost::multi_array<double, 2> stored;
    ar["/data/multi_array"] >> stored;
    ASSERT_EQ(stored.shape()[0], 2u);
    ASSERT_EQ(stored.shape()[1], 3u);
    for (std::size_t i = 0; i != stored.num_elements(); ++i)
        EXPECT_EQ(stored.data()[i], alps_fixture::multi_array_element(i)) << i;
}

TEST(Hdf5AlpscoreFixture, optional) {
    alps::hdf5::archive ar(fixture, "r");
    boost::optional<int> full, empty(5);
    ar["/data/optional_full"] >> full;
    ar["/data/optional_empty"] >> empty;
    ASSERT_TRUE(full);
    EXPECT_EQ(*full, 17);
    EXPECT_FALSE(empty);
}

TEST(Hdf5AlpscoreFixture, attributes_in_both_spellings) {
    alps::hdf5::archive ar(fixture, "r");
    std::string unit;
    int version = 0;
    ar["/data/int_scalar@unit"] >> unit;
    ar["/data@version"] >> version;
    EXPECT_EQ(unit, "meV");
    EXPECT_EQ(version, 3);
    unit.clear();
    version = 0;
    ar["/data/int_scalar/@unit"] >> unit;
    ar["/data/@version"] >> version;
    EXPECT_EQ(unit, "meV");
    EXPECT_EQ(version, 3);
}

TEST(Hdf5AlpscoreFixture, untagged_signed_bytes) {
    alps::hdf5::archive ar(fixture, "r");
    EXPECT_FALSE(ar.is_attribute("/data/bool_true/@__alps_type__"));
    EXPECT_FALSE(ar.is_attribute("/data/int8_scalar/@__alps_type__"));
}
