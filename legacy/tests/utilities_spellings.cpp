// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// Code written against ALPS 3.0.0 keeps compiling through the legacy layer:
// each deprecated include spelling forwards to its canonical header, and the
// renamed alps::ngs names still refer to the canonical entities.
#include <gtest/gtest.h>

#include <alps/ngs/cast.hpp>
#include <alps/ngs/detail/remove_cvr.hpp>
#include <alps/ngs/detail/type_wrapper.hpp>
#include <alps/ngs/short_print.hpp>
#include <alps/ngs/signal.hpp>
#include <alps/ngs/stacktrace.hpp>
#include <alps/ngs/stringify.hpp>
#include <alps/utility/citations.hpp>
#include <alps/utility/copyright.hpp>
#include <alps/utility/data.hpp>
#include <alps/utility/encode.hpp>
#include <alps/utility/make_copy.hpp>
#include <alps/utility/os.hpp>
#include <alps/utility/resize.hpp>
#include <alps/utility/size.hpp>
#include <alps/utility/temporary_filename.hpp>
#include <alps/utility/vectorio.hpp>
#include <alps/utility/vmusage.hpp>

#include <cstring>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

// The deprecated names are the subject of this test.
#if defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

TEST(LegacyUtilitiesSpellings, NgsHeadersForwardToUtilities)
{
    EXPECT_EQ(alps::cast<int>(std::string("42")), 42);
    EXPECT_TRUE((std::is_same<alps::detail::remove_cvr<int const &>::type, int>::value));
    EXPECT_TRUE((std::is_same<alps::detail::type_wrapper<double>::type, double>::value));
    std::ostringstream out;
    out << alps::short_print(std::vector<int>{1, 2});
    EXPECT_EQ(out.str(), "[1,2]");
}

TEST(LegacyUtilitiesSpellings, UtilityHeadersForwardToUtilities)
{
    EXPECT_EQ(alps::make_copy(3), 3);
    EXPECT_EQ(alps::size(std::vector<int>(4)), 4u);
    EXPECT_EQ(alps::size(1.5), 1u);
    EXPECT_FALSE(alps::temporary_filename("legacy_spelling_").empty());
    EXPECT_FALSE(alps::citation_text().empty());
    EXPECT_FALSE(alps::version().empty());
    std::vector<int> storage{5, 6};
    EXPECT_EQ(alps::data(storage), storage.data());
    EXPECT_EQ(alps::hdf5_name_decode(alps::hdf5_name_encode("a/b")), "a/b");
    EXPECT_FALSE(alps::temp_directory_path().empty());
    std::vector<double> resized;
    alps::resize_same_as(resized, storage);
    EXPECT_EQ(resized.size(), 2u);
    std::vector<int> parsed;
    alps::read_vector_resize(std::string("1 2 3"), parsed);
    EXPECT_EQ(parsed, (std::vector<int>{1, 2, 3}));
    (void)alps::vmusage();
}

TEST(LegacyUtilitiesSpellings, NgsNamesReferToCanonicalEntities)
{
    EXPECT_TRUE((std::is_same<alps::ngs::signal, alps::signal>::value));
    EXPECT_EQ(std::strcmp(ALPS_NGS_STRINGIFY(alps), "alps"), 0);
    EXPECT_EQ(std::strcmp(ALPS_NGS_STRINGIFY(alps), ALPS_STRINGIFY(alps)), 0);
    // Both names reach the same stack walker; frames differ by the wrapper.
    std::string trace = alps::ngs::stacktrace();
    (void)trace;
    alps::ngs::signal::slot(0);
    alps::signal handler;
    EXPECT_FALSE(handler.empty());
    EXPECT_EQ(handler.top(), 0);
    handler.pop();
}
