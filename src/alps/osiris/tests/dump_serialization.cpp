/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <gtest/gtest.h>
#include <alps/osiris.h>
#include <alps/osiris/archivedump.h>
#include <alps/testing/temporary_directory.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <fstream>
#include <iterator>

namespace {
// The legacy fixture has eleven fields. The native extended record additionally
// exercises a positive 64-bit integer and a complex value in their original order.
template<bool Extended = false, class Output>
void write_record(Output& output) {
  output << bool(false) << int8_t(63) << uint8_t(201) << int16_t(-699)
         << uint16_t(43299) << int32_t(847229) << uint32_t(4294967295u)
         << int64_t(-1152921504606846976ll);
  if constexpr (Extended) output << int64_t(343434545665ll);
  output << uint64_t(18446744073709551614ull) << double(3.14159265358979323846)
         << std::string("test string");
  if constexpr (Extended) output << std::complex<double>(1, 2);
}

// Keep both reading APIs explicit: native get/conversion operators and the
// Boost archive adapter's free get function are distinct compatibility paths.
void expect_record(alps::IDump& input, bool extended = false) {
  EXPECT_EQ(input.get<bool>(), false);
  EXPECT_EQ(input.get<int8_t>(), 63);
  EXPECT_EQ(input.get<uint8_t>(), 201);
  EXPECT_EQ(input.get<int16_t>(), -699);
  EXPECT_EQ(input.get<uint16_t>(), 43299);
  EXPECT_EQ(static_cast<int32_t>(input), 847229);
  EXPECT_EQ(static_cast<uint32_t>(input), 4294967295u);
  EXPECT_EQ(static_cast<int64_t>(input), -1152921504606846976ll);
  if (extended) EXPECT_EQ(static_cast<int64_t>(input), 343434545665ll);
  EXPECT_EQ(static_cast<uint64_t>(input), 18446744073709551614ull);
  // Full double precision is part of the contract for both wire and text archives.
  EXPECT_EQ(static_cast<double>(input), 3.14159265358979323846);
  std::string text;
  input >> text;
  EXPECT_EQ(text, "test string");
  if (extended) {
    std::complex<double> number;
    input >> number;
    EXPECT_EQ(number, std::complex<double>(1, 2));
  }
}

void expect_record(alps::idump_archive& input) {
  EXPECT_EQ(alps::get<bool>(input), false);
  EXPECT_EQ(alps::get<int8_t>(input), 63);
  EXPECT_EQ(alps::get<uint8_t>(input), 201);
  EXPECT_EQ(alps::get<int16_t>(input), -699);
  EXPECT_EQ(alps::get<uint16_t>(input), 43299);
  EXPECT_EQ(alps::get<int32_t>(input), 847229);
  EXPECT_EQ(alps::get<uint32_t>(input), 4294967295u);
  EXPECT_EQ(alps::get<int64_t>(input), -1152921504606846976ll);
  EXPECT_EQ(alps::get<uint64_t>(input), 18446744073709551614ull);
  EXPECT_EQ(alps::get<double>(input), 3.14159265358979323846);
  std::string text;
  input >> text;
  EXPECT_EQ(text, "test string");
}
} // namespace

TEST(OsirisSerialization, PreservesScalarAndComplexValues) {
  alps::testing::TemporaryDirectory directory;
  const auto path = directory.path() / "record.dump";
  {
    alps::OXDRFileDump output(boost::filesystem::path(path.string()));
    write_record<true>(output);
  }
  alps::IXDRFileDump input(boost::filesystem::path(path.string()));
  expect_record(input, true);
}

TEST(OsirisSerialization, ReadsHistoricalXdrRecord) {
  // These historical readers are independent of the current writer.
  alps::IXDRFileDump input(boost::filesystem::path(ALPS_TEST_SOURCE_DIR "/xdrdump2.dump"));
  expect_record(input);
}

TEST(OsirisSerialization, ArchiveAdapterReadsHistoricalXdrRecord) {
  alps::IXDRFileDump input(boost::filesystem::path(ALPS_TEST_SOURCE_DIR "/xdrdump2.dump"));
  alps::idump_archive archive(input);
  expect_record(archive);
}

TEST(OsirisSerialization, TextArchivePreservesScalarValues) {
  alps::testing::TemporaryDirectory directory;
  const auto path = directory.path() / "record.dump";
  {
    std::ofstream stream(path);
    ASSERT_TRUE(stream.good());
    boost::archive::text_oarchive archive(stream);
    alps::archive_odump<boost::archive::text_oarchive> output(archive);
    write_record(output);
  }
  std::ifstream stream(path);
  ASSERT_TRUE(stream.good());
  boost::archive::text_iarchive archive(stream);
  alps::archive_idump<boost::archive::text_iarchive> input(archive);
  expect_record(input);
}

TEST(OsirisSerialization, ArchiveAdapterReadsNativeDump) {
  alps::testing::TemporaryDirectory directory;
  const auto path = directory.path() / "record.dump";
  {
    alps::OXDRFileDump output(boost::filesystem::path(path.string()));
    write_record(output);
  }
  alps::IXDRFileDump input(boost::filesystem::path(path.string()));
  alps::idump_archive archive(input);
  expect_record(archive);
}

TEST(OsirisSerialization, NativeDumpReadsArchiveAdapter) {
  alps::testing::TemporaryDirectory directory;
  const auto path = directory.path() / "record.dump";
  {
    alps::OXDRFileDump dump(boost::filesystem::path(path.string()));
    alps::odump_archive output(dump);
    write_record(output);
  }
  alps::IXDRFileDump input(boost::filesystem::path(path.string()));
  expect_record(input);

  // A matching writer/reader bug must not silently change the checkpoint format.
  std::ifstream expected(ALPS_TEST_SOURCE_DIR "/xdrdump2.dump", std::ios::binary);
  std::ifstream actual(path, std::ios::binary);
  ASSERT_TRUE(expected.good());
  ASSERT_TRUE(actual.good());
  const std::string expected_bytes{std::istreambuf_iterator<char>(expected), {}};
  const std::string actual_bytes{std::istreambuf_iterator<char>(actual), {}};
  EXPECT_EQ(actual_bytes, expected_bytes);
}
