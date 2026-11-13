// SPDX-License-Identifier: MIT
#include <alps/parapack/option.h>
#include <gtest/gtest.h>

#include <initializer_list>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

alps::parapack::option parse(std::initializer_list<char const*> args) {
  std::vector<std::string> storage{"parapack"};
  for (auto const* arg : args) storage.emplace_back(arg);
  std::vector<char*> argv;
  for (auto& arg : storage) argv.push_back(&arg[0]);
  argv.push_back(nullptr);
  // A rejected command line prints its reason; keep it out of the test log,
  // and restore std::cerr even if the constructor throws.
  struct silence_cerr {
    std::ostringstream sink;
    std::streambuf* previous = std::cerr.rdbuf(sink.rdbuf());
    ~silence_cerr() { std::cerr.rdbuf(previous); }
  } quiet;
  return alps::parapack::option(static_cast<int>(storage.size()), argv.data());
}

}

// pyalps runs parapack applications with --evaluate, --write-xml and --mpi.
TEST(ParapackOption, EvaluateAndWriteXml) {
  auto opt = parse({"--evaluate", "--write-xml", "job.in.xml"});
  ASSERT_TRUE(opt.valid);
  EXPECT_TRUE(opt.evaluate_only);
  EXPECT_TRUE(opt.write_xml);
  EXPECT_EQ(opt.jobfiles, std::vector<std::string>{"job.in.xml"});

  auto defaults = parse({"first.in.xml", "second.in.xml"});
  ASSERT_TRUE(defaults.valid);
  EXPECT_FALSE(defaults.evaluate_only || defaults.write_xml);
  EXPECT_EQ(defaults.jobfiles.size(), 2u);
}

TEST(ParapackOption, TotalThreads) {
  auto automatic = parse({"--total-threads", "auto", "job.in.xml"});
  ASSERT_TRUE(automatic.valid);
  EXPECT_TRUE(automatic.auto_total_threads && !automatic.default_total_threads);

  auto counted = parse({"-r", "4", "-p2", "job.in.xml"});
  ASSERT_TRUE(counted.valid);
  EXPECT_FALSE(counted.auto_total_threads || counted.default_total_threads);
  EXPECT_EQ(counted.num_total_threads, 4);
  EXPECT_EQ(counted.threads_per_clone, 2);

  EXPECT_TRUE(parse({"job.in.xml"}).default_total_threads);
}

TEST(ParapackOption, RejectsBadCommandLines) {
  EXPECT_TRUE(parse({"--dump-format", "xdr", "job.in.xml"}).valid);
  EXPECT_FALSE(parse({"--dump-format", "text", "job.in.xml"}).valid);
  EXPECT_FALSE(parse({"--time-limit", "soon", "job.in.xml"}).valid);
  EXPECT_FALSE(parse({"--no-such-option", "job.in.xml"}).valid);
}
