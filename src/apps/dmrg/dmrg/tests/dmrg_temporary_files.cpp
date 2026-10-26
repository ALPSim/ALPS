#include <gtest/gtest.h>
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT
#include "dmtk/filelist.h"
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

namespace fs = boost::filesystem;

TEST(DmrgScratchFiles, LifetimesAndCleanup)
{
  const fs::path root = fs::temp_directory_path() / fs::unique_path("alps-dmrg-%%%%-%%%%");
  fs::create_directories(root / "first");
  fs::create_directories(root / "second");
  const fs::path original_directory = fs::current_path();
  struct cleanup {
    fs::path previous, directory;
    ~cleanup() {
      boost::system::error_code ignored;
      fs::current_path(previous, ignored);
      fs::remove_all(directory, ignored);
    }
  } guard{original_directory, root};
  try {
    const fs::path unrelated = root / "first" / "block_ALPS_unrelated";
    std::ofstream(unrelated.string()) << "another run";
    {
      dmtk::FileList files((root / "first").string().c_str());
      const char* first = files.get_filename("block_ALPS_1.dat");
      const std::string first_name(first);
      const std::string second = files.get_filename("rho_ALPS_1.dat");
      ASSERT_TRUE((first_name == first)) << "a later allocation invalidated a filename pointer";
      ASSERT_TRUE((first_name == files.get_filename("block_ALPS_1.dat"))) << "filename was not reused";
      ASSERT_TRUE((fs::exists(first_name) && fs::exists(second))) << "scratch files missing while active";
      files.cleanup();
      ASSERT_TRUE((!fs::exists(first_name) && !fs::exists(second))) << "explicit cleanup leaked files";
      ASSERT_TRUE((fs::exists(unrelated))) << "cleanup removed an unrelated file";
      files.cleanup(); // Idempotent; also tolerate a file removed by a caller.
      files.set_temp_dir((root / "second").string().c_str());
      const std::string next = files.get_filename("block_ALPS_1.dat");
      ASSERT_TRUE((fs::path(next).parent_path() == root / "second")) << "next task reused the old directory";
      fs::remove(next);
      files.get_filename("gs_ALPS_1.dat");
    }
    ASSERT_TRUE((fs::is_empty(root / "second"))) << "destructor leaked files";
    try {
      fs::current_path(root);
      dmtk::FileList files("second");
      files.get_filename("system.dat");
      fs::current_path(original_directory);
      throw std::runtime_error("simulate calculation failure");
    } catch (const std::runtime_error&) {}
    ASSERT_TRUE((fs::is_empty(root / "second"))) << "exception or directory change prevented cleanup";

    // Relative scratch paths must still work when their absolute form is
    // longer than the old temporary_filename buffer. Each path component
    // stays short enough for ordinary filesystem limits.
    fs::path long_directory = root / "long";
    while (long_directory.string().size() < 300)
      long_directory /= std::string(60, 'd');
    fs::create_directories(long_directory / "scratch");
    fs::current_path(long_directory);
    const fs::path scratch = fs::current_path() / "scratch";
    const fs::path sentinel = scratch / "unrelated";
    std::ofstream(sentinel.string()) << "keep me";
    {
      dmtk::FileList files("scratch");
      const fs::path block = files.get_filename("block_ALPS_1.dat");
      const fs::path rho = files.get_filename("rho_ALPS_1.dat");
      ASSERT_TRUE((block.parent_path() == scratch && rho.parent_path() == scratch)) << "long scratch path was truncated or escaped its directory";
      ASSERT_TRUE((block != rho && fs::exists(block) && fs::exists(rho))) << "long scratch paths did not reserve distinct files";
      fs::current_path(original_directory);
    }
    ASSERT_TRUE((std::distance(fs::directory_iterator(scratch), fs::directory_iterator()) == 1
            && fs::exists(sentinel))) << "long-path cleanup leaked files or removed an unrelated file";

    // The first available descriptor must remain available after repeated calls.
    int before = open("/dev/null", O_RDONLY);
    ASSERT_TRUE((before >= 0)) << "could not open descriptor probe";
    close(before);
    std::vector<std::string> names;
    for (int i = 0; i < 64; ++i)
      names.push_back(alps::temporary_filename((root / "descriptor").string()));
    int after = open("/dev/null", O_RDONLY);
    ASSERT_TRUE((after >= 0)) << "temporary_filename exhausted descriptors";
    close(after);
    ASSERT_TRUE((after == before)) << "temporary_filename leaked a file descriptor";
    for (const auto& name : names) fs::remove(name);
    fs::remove(unrelated);
    fs::remove_all(root);
  } catch (...) {
    fs::current_path(original_directory);
    fs::remove_all(root);
    throw;
  }
}

namespace {
std::string slurp(const fs::path& p)
{
  std::ifstream in(p.string().c_str());
  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}
}

TEST(DmrgScratchFiles, CheckpointEntryCommits)
{
  const fs::path entry = fs::temp_directory_path() / fs::unique_path("alps-dmrg-entry-%%%%-%%%%");
  fs::create_directories(entry);
  struct cleanup {
    fs::path directory;
    ~cleanup() { boost::system::error_code ignored; fs::remove_all(directory, ignored); }
  } guard{entry};
  const fs::path block = entry / "block_ALPS_2_l.dat";
  {
    dmtk::FileList files;
    ASSERT_EQ(files.open_entry(entry.string().c_str()), "") << "empty entry has a record";
    std::ofstream(files.get_write_filename("block_ALPS_2_l.dat")) << "first";
    ASSERT_EQ(slurp(files.get_read_filename("block_ALPS_2_l.dat")), "first") << "reads miss staged writes";
    ASSERT_FALSE(fs::exists(block)) << "write visible before commit";
    files.commit("record 1\n");
    ASSERT_EQ(slurp(block), "first") << "commit did not promote the staged file";

    // A rerun after a crash before the next commit sees the committed state.
    std::ofstream(files.get_write_filename("block_ALPS_2_l.dat")) << "second";
    dmtk::FileList rerun;
    ASSERT_EQ(rerun.open_entry(entry.string().c_str()), "record 1\n");
    ASSERT_EQ(slurp(block), "first") << "uncommitted write replaced the committed file";
    for (fs::directory_iterator it(entry), end; it != end; ++it)
      ASSERT_EQ(it->path().filename().string().find(".new."), std::string::npos) << "uncommitted file kept";
    rerun.cleanup();
    files.cleanup();
    ASSERT_EQ(slurp(block), "first") << "cleanup removed a committed checkpoint file";
  }

  // A crash after the record was written but before the renames finished:
  // the files it lists are promoted.
  {
    std::ofstream(fs::path(entry / "progress").string().c_str())
      << "ALPS DMRG progress 1\ngeneration 7\nstaged 1\nblock_ALPS_2_l.dat\nrecord 7\n";
    std::ofstream((entry / "block_ALPS_2_l.dat.new.7").string().c_str()) << "seventh";
    std::ofstream((entry / "block_ALPS_2_l.dat.new.8").string().c_str()) << "uncommitted";
    dmtk::FileList files;
    ASSERT_EQ(files.open_entry(entry.string().c_str()), "record 7\n");
    ASSERT_EQ(slurp(block), "seventh") << "interrupted commit was not completed";
    ASSERT_FALSE(fs::exists(entry / "block_ALPS_2_l.dat.new.8")) << "uncommitted file kept";
    ASSERT_NE(std::string(files.get_write_filename("block_ALPS_2_l.dat")).find(".new.8"), std::string::npos)
      << "generation did not continue after the recovered commit";
  }

  // Failed writes and unreadable records are errors, not silent corruption.
  {
    dmtk::FileList files;
    files.open_entry(entry.string().c_str());
    files.write_failed("rho_ALPS_2_l.dat");
    ASSERT_THROW(files.commit("record\n"), std::runtime_error);
  }
  std::ofstream(fs::path(entry / "progress").string().c_str()) << "garbage\n";
  dmtk::FileList files;
  ASSERT_THROW(files.open_entry(entry.string().c_str()), std::runtime_error);
}
