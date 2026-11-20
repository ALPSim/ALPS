#include <gtest/gtest.h>
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT
#include "checkpoint_cache.h"
#include <filesystem>
#include <fstream>
#include <iterator>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <system_error>

namespace fs = std::filesystem;

namespace {
std::string random_name(std::string name)
{
  std::mt19937 engine(std::random_device{}());
  for (int i = 0; i < 8; ++i)
    name += "0123456789abcdef"[std::uniform_int_distribution<int>(0, 15)(engine)];
  return name;
}

std::string slurp(const fs::path& p)
{
  std::ifstream in(p.string().c_str());
  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}
}

// DMTK reaches the checkpoint files only through the FileList they are attached to.
TEST(DmrgScratchFiles, CheckpointEntryCommits)
{
  const fs::path entry = fs::temp_directory_path() / random_name("alps-dmrg-entry-");
  fs::create_directories(entry);
  struct cleanup {
    fs::path directory;
    ~cleanup() { std::error_code ignored; fs::remove_all(directory, ignored); }
  } guard{entry};
  const fs::path block = entry / "block_ALPS_2_l.dat";
  {
    dmtk::FileList files(entry.string().c_str());
    dmrg_checkpoint::Files store;
    ASSERT_EQ(store.open(entry), "") << "empty entry has a record";
    store.attach(files);
    std::ofstream(files.get_write_filename("block_ALPS_2_l.dat")) << "first";
    ASSERT_EQ(slurp(files.get_read_filename("block_ALPS_2_l.dat")), "first") << "reads miss staged writes";
    ASSERT_FALSE(fs::exists(block)) << "write visible before commit";
    ASSERT_EQ(fs::path(files.get_filename("system.dat")), entry / "system.dat") << "file not kept by name";
    store.commit("record 1\n");
    ASSERT_EQ(slurp(block), "first") << "commit did not promote the staged file";

    // A rerun after a crash before the next commit sees the committed state.
    std::ofstream(files.get_write_filename("block_ALPS_2_l.dat")) << "second";
    dmrg_checkpoint::Files rerun;
    ASSERT_EQ(rerun.open(entry), "record 1\n");
    ASSERT_EQ(slurp(block), "first") << "uncommitted write replaced the committed file";
    for (fs::directory_iterator it(entry), end; it != end; ++it)
      ASSERT_EQ(it->path().filename().string().find(".new."), std::string::npos) << "uncommitted file kept";
    rerun.close();
    files.cleanup();
    ASSERT_EQ(slurp(block), "first") << "cleanup removed a committed checkpoint file";

    // cleanup() detached the store: FileList is back to temporary files.
    std::string temp = files.get_write_filename("block_ALPS_2_l.dat");
    ASSERT_EQ(temp.find("block_ALPS_2_l.dat"), std::string::npos) << "store still attached after cleanup";
    files.cleanup();
  }

  // A crash after the record was written but before the renames finished:
  // the files it lists are promoted.
  {
    std::ofstream(fs::path(entry / "progress").string().c_str())
      << "ALPS DMRG progress 1\ngeneration 7\nstaged 1\nblock_ALPS_2_l.dat\n"
      << "files 1\nblock_ALPS_2_l.dat 7\nrecord 7\n";
    std::ofstream((entry / "block_ALPS_2_l.dat.new.7").string().c_str()) << "seventh";
    std::ofstream((entry / "block_ALPS_2_l.dat.new.8").string().c_str()) << "uncommitted";
    dmrg_checkpoint::Files store;
    ASSERT_EQ(store.open(entry), "record 7\n");
    ASSERT_EQ(slurp(block), "seventh") << "interrupted commit was not completed";
    ASSERT_FALSE(fs::exists(entry / "block_ALPS_2_l.dat.new.8")) << "uncommitted file kept";
    ASSERT_NE(std::string(store.get_write_filename("block_ALPS_2_l.dat")).find(".new.8"), std::string::npos)
      << "generation did not continue after the recovered commit";
  }

  // A committed file that is missing or truncated is an error: DMTK would
  // otherwise go on with an empty block.
  {
    dmrg_checkpoint::Files store;
    store.open(entry);
    std::ofstream(store.get_write_filename("rho_ALPS_2_l.dat")) << "rho";
    store.commit("record 9\n");
    store.close();
    dmrg_checkpoint::Files rerun;
    ASSERT_EQ(rerun.open(entry), "record 9\n") << "manifest not accepted";
    rerun.close();
    std::ofstream((entry / "block_ALPS_2_l.dat").string().c_str()) << "short";
    dmrg_checkpoint::Files truncated;
    ASSERT_THROW(truncated.open(entry), std::runtime_error) << "truncated file accepted";
    fs::remove(entry / "rho_ALPS_2_l.dat");
    std::ofstream((entry / "block_ALPS_2_l.dat").string().c_str()) << "seventh";
    dmrg_checkpoint::Files missing;
    ASSERT_THROW(missing.open(entry), std::runtime_error) << "missing file accepted";
    fs::remove(entry / "progress");
  }

  // Failed writes and unreadable records are errors, not silent corruption.
  {
    dmtk::FileList files(entry.string().c_str());
    dmrg_checkpoint::Files store;
    store.open(entry);
    store.attach(files);
    files.write_failed("rho_ALPS_2_l.dat");
    ASSERT_THROW(store.commit("record\n"), std::runtime_error);
  }
  std::ofstream(fs::path(entry / "progress").string().c_str()) << "garbage\n";
  dmrg_checkpoint::Files store;
  ASSERT_THROW(store.open(entry), std::runtime_error);
}

// A store destroyed while attached detaches itself and removes its staged files.
TEST(DmrgScratchFiles, CheckpointStoreOutlivedByFileList)
{
  const fs::path entry = fs::temp_directory_path() / random_name("alps-dmrg-entry-");
  fs::create_directories(entry);
  struct cleanup {
    fs::path directory;
    ~cleanup() { std::error_code ignored; fs::remove_all(directory, ignored); }
  } guard{entry};
  dmtk::FileList files(entry.string().c_str());
  std::string staged;
  {
    dmrg_checkpoint::Files store;
    store.open(entry);
    store.attach(files);
    staged = files.get_write_filename("block_ALPS_2_l.dat");
    std::ofstream(staged.c_str()) << "uncommitted";
  }
  ASSERT_FALSE(fs::exists(staged)) << "staged file kept after the store was destroyed";
  std::string temp = files.get_write_filename("block_ALPS_2_l.dat");
  ASSERT_EQ(temp.find("block_ALPS_2_l.dat"), std::string::npos) << "FileList still uses a destroyed store";
}
