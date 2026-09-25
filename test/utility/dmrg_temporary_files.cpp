// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT
#include "applications/dmrg/dmrg/dmtk/filelist.h"
#include <fstream>
#include <stdexcept>
#include <vector>
#ifndef BOOST_MSVC
#include <fcntl.h>
#include <unistd.h>
#endif

namespace fs = boost::filesystem;

void require(bool condition, const char* message)
{
  if (!condition) throw std::runtime_error(message);
}

int main()
{
  const fs::path root = fs::temp_directory_path() / fs::unique_path("alps-dmrg-%%%%-%%%%");
  fs::create_directories(root / "first");
  fs::create_directories(root / "second");
  const fs::path original_directory = fs::current_path();
  try {
    const fs::path unrelated = root / "first" / "block_ALPS_unrelated";
    std::ofstream(unrelated.string()) << "another run";
    {
      dmtk::FileList files((root / "first").string().c_str());
      const char* first = files.get_filename("block_ALPS_1.dat");
      const std::string first_name(first);
      const std::string second = files.get_filename("rho_ALPS_1.dat");
      require(first_name == first, "a later allocation invalidated a filename pointer");
      require(first_name == files.get_filename("block_ALPS_1.dat"), "filename was not reused");
      require(fs::exists(first_name) && fs::exists(second), "scratch files missing while active");
      files.cleanup();
      require(!fs::exists(first_name) && !fs::exists(second), "explicit cleanup leaked files");
      require(fs::exists(unrelated), "cleanup removed an unrelated file");
      files.cleanup(); // Idempotent; also tolerate a file removed by a caller.
      files.set_temp_dir((root / "second").string().c_str());
      const std::string next = files.get_filename("block_ALPS_1.dat");
      require(fs::path(next).parent_path() == root / "second", "next task reused the old directory");
      fs::remove(next);
      files.get_filename("gs_ALPS_1.dat");
    }
    require(fs::is_empty(root / "second"), "destructor leaked files");
    try {
      fs::current_path(root);
      dmtk::FileList files("second");
      files.get_filename("system.dat");
      fs::current_path(original_directory);
      throw std::runtime_error("simulate calculation failure");
    } catch (const std::runtime_error&) {}
    require(fs::is_empty(root / "second"), "exception or directory change prevented cleanup");

#ifndef BOOST_MSVC
    // The first available descriptor must remain available after repeated calls.
    int before = open("/dev/null", O_RDONLY);
    require(before >= 0, "could not open descriptor probe");
    close(before);
    std::vector<std::string> names;
    for (int i = 0; i < 64; ++i)
      names.push_back(alps::temporary_filename((root / "descriptor").string()));
    int after = open("/dev/null", O_RDONLY);
    require(after >= 0, "temporary_filename exhausted descriptors");
    close(after);
    require(after == before, "temporary_filename leaked a file descriptor");
    for (const auto& name : names) fs::remove(name);
#endif
    fs::remove(unrelated);
    fs::remove_all(root);
  } catch (...) {
    fs::current_path(original_directory);
    fs::remove_all(root);
    throw;
  }
}
