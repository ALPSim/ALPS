// Copyright (C) 2006-2010 Adrian Feiguin
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT
#ifndef ALPS_DMTK_FILELIST_H
#define ALPS_DMTK_FILELIST_H

#include <alps/utility/temporary_filename.hpp>
#include <boost/filesystem.hpp>
#include <boost/throw_exception.hpp>
#include <iostream>
#include <map>
#include <string>
#include <stdexcept>

namespace dmtk {

// Own only the files allocated by this registry, never the user's scratch
// directory or files belonging to another calculation.
class FileList {
  public:
    FileList() : temp_dir(boost::filesystem::current_path()) {}
    explicit FileList(const char* dir) { set_temp_dir(dir); }
    FileList(const FileList&) = delete;
    FileList& operator=(const FileList&) = delete;
    ~FileList() { cleanup(); }

    void set_temp_dir(const char* dir)
    {
      // Store absolute paths so a later working-directory change cannot make
      // cleanup remove a different file or lose track of the original one.
      boost::filesystem::path path = boost::filesystem::absolute(dir);
      if (!boost::filesystem::is_directory(path))
        boost::throw_exception(std::runtime_error("ALPS DMRG temporary directory does not exist: " + path.string()));
      temp_dir = path;
      std::cout << "ALPS DMRG temporary files will be written to " << temp_dir.string() << std::endl;
    }

    const char* get_filename(const char* input)
    {
      auto old = _tmp_filenames.find(input);
      if (old != _tmp_filenames.end()) return old->second.c_str();
      std::string name(input);
      name = name.substr(0, name.find_first_of('.'));
      std::string filename = alps::temporary_filename((temp_dir / name).string());
      auto entry = _tmp_filenames.emplace(input, filename);
      std::cout << "Creating temp file " << filename << std::endl;
      return entry.first->second.c_str();
    }

    void cleanup() noexcept
    {
      for (const auto& entry : _tmp_filenames) {
        boost::system::error_code error;
        boost::filesystem::remove(entry.second, error);
        if (error)
          std::cerr << "ALPS DMRG could not remove temporary file " << entry.second
                    << ": " << error.message() << std::endl;
      }
      _tmp_filenames.clear();
    }

  private:
    std::map<std::string, std::string> _tmp_filenames;
    boost::filesystem::path temp_dir;
};

} // namespace dmtk
#endif
