// Copyright (C) 2006-2010 Adrian Feiguin
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT
#ifndef ALPS_DMTK_FILELIST_H
#define ALPS_DMTK_FILELIST_H

#include <alps/utilities/temporary_filename.hpp>
#include <filesystem>
#include <boost/throw_exception.hpp>
#include <iostream>
#include <map>
#include <string>
#include <stdexcept>

namespace dmtk {

// Keeps the files under their own names elsewhere, e.g. the DMRG checkpoints
// (checkpoint_cache.h). While one is attached, FileList passes every file to it.
class FileStore {
  public:
    virtual const char* get_read_filename(const char* input) = 0;
    virtual const char* get_write_filename(const char* input) = 0;
    virtual void write_failed(const char* input) = 0;
    // Called when FileList detaches the store.
    virtual void cleanup() noexcept = 0;
  protected:
    ~FileStore() {}
};

// Own only the files allocated by this registry, never the user's scratch
// directory or files belonging to another calculation.
class FileList {
  public:
    FileList() : temp_dir(std::filesystem::current_path()) {}
    explicit FileList(const char* dir) { set_temp_dir(dir); }
    FileList(const FileList&) = delete;
    FileList& operator=(const FileList&) = delete;
    ~FileList() { cleanup(); }

    void set_temp_dir(const char* dir)
    {
      if (_store) cleanup();
      // Store absolute paths so a later working-directory change cannot make
      // cleanup remove a different file or lose track of the original one.
      std::filesystem::path path = std::filesystem::absolute(dir);
      if (!std::filesystem::is_directory(path))
        boost::throw_exception(std::runtime_error("ALPS DMRG temporary directory does not exist: " + path.string()));
      temp_dir = path;
      std::cout << "ALPS DMRG temporary files will be written to " << temp_dir.string() << std::endl;
    }

    // Attached until the next cleanup().
    void attach(FileStore* store) { cleanup(); _store = store; }

    const char* get_filename(const char* input)
    {
      if (_store) return _store->get_read_filename(input);
      auto old = _tmp_filenames.find(input);
      if (old != _tmp_filenames.end()) return old->second.c_str();
      std::string name(input);
      name = name.substr(0, name.find_first_of('.'));
      std::string filename = alps::temporary_filename((temp_dir / name).string());
      auto entry = _tmp_filenames.emplace(input, filename);
      std::cout << "Creating temp file " << filename << std::endl;
      return entry.first->second.c_str();
    }

    const char* get_read_filename(const char* input) { return get_filename(input); }
    const char* get_write_filename(const char* input)
      { return _store ? _store->get_write_filename(input) : get_filename(input); }
    void write_failed(const char* input) { if (_store) _store->write_failed(input); }

    void cleanup() noexcept
    {
      FileStore* store = _store;
      _store = nullptr;
      if (store) store->cleanup();
      for (const auto& entry : _tmp_filenames) {
        std::error_code error;
        std::filesystem::remove(entry.second, error);
        if (error)
          std::cerr << "ALPS DMRG could not remove temporary file " << entry.second
                    << ": " << error.message() << std::endl;
      }
      _tmp_filenames.clear();
    }

  private:
    std::map<std::string, std::string> _tmp_filenames;
    std::filesystem::path temp_dir;
    FileStore* _store = nullptr;
};

} // namespace dmtk
#endif
