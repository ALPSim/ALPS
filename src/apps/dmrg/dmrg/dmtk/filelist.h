// Copyright (C) 2006-2010 Adrian Feiguin
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT
#ifndef ALPS_DMTK_FILELIST_H
#define ALPS_DMTK_FILELIST_H

#include <alps/utilities/temporary_filename.hpp>
#include <boost/filesystem.hpp>
#include <boost/throw_exception.hpp>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <map>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

namespace dmtk {

namespace detail {

// Flush a file (or directory) to disk; ofstream gives no access to its descriptor.
inline bool sync_path(const boost::filesystem::path& path)
{
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0) return false;
  bool ok = ::fsync(fd) == 0;
  ::close(fd);
  return ok;
}

// Replace path by contents so that a crash leaves either the old or the new file.
inline void write_atomically(const boost::filesystem::path& path, const std::string& contents)
{
  boost::filesystem::path tmp = path;
  tmp += ".tmp";
  {
    std::ofstream out(tmp.string().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    out << contents;
    out.close();
    if (!out || !sync_path(tmp))
      boost::throw_exception(std::runtime_error("could not write " + tmp.string()));
  }
  boost::filesystem::rename(tmp, path);
  sync_path(path.parent_path()); // not supported on every filesystem
}

} // namespace detail

// Own only the files allocated by this registry, never the user's scratch
// directory or files belonging to another calculation.
//
// In checkpoint mode (open_entry) files keep their own names in the entry
// directory so that a rerun can find them. Writes go to
// <name>.new.<generation> and replace the committed files only in commit(),
// so a run killed at any point leaves the last committed state intact.
class FileList {
  public:
    FileList() : temp_dir(boost::filesystem::current_path()) {}
    explicit FileList(const char* dir) { set_temp_dir(dir); }
    FileList(const FileList&) = delete;
    FileList& operator=(const FileList&) = delete;
    ~FileList() { cleanup(); }

    void set_temp_dir(const char* dir)
    {
      if (in_entry) cleanup();
      // Store absolute paths so a later working-directory change cannot make
      // cleanup remove a different file or lose track of the original one.
      boost::filesystem::path path = boost::filesystem::absolute(dir);
      if (!boost::filesystem::is_directory(path))
        boost::throw_exception(std::runtime_error("ALPS DMRG temporary directory does not exist: " + path.string()));
      temp_dir = path;
      std::cout << "ALPS DMRG temporary files will be written to " << temp_dir.string() << std::endl;
    }

    // Use the checkpoint entry dir and return the caller's part of its last
    // progress record, empty if nothing was committed. Staged files of a
    // commit whose record was written are promoted, all others discarded.
    // Throws if a committed file is missing or has the wrong size: DMTK only
    // prints an error for a file it cannot read and goes on with an empty block.
    std::string open_entry(const char* dir)
    {
      cleanup();
      temp_dir = boost::filesystem::absolute(dir);
      in_entry = true;
      generation = 0;
      _committed.clear();
      std::string record;
      const boost::filesystem::path progress = temp_dir / "progress";
      if (boost::filesystem::exists(progress)) {
        std::ifstream in(progress.string().c_str(), std::ios::in | std::ios::binary);
        std::string magic, word;
        std::size_t count = 0;
        std::getline(in, magic);
        if (magic != progress_magic() ||!(in >> word >> generation) || word != "generation"
            || !(in >> word >> count) || word != "staged")
          boost::throw_exception(std::runtime_error("unreadable progress record " + progress.string()));
        for (std::size_t i = 0; i < count; ++i) {
          std::string name;
          in >> name;
          boost::filesystem::path staged = temp_dir / (name + suffix(generation));
          if (boost::filesystem::exists(staged))
            boost::filesystem::rename(staged, temp_dir / name);
        }
        if (!in || !(in >> word >> count) || word != "files")
          boost::throw_exception(std::runtime_error("unreadable progress record " + progress.string()));
        for (std::size_t i = 0; i < count; ++i) {
          std::string name;
          std::uintmax_t size;
          if (!(in >> name >> size))
            boost::throw_exception(std::runtime_error("unreadable progress record " + progress.string()));
          boost::system::error_code error;
          std::uintmax_t actual = boost::filesystem::file_size(temp_dir / name, error);
          if (error)
            boost::throw_exception(std::runtime_error("checkpoint file " + name + " is missing"));
          if (actual != size)
            boost::throw_exception(std::runtime_error("checkpoint file " + name + " has "
              + std::to_string(actual) + " bytes instead of " + std::to_string(size)));
          _committed[name] = size;
        }
        in.ignore(1);
        std::ostringstream rest;
        rest << in.rdbuf();
        record = rest.str();
        ++generation;
      }
      for (boost::filesystem::directory_iterator it(temp_dir), end; it != end; ++it)
        if (it->path().filename().string().find(".new.") != std::string::npos)
          boost::filesystem::remove(it->path());
      return record;
    }

    bool checkpointing() const { return in_entry; }

    const char* get_filename(const char* input)
    {
      auto old = _tmp_filenames.find(input);
      if (old != _tmp_filenames.end()) return old->second.c_str();
      std::string filename;
      if (in_entry) {
        filename = (temp_dir / input).string();
      } else {
        std::string name(input);
        name = name.substr(0, name.find_first_of('.'));
        filename = alps::temporary_filename((temp_dir / name).string());
        std::cout << "Creating temp file " << filename << std::endl;
      }
      auto entry = _tmp_filenames.emplace(input, filename);
      return entry.first->second.c_str();
    }

    // Checkpointed state: reads see this run's uncommitted writes.
    const char* get_read_filename(const char* input)
    {
      auto staged = _staged.find(input);
      if (staged != _staged.end()) return staged->second.c_str();
      return get_filename(input);
    }

    const char* get_write_filename(const char* input)
    {
      if (!in_entry) return get_filename(input);
      auto staged = _staged.find(input);
      if (staged != _staged.end()) return staged->second.c_str();
      auto entry = _staged.emplace(input, (temp_dir / (input + suffix(generation))).string());
      return entry.first->second.c_str();
    }

    void write_failed(const char* input) { _failed_writes.push_back(input); }

    // Make all staged files and the caller's record the state a rerun resumes from.
    void commit(const std::string& record)
    {
      if (!_failed_writes.empty())
        boost::throw_exception(std::runtime_error("could not write checkpoint file " + _failed_writes.front()
          + " in " + temp_dir.string() + " (disk full or over quota?)"));
      std::ostringstream out;
      out << progress_magic() <<"\ngeneration " << generation << "\nstaged " << _staged.size() << "\n";
      std::map<std::string, std::uintmax_t> committed = _committed;
      for (const auto& entry : _staged) {
        boost::system::error_code error;
        committed[entry.first] = boost::filesystem::file_size(entry.second, error);
        if (error || !detail::sync_path(entry.second))
          boost::throw_exception(std::runtime_error("could not write checkpoint file " + entry.second));
        out << entry.first << "\n";
      }
      out << "files " << committed.size() << "\n";
      for (const auto& entry : committed)
        out << entry.first << " " << entry.second << "\n";
      out << record;
      detail::write_atomically(temp_dir / "progress", out.str());
      for (const auto& entry : _staged)
        boost::filesystem::rename(entry.second, temp_dir / entry.first);
      _staged.clear();
      _committed.swap(committed);
      ++generation;
    }

    // Generation of the last commit, or of the first one if there was none.
    unsigned long last_generation() const { return generation == 0 ? 0 : generation - 1; }

    void cleanup() noexcept
    {
      // A checkpoint entry outlives the run; only uncommitted writes go.
      for (const auto& entry : in_entry ? _staged : _tmp_filenames) {
        boost::system::error_code error;
        boost::filesystem::remove(entry.second, error);
        if (error)
          std::cerr << "ALPS DMRG could not remove temporary file " << entry.second
                    << ": " << error.message() << std::endl;
      }
      _tmp_filenames.clear();
      _staged.clear();
      _failed_writes.clear();
      in_entry = false;
    }

  private:
    static const char* progress_magic() { return "ALPS DMRG progress 1"; }
    static std::string suffix(unsigned long generation) { return ".new." + std::to_string(generation); }

    std::map<std::string, std::string> _tmp_filenames;
    std::map<std::string, std::string> _staged;
    // name and size of every file the last commit left in the entry
    std::map<std::string, std::uintmax_t> _committed;
    std::vector<std::string> _failed_writes;
    boost::filesystem::path temp_dir;
    bool in_entry = false;
    unsigned long generation = 0;
};

} // namespace dmtk
#endif
