// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT
#ifndef ALPS_DMRG_CHECKPOINT_CACHE_H
#define ALPS_DMRG_CHECKPOINT_CACHE_H

// Checkpoints of DMRG runs. A task writes its files to <task>.checkpoint/
// next to its task files from the start; a run that fails, by any means,
// leaves the directory behind and rerunning the job resumes from it. It is
// removed after a successful run.

#include "dmtk/filelist.h"
#include <alps/parameter.h>
#include <alps/params.hpp>
#include <alps/scheduler/signal.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <cerrno>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace dmrg_checkpoint {

namespace fs = boost::filesystem;

const char* const meta_magic = "ALPS DMRG checkpoint 1";
// reruns that may fail to commit anything after resuming before the
// checkpoint is taken to be what makes them fail
const int max_resume_attempts = 2;

// Last committed iteration of a regular sweep, and the state of the Lanczos
// start vector generator after it.
struct Position { int sweep; int dir; int iter; long random_state; };

typedef std::map<std::string, std::vector<double> > Series;

// Everything that makes a checkpoint invalid for another run. Parameters that
// only affect where or how results are written are left out, and so is SEED:
// parameter2xml and pyalps.writeInputFiles choose a new one on each call. A
// resumed run continues with the generator state of the run it resumes.
inline std::string fingerprint(const alps::Parameters& parms, const std::string& hamiltonian, const char* version)
{
  static const char* const ignored[] = {"CHECKPOINT", "TEMP_DIRECTORY", "VERBOSE", "SAVE_MPS", "SEED"};
  std::vector<std::pair<std::string, std::string> > values;
  for (alps::Parameters::const_iterator it = parms.begin(); it != parms.end(); ++it) {
    if (std::find(std::begin(ignored), std::end(ignored), it->key()) != std::end(ignored))
      continue;
    values.push_back(std::make_pair(it->key(), static_cast<std::string>(it->value())));
  }
  std::sort(values.begin(), values.end());
  std::ostringstream out;
  out << "dmrg " << version << "\n";
  for (const auto& v : values)
    out << v.first << "=" << v.second << "\n";
  out << "hamiltonian\n" << hamiltonian;
  return out.str();
}

inline std::string key(const std::string& fingerprint)
{
  std::uint64_t h = 14695981039346656037ull;
  for (unsigned char c : fingerprint) {
    h ^= c;
    h *= 1099511628211ull;
  }
  char buf[17];
  std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(h));
  return buf;
}

// <dir>/<base>.checkpoint for the task file <dir>/<base>.in.xml or
// <base>.out.xml; a task without a file uses the working directory.
inline fs::path location(const fs::path& task_file, const std::string& fingerprint)
{
  if (task_file.empty())
    return fs::current_path() / ("dmrg-" + key(fingerprint) + ".checkpoint");
  fs::path file = fs::absolute(task_file);
  std::string base = file.filename().string();
  for (const char* suffix : {".xml", ".h5", ".in", ".out"}) {
    std::size_t n = std::strlen(suffix);
    if (base.size() > n && base.compare(base.size() - n, n, suffix) == 0)
      base.erase(base.size() - n);
  }
  return file.parent_path() / (base + ".checkpoint");
}

inline std::string write_progress(const Position& at, const std::string& note, const Series& series)
{
  std::ostringstream out;
  out << "sweep " << at.sweep << "\ndir " << at.dir << "\niter " << at.iter
      << "\nrandom " << at.random_state << "\nnote " << note << "\nseries " << series.size() << "\n";
  char buf[32];
  for (const auto& s : series) {
    out << s.first << "\n" << s.second.size();
    for (double v : s.second) {
      std::snprintf(buf, sizeof(buf), "%.17g", v);
      out << " " << buf;
    }
    out << "\n";
  }
  return out.str();
}

inline bool read_progress(const std::string& record, Position& at, std::string& note, Series& series)
{
  std::istringstream in(record);
  std::string word;
  std::size_t count;
  if (!(in >> word >> at.sweep) || word != "sweep" || !(in >> word >> at.dir) || word != "dir"
      || !(in >> word >> at.iter) || word != "iter" || !(in >> word >> at.random_state) || word != "random"
      || !(in >> word) || word != "note")
    return false;
  in.ignore(1);
  std::getline(in, note);
  if (!(in >> word >> count) || word != "series")
    return false;
  Series result;
  for (std::size_t i = 0; i < count; ++i) {
    std::string name;
    std::size_t n;
    if (!std::getline(in >> std::ws, name) || !(in >> n))
      return false;
    std::vector<double>& values = result[name];
    for (std::size_t j = 0; j < n; ++j) {
      // strtod also reads the nan and inf that %g writes
      if (!(in >> word))
        return false;
      values.push_back(std::strtod(word.c_str(), 0));
    }
  }
  series.swap(result);
  return true;
}

// The fingerprint of the run that wrote the checkpoint; false if there is no
// readable meta file of this version.
inline bool read_meta(const fs::path& entry, std::string& fingerprint)
{
  std::ifstream in((entry / "meta").string().c_str(), std::ios::in | std::ios::binary);
  std::string line;
  if (!std::getline(in, line) || line != meta_magic || !std::getline(in, line) || line != "fingerprint")
    return false;
  std::ostringstream rest;
  rest << in.rdbuf();
  fingerprint = rest.str();
  return true;
}

class Cache {
public:
  Cache() {}
  Cache(const Cache&) = delete;
  Cache& operator=(const Cache&) = delete;
  ~Cache() { release(); }

  void configure(const alps::params& settings) { enabled_ = settings["CHECKPOINT"].as<bool>(); }

  bool enabled() const { return enabled_; }
  bool active() const { return active_; }
  const fs::path& directory() const { return entry_; }
  // Progress record of the last commit, empty if none.
  const std::string& record() const { return record_; }

  // Take the checkpoint of the task in task_file and point dmtk::tmp_files at
  // it. False if the run has to go on without checkpoints; a checkpoint
  // problem never fails a run.
  bool open(const std::string& fingerprint, const fs::path& task_file)
  {
    fingerprint_ = fingerprint;
    try {
      entry_ = location(task_file, fingerprint_);
      fs::create_directories(entry_);
      lock_fd_ = ::open((entry_ / "lock").c_str(), O_RDWR | O_CREAT, 0644);
      if (lock_fd_ < 0)
        throw std::runtime_error(std::string("cannot create the lock file: ") + std::strerror(errno));
      if (::flock(lock_fd_, LOCK_EX | LOCK_NB) != 0) {
        if (errno == EWOULDBLOCK) {
          std::cerr << "ALPS DMRG: another run of this task is using " << entry_.string()
                    << "; this run continues without checkpoints." << std::endl;
          release();
          return false;
        }
        // some network filesystems have no flock
        std::cerr << "ALPS DMRG: cannot lock " << (entry_ / "lock").string() << " (" << std::strerror(errno)
                  << "); do not run this task twice at the same time." << std::endl;
      }
      std::string written_by;
      bool readable = read_meta(entry_, written_by);
      if (!readable || written_by != fingerprint_) {
        if (fs::exists(entry_ / "meta"))
          std::cout << "ALPS DMRG: discarding checkpoint " << entry_.string() << ": "
                    << (readable ? "it belongs to other parameters" : "unreadable or from another version") << std::endl;
        clear_entry();
      }
      write_meta();
      try {
        record_ = dmtk::tmp_files.open_entry(entry_.c_str());
        if (!record_.empty() && !count_attempt()) {
          std::cout << "ALPS DMRG: discarding checkpoint " << entry_.string() << ": the last "
                    << max_resume_attempts << " runs that resumed from it failed before the next checkpoint" << std::endl;
          clear_entry();
          write_meta();
          record_ = dmtk::tmp_files.open_entry(entry_.c_str());
        }
      } catch (const std::exception& e) {
        std::cout << "ALPS DMRG: discarding checkpoint " << entry_.string() << ": " << e.what() << std::endl;
        clear_entry();
        write_meta();
        record_ = dmtk::tmp_files.open_entry(entry_.c_str());
      }
    } catch (const std::exception& e) {
      std::cerr << "ALPS DMRG: cannot use checkpoint directory " << entry_.string() << " (" << e.what()
                << "); this run continues without checkpoints." << std::endl;
      dmtk::tmp_files.cleanup();
      release();
      return false;
    }
    std::cout << "ALPS DMRG checkpoint files will be written to " << entry_.string() << std::endl;
    active_ = true;
    return true;
  }

  // Called after every regular sweep iteration.
  void commit(const std::string& record)
  {
    ++commits_;
    simulate_failure(true);
    dmtk::tmp_files.commit(record);
    record_ = record;
    simulate_failure(false);
    // The scheduler only looks at signals between tasks. SIGUSR1 and SIGUSR2
    // ask for a checkpoint, which the commit just made.
    typedef alps::scheduler::SignalHandler Signals;
    for (Signals::SignalInfo signal = Signals()(); signal != Signals::NOSIGNAL; signal = Signals()())
      if (signal == Signals::TERMINATE)
        throw std::runtime_error("DMRG run interrupted by a signal");
  }

  void succeeded()
  {
    if (!active_) return;
    active_ = false;
    boost::system::error_code error;
    fs::remove_all(entry_, error);
    if (error)
      std::cerr << "ALPS DMRG: could not remove checkpoint " << entry_.string() << ": " << error.message() << std::endl;
    release();
  }

  void failed()
  {
    if (!active_) return;
    active_ = false;
    release();
    Position at;
    std::string note;
    Series series;
    if (!read_progress(record_, at, note, series)) {
      boost::system::error_code error;
      fs::remove_all(entry_, error);
      std::cerr << "DMRG run failed before its first checkpoint; there is nothing to resume." << std::endl;
      return;
    }
    std::cerr << "DMRG run failed. Checkpoint kept in " << entry_.string() << " (" << note << ").\n"
              << "Rerunning this job with the same parameters resumes from it; "
                 "delete the directory to start over." << std::endl;
  }

private:
  // Count a run resuming from the last commit; false once max_resume_attempts
  // runs have done so without committing again.
  bool count_attempt()
  {
    const fs::path file = entry_ / "attempts";
    const unsigned long generation = dmtk::tmp_files.last_generation();
    unsigned long seen = 0;
    int count = 0;
    std::ifstream in(file.string().c_str());
    if (!(in >> seen >> count) || seen != generation)
      count = 0;
    in.close();
    if (count >= max_resume_attempts)
      return false;
    dmtk::detail::write_atomically(file, std::to_string(generation) + " " + std::to_string(count + 1) + "\n");
    return true;
  }

  // Keep the lock file: removing it would let another run lock a new one.
  void clear_entry()
  {
    std::vector<fs::path> contents;
    for (fs::directory_iterator it(entry_), end; it != end; ++it)
      if (it->path().filename() != "lock")
        contents.push_back(it->path());
    for (const fs::path& p : contents)
      fs::remove_all(p);
  }

  void write_meta()
  {
    dmtk::detail::write_atomically(entry_ / "meta", std::string(meta_magic) + "\nfingerprint\n" + fingerprint_);
  }

  void release()
  {
    if (lock_fd_ >= 0) {
      ::flock(lock_fd_, LOCK_UN);
      ::close(lock_fd_);
      lock_fd_ = -1;
    }
  }

  // ALPS_DMRG_CHECKPOINT_TEST=<action>:<n> fails the n-th commit of this
  // process for the tests: kill-staged (SIGKILL before the commit), or kill,
  // throw, term (SIGTERM) or usr1 (SIGUSR1) after it; sleep holds the
  // checkpoint for 3 seconds.
  void simulate_failure(bool before_commit) const
  {
    const char* test = std::getenv("ALPS_DMRG_CHECKPOINT_TEST");
    if (!test) return;
    std::string spec(test);
    std::string::size_type colon = spec.find(':');
    if (colon == std::string::npos || std::atoi(spec.c_str() + colon + 1) != commits_) return;
    std::string action = spec.substr(0, colon);
    if (before_commit) {
      if (action == "kill-staged") std::raise(SIGKILL);
      return;
    }
    if (action == "kill") std::raise(SIGKILL);
    if (action == "throw") throw std::runtime_error("simulated DMRG failure");
    if (action == "term") std::raise(SIGTERM);
    if (action == "usr1") std::raise(SIGUSR1);
    if (action == "sleep") ::sleep(3);
  }

  bool enabled_ = true;
  bool active_ = false;
  fs::path entry_;
  std::string fingerprint_;
  std::string record_;
  int lock_fd_ = -1;
  int commits_ = 0;
};

} // namespace dmrg_checkpoint
#endif
