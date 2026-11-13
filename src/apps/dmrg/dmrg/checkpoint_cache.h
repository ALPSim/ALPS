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
#include <boost/throw_exception.hpp>
#include <filesystem>
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
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace dmtk { extern FileList tmp_files; } // defined in dmtk/system.h

namespace dmrg_checkpoint {

namespace fs = std::filesystem;

const char* const meta_magic = "ALPS DMRG checkpoint 1";
// reruns that may fail to commit anything after resuming before the
// checkpoint is taken to be what makes them fail
const int max_resume_attempts = 2;

// Last committed iteration of a regular sweep, the state of the Lanczos
// start vector generator after it, and the SEED the run started with.
struct Position { int sweep; int dir; int iter; long random_state; long seed; };

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
      << "\nrandom " << at.random_state << "\nseed " << at.seed << "\nnote " << note << "\nseries " << series.size() << "\n";
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
      || !(in >> word >> at.seed) || word != "seed" || !(in >> word) || word != "note")
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

// Flush a file (or directory) to disk; ofstream gives no access to its descriptor.
inline bool sync_path(const fs::path& path)
{
  int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0) return false;
  bool ok = ::fsync(fd) == 0;
  ::close(fd);
  return ok;
}

// Replace path by contents so that a crash leaves either the old or the new file.
inline void write_atomically(const fs::path& path, const std::string& contents)
{
  fs::path tmp = path;
  tmp += ".tmp";
  {
    std::ofstream out(tmp.string().c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    out << contents;
    out.close();
    if (!out || !sync_path(tmp))
      boost::throw_exception(std::runtime_error("could not write " + tmp.string()));
  }
  fs::rename(tmp, path);
  sync_path(path.parent_path()); // not supported on every filesystem
}

// The files of a checkpoint entry, attached to a dmtk::FileList. Files keep
// their own names in the entry directory so that a rerun can find them.
// Writes go to <name>.new.<generation> and replace the committed files only
// in commit(), so a run killed at any point leaves the last committed state
// intact.
class Files : public dmtk::FileStore {
public:
  Files() {}
  Files(const Files&) = delete;
  Files& operator=(const Files&) = delete;
  ~Files() { close(); }

  // Use the entry dir and return the caller's part of its last progress
  // record, empty if nothing was committed. Staged files of a commit whose
  // record was written are promoted, all others discarded. Throws if a
  // committed file is missing or has the wrong size: DMTK only prints an
  // error for a file it cannot read and goes on with an empty block.
  std::string open(const fs::path& dir)
  {
    close();
    dir_ = fs::absolute(dir);
    generation_ = 0;
    committed_.clear();
    std::string record;
    const fs::path progress = dir_ / "progress";
    if (fs::exists(progress)) {
      std::ifstream in(progress.string().c_str(), std::ios::in | std::ios::binary);
      std::string magic, word;
      std::size_t count = 0;
      std::getline(in, magic);
      if (magic != progress_magic() || !(in >> word >> generation_) || word != "generation"
          || !(in >> word >> count) || word != "staged")
        boost::throw_exception(std::runtime_error("unreadable progress record " + progress.string()));
      for (std::size_t i = 0; i < count; ++i) {
        std::string name;
        in >> name;
        fs::path staged = dir_ / (name + suffix(generation_));
        if (fs::exists(staged))
          fs::rename(staged, dir_ / name);
      }
      if (!in || !(in >> word >> count) || word != "files")
        boost::throw_exception(std::runtime_error("unreadable progress record " + progress.string()));
      for (std::size_t i = 0; i < count; ++i) {
        std::string name;
        std::uintmax_t size;
        if (!(in >> name >> size))
          boost::throw_exception(std::runtime_error("unreadable progress record " + progress.string()));
        std::error_code error;
        std::uintmax_t actual = fs::file_size(dir_ / name, error);
        if (error)
          boost::throw_exception(std::runtime_error("checkpoint file " + name + " is missing"));
        if (actual != size)
          boost::throw_exception(std::runtime_error("checkpoint file " + name + " has "
            + std::to_string(actual) + " bytes instead of " + std::to_string(size)));
        committed_[name] = size;
      }
      in.ignore(1);
      std::ostringstream rest;
      rest << in.rdbuf();
      record = rest.str();
      ++generation_;
    }
    for (fs::directory_iterator it(dir_), end; it != end; ++it)
      if (it->path().filename().string().find(".new.") != std::string::npos)
        fs::remove(it->path());
    return record;
  }

  // Until list.cleanup() or close().
  void attach(dmtk::FileList& list)
  {
    list.attach(this);
    list_ = &list;
  }

  // Reads see this run's uncommitted writes.
  const char* get_read_filename(const char* input) override
  {
    auto staged = staged_.find(input);
    if (staged != staged_.end()) return staged->second.c_str();
    auto entry = paths_.emplace(input, (dir_ / input).string());
    return entry.first->second.c_str();
  }

  const char* get_write_filename(const char* input) override
  {
    auto staged = staged_.find(input);
    if (staged != staged_.end()) return staged->second.c_str();
    auto entry = staged_.emplace(input, (dir_ / (input + suffix(generation_))).string());
    return entry.first->second.c_str();
  }

  void write_failed(const char* input) override { failed_writes_.push_back(input); }

  // Make all staged files and the caller's record the state a rerun resumes from.
  void commit(const std::string& record)
  {
    if (!failed_writes_.empty())
      boost::throw_exception(std::runtime_error("could not write checkpoint file " + failed_writes_.front()
        + " in " + dir_.string() + " (disk full or over quota?)"));
    std::ostringstream out;
    out << progress_magic() << "\ngeneration " << generation_ << "\nstaged " << staged_.size() << "\n";
    std::map<std::string, std::uintmax_t> committed = committed_;
    for (const auto& entry : staged_) {
      std::error_code error;
      committed[entry.first] = fs::file_size(entry.second, error);
      if (error || !sync_path(entry.second))
        boost::throw_exception(std::runtime_error("could not write checkpoint file " + entry.second));
      out << entry.first << "\n";
    }
    out << "files " << committed.size() << "\n";
    for (const auto& entry : committed)
      out << entry.first << " " << entry.second << "\n";
    out << record;
    write_atomically(dir_ / "progress", out.str());
    for (const auto& entry : staged_)
      fs::rename(entry.second, dir_ / entry.first);
    staged_.clear();
    committed_.swap(committed);
    ++generation_;
  }

  // Generation of the last commit, or of the first one if there was none.
  unsigned long last_generation() const { return generation_ == 0 ? 0 : generation_ - 1; }

  // The entry outlives the run; only uncommitted writes go.
  void cleanup() noexcept override
  {
    for (const auto& entry : staged_) {
      std::error_code error;
      fs::remove(entry.second, error);
      if (error)
        std::cerr << "ALPS DMRG could not remove temporary file " << entry.second
                  << ": " << error.message() << std::endl;
    }
    paths_.clear();
    staged_.clear();
    failed_writes_.clear();
    list_ = nullptr;
  }

  // Detach from the FileList, which then calls cleanup().
  void close() noexcept
  {
    if (list_) list_->cleanup();
    else cleanup();
  }

private:
  static const char* progress_magic() { return "ALPS DMRG progress 1"; }
  static std::string suffix(unsigned long generation) { return ".new." + std::to_string(generation); }

  std::map<std::string, std::string> paths_;
  std::map<std::string, std::string> staged_;
  // name and size of every file the last commit left in the entry
  std::map<std::string, std::uintmax_t> committed_;
  std::vector<std::string> failed_writes_;
  fs::path dir_;
  unsigned long generation_ = 0;
  dmtk::FileList* list_ = nullptr;
};

// The ALPS handlers of SIGUSR1 and SIGUSR2 are one-shot: a second signal
// would kill the run. While checkpoints are kept, a handler that is not
// one-shot passes each of them on to the ALPS one.
typedef void (*signal_handler)(int);

inline signal_handler& alps_user_handler(int signal)
{
  static signal_handler handlers[2] = {0, 0};
  return handlers[signal == SIGUSR1 ? 0 : 1];
}

extern "C" inline void forward_user_signal(int signal)
{
  signal_handler handler = alps_user_handler(signal);
  if (handler) handler(signal);
}

class Cache {
public:
  Cache() {}
  Cache(const Cache&) = delete;
  Cache& operator=(const Cache&) = delete;
  ~Cache() { files_.close(); release(); }

  void configure(const alps::params& settings) { enabled_ = settings["CHECKPOINT"].as<bool>(); }

  bool enabled() const { return enabled_; }
  bool active() const { return active_; }
  const fs::path& directory() const { return entry_; }
  // Progress record of the last commit, empty if none.
  const std::string& record() const { return record_; }

  // Take the checkpoint of the task in task_file and attach its files to
  // dmtk::tmp_files. False if the run has to go on without checkpoints; a
  // checkpoint problem never fails a run.
  bool open(const std::string& fingerprint, const fs::path& task_file)
  {
    fingerprint_ = fingerprint;
    try {
      entry_ = location(task_file, fingerprint_);
      if (!ours(entry_)) {
        std::cerr << "ALPS DMRG: " << entry_.string() << " exists and is not a DMRG checkpoint; "
                     "this run continues without checkpoints." << std::endl;
        return false;
      }
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
        record_ = files_.open(entry_);
        if (!record_.empty() && !count_attempt()) {
          std::cout << "ALPS DMRG: discarding checkpoint " << entry_.string() << ": the last "
                    << max_resume_attempts << " runs that resumed from it failed before the next checkpoint" << std::endl;
          clear_entry();
          write_meta();
          record_ = files_.open(entry_);
        }
      } catch (const std::exception& e) {
        std::cout << "ALPS DMRG: discarding checkpoint " << entry_.string() << ": " << e.what() << std::endl;
        clear_entry();
        write_meta();
        record_ = files_.open(entry_);
      }
      files_.attach(dmtk::tmp_files);
    } catch (const std::exception& e) {
      std::cerr << "ALPS DMRG: cannot use checkpoint directory " << entry_.string() << " (" << e.what()
                << "); this run continues without checkpoints." << std::endl;
      dmtk::tmp_files.cleanup();
      release();
      return false;
    }
    std::cout << "ALPS DMRG checkpoint files will be written to " << entry_.string() << std::endl;
    active_ = true;
    forward_signals();
    return true;
  }

  // Called after every regular sweep iteration.
  void commit(const std::string& record)
  {
    ++commits_;
    simulate_failure(true);
    files_.commit(record);
    record_ = record;
    simulate_failure(false);
    simulate_signals();
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
    std::error_code error;
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
      std::error_code error;
      fs::remove_all(entry_, error);
      std::cerr << "DMRG run failed before its first checkpoint; there is nothing to resume." << std::endl;
      return;
    }
    std::cerr << "DMRG run failed. Checkpoint kept in " << entry_.string() << " (" << note << ").\n"
              << "Rerunning this job with the same parameters resumes from it; "
                 "delete the directory to start over." << std::endl;
  }

private:
  // A directory this code may clear: absent, or holding a meta file, or only
  // what open() writes before the meta file.
  static bool ours(const fs::path& entry)
  {
    if (!fs::exists(entry) || fs::exists(entry / "meta"))
      return true;
    if (!fs::is_directory(entry))
      return false;
    for (fs::directory_iterator it(entry), end; it != end; ++it) {
      std::string name = it->path().filename().string();
      if (name != "lock" && name != "meta.tmp")
        return false;
    }
    return true;
  }

  static const int* user_signals() { static const int signals[2] = {SIGUSR1, SIGUSR2}; return signals; }

  void forward_signals()
  {
    for (int i = 0; i < 2; ++i) {
      int signal = user_signals()[i];
      saved_[i] = ::sigaction(signal, 0, &actions_[i]) == 0 && !(actions_[i].sa_flags & SA_SIGINFO)
        && actions_[i].sa_handler != SIG_DFL && actions_[i].sa_handler != SIG_IGN;
      if (!saved_[i]) continue;
      alps_user_handler(signal) = actions_[i].sa_handler;
      struct sigaction forward;
      std::memset(&forward, 0, sizeof(forward));
      sigemptyset(&forward.sa_mask);
      forward.sa_flags = SA_RESTART;
      forward.sa_handler = &forward_user_signal;
      ::sigaction(signal, &forward, 0);
    }
  }

  // Back to the ALPS handlers, armed again.
  void restore_signals()
  {
    for (int i = 0; i < 2; ++i)
      if (saved_[i]) {
        ::sigaction(user_signals()[i], &actions_[i], 0);
        saved_[i] = false;
      }
  }

  // Count a run resuming from the last commit; false once max_resume_attempts
  // runs have done so without committing again.
  bool count_attempt()
  {
    const fs::path file = entry_ / "attempts";
    const unsigned long generation = files_.last_generation();
    unsigned long seen = 0;
    int count = 0;
    std::ifstream in(file.string().c_str());
    if (!(in >> seen >> count) || seen != generation)
      count = 0;
    in.close();
    if (count >= max_resume_attempts)
      return false;
    write_atomically(file, std::to_string(generation) + " " + std::to_string(count + 1) + "\n");
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
    write_atomically(entry_ / "meta", std::string(meta_magic) + "\nfingerprint\n" + fingerprint_);
  }

  void release()
  {
    restore_signals();
    if (lock_fd_ >= 0) {
      ::flock(lock_fd_, LOCK_UN);
      ::close(lock_fd_);
      lock_fd_ = -1;
    }
  }

  // ALPS_DMRG_CHECKPOINT_TEST=<action>:<n> fails the n-th commit of this
  // process for the tests: kill-staged (SIGKILL before the commit), or kill,
  // throw, term (SIGTERM) or usr1 (SIGUSR1) after it; sleep holds the
  // checkpoint for 3 seconds. See also simulate_signals.
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

  // usr1-twice:<n> raises SIGUSR1 twice after commit n, before either is seen.
  void simulate_signals() const
  {
    const char* test = std::getenv("ALPS_DMRG_CHECKPOINT_TEST");
    if (!test || std::strncmp(test, "usr1-twice:", 11) != 0 || std::atoi(test + 11) != commits_) return;
    std::raise(SIGUSR1);
    std::raise(SIGUSR1);
  }

  Files files_;
  bool enabled_ = true;
  bool active_ = false;
  fs::path entry_;
  std::string fingerprint_;
  std::string record_;
  int lock_fd_ = -1;
  int commits_ = 0;
  struct sigaction actions_[2];
  bool saved_[2] = {false, false};
};

} // namespace dmrg_checkpoint
#endif
