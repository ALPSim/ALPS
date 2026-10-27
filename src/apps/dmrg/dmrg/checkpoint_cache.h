// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT
#ifndef ALPS_DMRG_CHECKPOINT_CACHE_H
#define ALPS_DMRG_CHECKPOINT_CACHE_H

// Checkpoint cache for DMRG runs. Each run writes its files to
// <CHECKPOINT_DIRECTORY>/<key>/ from the start; a run that fails, by any
// means, leaves the entry behind and a rerun with the same parameters resumes
// from it. Entries are deleted after a successful run, 72 hours after their
// last update, or oldest first when the cache exceeds CHECKPOINT_MAX_SIZE.

#include "dmtk/filelist.h"
#include <alps/parameter.h>
#include <alps/params.hpp>
#include <alps/scheduler/signal.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <condition_variable>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#include <fcntl.h>
#include <poll.h>
#include <sys/file.h>
#include <sys/statvfs.h>
#include <unistd.h>

namespace dmrg_checkpoint {

namespace fs = boost::filesystem;

const char* const meta_magic = "ALPS DMRG checkpoint 1";
const double ttl_hours = 72;
const double default_max_size_gb = 10;
const double min_free_gb = 1;
// flock may not see a lock taken on another host of a network filesystem, so
// a run touches its meta file this often and an entry of another host counts
// as live until it has not been touched for other_host_live_hours
const int heartbeat_seconds = 300;
const double other_host_live_hours = 1;
// reruns that may fail to commit anything after resuming before the
// checkpoint is taken to be what makes them fail
const int max_resume_attempts = 2;
// how long the prompt for a checkpoint location waits for an answer
const int prompt_timeout_seconds = 120;

// Last committed iteration of a regular sweep, and the state of the Lanczos
// start vector generator after it.
struct Position { int sweep; int dir; int iter; long random_state; };

typedef std::map<std::string, std::vector<double> > Series;

// Host name, with the machine id where there is one: containers on different
// nodes can share a host name.
inline std::string hostname()
{
  char buf[256] = {0};
  gethostname(buf, sizeof(buf) - 1);
  std::string host(buf);
  for (const char* file : {"/etc/machine-id", "/var/lib/dbus/machine-id"}) {
    std::ifstream in(file);
    std::string id;
    if (in >> id && !id.empty())
      return host + ":" + id;
  }
  return host;
}

// Expand a leading ~ and $NAME or ${NAME}, which a parameter file passes on
// literally. False if a variable is not set.
inline bool expand_path(const std::string& path, std::string& result, std::string& error)
{
  result.clear();
  std::string::size_type i = 0;
  if (!path.empty() && path[0] == '~' && (path.size() == 1 || path[1] == '/')) {
    const char* home = std::getenv("HOME");
    if (!home || !*home) {
      error = "HOME is not set";
      return false;
    }
    result = home;
    i = 1;
  }
  while (i < path.size()) {
    if (path[i] != '$') {
      result += path[i++];
      continue;
    }
    std::string name;
    std::string::size_type j = i + 1;
    if (j < path.size() && path[j] == '{') {
      std::string::size_type close = path.find('}', j);
      if (close == std::string::npos) {
        error = "unterminated ${ in " + path;
        return false;
      }
      name = path.substr(j + 1, close - j - 1);
      j = close + 1;
    } else {
      while (j < path.size() && (std::isalnum(static_cast<unsigned char>(path[j])) || path[j] == '_'))
        name += path[j++];
    }
    if (name.empty()) {
      result += path[i++];
      continue;
    }
    const char* value = std::getenv(name.c_str());
    if (!value) {
      error = "environment variable " + name + " is not set";
      return false;
    }
    result += value;
    i = j;
  }
  return true;
}

// Everything that makes a checkpoint invalid for another run. Parameters that
// only affect where or how results are written are left out, and so is SEED:
// parameter2xml and pyalps.writeInputFiles choose a new one on each call. A
// resumed run continues with the generator state of the run it resumes.
inline std::string fingerprint(const alps::Parameters& parms, const std::string& hamiltonian, const char* version)
{
  static const char* const ignored[] = {"CHECKPOINT", "CHECKPOINT_DIRECTORY", "CHECKPOINT_MAX_SIZE",
                                        "TEMP_DIRECTORY", "VERBOSE", "SAVE_MPS", "SEED"};
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

// FNV-1a; the full fingerprint stored in the entry is what is compared.
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

struct Meta {
  bool valid = false;
  std::string host;
  long pid = 0;
  std::string fingerprint;
};

// The first line marks a directory as an entry; nothing else is ever deleted.
inline Meta read_meta(const fs::path& entry)
{
  Meta meta;
  std::ifstream in((entry / "meta").string().c_str(), std::ios::in | std::ios::binary);
  std::string line, word;
  if (!std::getline(in, line) || line != meta_magic)
    return meta;
  if (!(in >> word >> meta.host) || word != "host" || !(in >> word >> meta.pid) || word != "pid")
    return meta;
  in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  if (!std::getline(in, line) || line != "fingerprint")
    return meta;
  std::ostringstream rest;
  rest << in.rdbuf();
  meta.fingerprint = rest.str();
  meta.valid = true;
  return meta;
}

inline std::time_t heartbeat(const fs::path& entry)
{
  boost::system::error_code error;
  std::time_t t = fs::last_write_time(entry / "meta", error);
  if (error) t = 0;
  std::time_t progress = fs::last_write_time(entry / "progress", error);
  return error ? t : std::max(t, progress);
}

inline std::uintmax_t entry_size(const fs::path& entry)
{
  std::uintmax_t total = 0;
  boost::system::error_code error;
  for (fs::recursive_directory_iterator it(entry, error), end; !error && it != end; it.increment(error)) {
    std::uintmax_t size = fs::file_size(it->path(), error);
    if (!error) total += size;
    error.clear();
  }
  return total;
}

// Held by a running process, or that process is seen alive on this host.
inline bool live(const fs::path& entry, const Meta& meta, std::time_t now)
{
  bool local = meta.host == hostname();
  int fd = ::open((entry / "lock").c_str(), O_RDWR);
  if (fd >= 0) {
    int locked = ::flock(fd, LOCK_EX | LOCK_NB);
    int error = errno;
    ::close(fd);
    if (locked != 0 && error == EWOULDBLOCK)
      return true;
    if (locked == 0 && local)
      return false;
  }
  if (local)
    return meta.pid > 0 && (::kill(static_cast<pid_t>(meta.pid), 0) == 0 || errno == EPERM);
  return std::difftime(now, heartbeat(entry)) < other_host_live_hours * 3600;
}

inline std::string format_time(std::time_t t)
{
  char buf[32];
  std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M", std::localtime(&t));
  return buf;
}

inline fs::path default_directory()
{
  const char* xdg = std::getenv("XDG_CACHE_HOME");
  if (xdg && *xdg && fs::path(xdg).is_absolute())
    return fs::path(xdg) / "alps" / "dmrg";
  const char* home = std::getenv("HOME");
  if (home && *home)
    return fs::path(home) / ".cache" / "alps" / "dmrg";
  // per user: another user's directory is not writable
  return fs::temp_directory_path() / ("alps-dmrg-" + std::to_string(::getuid()));
}

class Cache {
public:
  Cache() {}
  Cache(const Cache&) = delete;
  Cache& operator=(const Cache&) = delete;
  ~Cache() { release(); }

  void configure(const alps::params& settings)
  {
    enabled_ = settings["CHECKPOINT"].as<bool>();
    root_error_.clear();
    if (settings.exists("CHECKPOINT_DIRECTORY")) {
      std::string dir;
      if (expand_path(settings["CHECKPOINT_DIRECTORY"].as<std::string>(), dir, root_error_) && !dir.empty())
        root_ = fs::absolute(dir);
      else if (root_error_.empty())
        root_error_ = "CHECKPOINT_DIRECTORY is empty";
    } else {
      root_ = default_directory();
    }
    max_size_ = settings["CHECKPOINT_MAX_SIZE"].as<double>() * 1e9;
  }

  bool enabled() const { return enabled_; }
  bool active() const { return active_; }
  // Progress record of the last commit, empty if none.
  const std::string& record() const { return record_; }

  // Take the entry for this run and point dmtk::tmp_files at it. False if the
  // run has to go on without checkpoints; a checkpoint problem never fails a run.
  bool open(const std::string& fingerprint)
  {
    fingerprint_ = fingerprint;
    if (!root_error_.empty()) {
      std::cerr << "ALPS DMRG: cannot use CHECKPOINT_DIRECTORY: " << root_error_ << std::endl;
      if (!choose_location())
        return false;
    } else if (!usable(root_) && !choose_location()) {
      return false;
    }
    try {
      collect_garbage();
    } catch (const std::exception& e) {
      std::cerr << "ALPS DMRG: checkpoint cleanup failed: " << e.what() << std::endl;
    }
    entry_ = root_ / key(fingerprint_);
    try {
      fs::create_directories(entry_);
      lock_fd_ = ::open((entry_ / "lock").c_str(), O_RDWR | O_CREAT, 0644);
      if (lock_fd_ < 0)
        throw std::runtime_error("cannot create " + (entry_ / "lock").string());
      if (::flock(lock_fd_, LOCK_EX | LOCK_NB) != 0 && errno == EWOULDBLOCK) {
        std::cerr << "ALPS DMRG: another run with the same parameters is using " << entry_.string()
                  << "; this run continues without checkpoints." << std::endl;
        release();
        return false;
      }
      Meta meta = read_meta(entry_);
      if (!meta.valid || meta.fingerprint != fingerprint_) {
        if (fs::exists(entry_ / "meta"))
          std::cout << "ALPS DMRG: discarding checkpoint " << entry_.string() << ": "
                    << (meta.valid ? "it belongs to other parameters" : "unreadable or from another version") << std::endl;
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
    start_heartbeat();
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
    stop_heartbeat();
    boost::system::error_code error;
    fs::remove_all(entry_, error);
    release();
  }

  void failed()
  {
    if (!active_) return;
    active_ = false;
    stop_heartbeat();
    Position at;
    std::string note;
    Series series;
    if (!read_progress(record_, at, note, series)) {
      boost::system::error_code error;
      fs::remove_all(entry_, error);
      release();
      std::cerr << "DMRG run failed before its first checkpoint; there is nothing to resume." << std::endl;
      return;
    }
    release();
    std::time_t expiry = heartbeat(entry_) + static_cast<std::time_t>(ttl_hours * 3600);
    std::cerr << "DMRG run failed. Checkpoint kept at " << entry_.string() << " (" << note << ").\n"
              << "Rerunning with the same parameters will resume automatically.\n"
              << "This checkpoint will be deleted after " << format_time(expiry)
              << " (72 hours after its last update)." << std::endl;
  }

private:
  bool usable(const fs::path& dir) const
  {
    boost::system::error_code error;
    fs::create_directories(dir, error);
    if (error) {
      std::cerr << "ALPS DMRG: cannot create checkpoint directory " << dir.string() << ": " << error.message() << std::endl;
      return false;
    }
    struct statvfs info;
    if (::statvfs(dir.c_str(), &info) == 0) {
      double free_gb = static_cast<double>(info.f_bavail) * info.f_frsize / 1e9;
      if (free_gb < min_free_gb) {
        std::cerr << "ALPS DMRG: only " << free_gb << " GB free in " << dir.string() << " for checkpoints." << std::endl;
        return false;
      }
    }
    return true;
  }

  static void continue_without()
  {
    std::cerr << "ALPS DMRG: continuing without checkpoints, so a failed run cannot be resumed. "
                 "Set CHECKPOINT_DIRECTORY to a location with more space (for example scratch), "
                 "or CHECKPOINT=0 to turn checkpoints off." << std::endl;
  }

  // A line from the terminal, or false after prompt_timeout_seconds or at its end.
  static bool read_answer(std::string& answer)
  {
    struct pollfd in = {STDIN_FILENO, POLLIN, 0};
    int ready;
    do {
      ready = ::poll(&in, 1, prompt_timeout_seconds * 1000);
    } while (ready < 0 && errno == EINTR);
    if (ready <= 0) {
      std::cerr << "\nALPS DMRG: no answer in " << prompt_timeout_seconds << " seconds." << std::endl;
      return false;
    }
    return static_cast<bool>(std::getline(std::cin, answer));
  }

  // Ask where checkpoints should go when someone is at the terminal: not in a
  // batch job, where nobody answers, nor in the background, where reading the
  // terminal stops the process.
  bool choose_location()
  {
    if (!::isatty(STDIN_FILENO) || !::isatty(STDERR_FILENO) || ::tcgetpgrp(STDIN_FILENO) != ::getpgrp()) {
      continue_without();
      return false;
    }
    for (;;) {
      std::cerr << "Where should DMRG checkpoints go?\n"
                << "  1) the current directory (" << fs::current_path().string() << ")\n"
                << "  2) another directory\n"
                << "  3) nowhere: continue without checkpoints\n> " << std::flush;
      std::string answer;
      if (!read_answer(answer) || answer == "3") {
        continue_without();
        return false;
      }
      fs::path dir;
      if (answer == "1") {
        dir = fs::current_path();
      } else if (answer == "2") {
        std::cerr << "Directory: " << std::flush;
        if (!read_answer(answer)) {
          continue_without();
          return false;
        }
        std::string expanded, error;
        if (answer.empty())
          continue;
        if (!expand_path(answer, expanded, error)) {
          std::cerr << "ALPS DMRG: " << error << std::endl;
          continue;
        }
        dir = fs::absolute(expanded);
      } else {
        continue;
      }
      if (usable(dir)) {
        root_ = dir;
        root_error_.clear();
        std::cerr << "To skip this question, set CHECKPOINT_DIRECTORY=\"" << dir.string() << "\"" << std::endl;
        return true;
      }
    }
  }

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

  void start_heartbeat()
  {
    stop_heartbeat();
    stop_ = false;
    const fs::path meta = entry_ / "meta";
    heartbeat_ = std::thread([this, meta]() {
      std::unique_lock<std::mutex> lock(mutex_);
      while (!wake_.wait_for(lock, std::chrono::seconds(heartbeat_seconds), [this]() { return stop_; })) {
        boost::system::error_code error;
        fs::last_write_time(meta, std::time(0), error);
      }
    });
  }

  void stop_heartbeat()
  {
    if (!heartbeat_.joinable()) return;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_ = true;
    }
    wake_.notify_all();
    heartbeat_.join();
  }

  void collect_garbage()
  {
    struct Candidate { std::time_t beat; std::uintmax_t size; fs::path path; };
    std::vector<fs::path> entries;
    for (fs::directory_iterator it(root_), end; it != end; ++it)
      if (fs::is_directory(it->path()) && read_meta(it->path()).valid)
        entries.push_back(it->path());

    const std::time_t now = std::time(0);
    const std::string current = key(fingerprint_);
    std::vector<Candidate> evictable;
    std::uintmax_t total = 0;
    for (const fs::path& entry : entries) {
      std::uintmax_t size = entry_size(entry);
      total += size;
      if (entry.filename().string() == current || live(entry, read_meta(entry), now))
        continue;
      std::time_t beat = heartbeat(entry);
      if (std::difftime(now, beat) > ttl_hours * 3600) {
        remove_entry(entry, "not updated for 72 hours");
        total -= size;
      } else {
        evictable.push_back(Candidate{beat, size, entry});
      }
    }
    std::sort(evictable.begin(), evictable.end(),
              [](const Candidate& a, const Candidate& b) { return a.beat < b.beat; });
    for (const Candidate& c : evictable) {
      if (total <= max_size_) break;
      remove_entry(c.path, "checkpoint cache larger than CHECKPOINT_MAX_SIZE");
      total -= c.size;
    }
  }

  static void remove_entry(const fs::path& entry, const char* reason)
  {
    boost::system::error_code error;
    fs::remove_all(entry, error);
    if (error)
      std::cerr << "ALPS DMRG: could not remove checkpoint " << entry.string() << ": " << error.message() << std::endl;
    else
      std::cout << "ALPS DMRG: removed checkpoint " << entry.string() << " (" << reason << ")" << std::endl;
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
    std::ostringstream out;
    out << meta_magic << "\nhost " << hostname() << "\npid " << ::getpid() << "\nfingerprint\n" << fingerprint_;
    dmtk::detail::write_atomically(entry_ / "meta", out.str());
  }

  void release()
  {
    stop_heartbeat();
    if (lock_fd_ >= 0) {
      ::flock(lock_fd_, LOCK_UN);
      ::close(lock_fd_);
      lock_fd_ = -1;
    }
  }

  // ALPS_DMRG_CHECKPOINT_TEST=<action>:<n> fails the n-th commit of this
  // process for the tests: kill-staged (SIGKILL before the commit), or kill,
  // throw, term (SIGTERM) or usr1 (SIGUSR1) after it.
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
  }

  bool enabled_ = true;
  bool active_ = false;
  fs::path root_;
  std::string root_error_;
  fs::path entry_;
  double max_size_ = default_max_size_gb * 1e9;
  std::string fingerprint_;
  std::string record_;
  int lock_fd_ = -1;
  int commits_ = 0;
  std::thread heartbeat_;
  std::mutex mutex_;
  std::condition_variable wake_;
  bool stop_ = false;
};

} // namespace dmrg_checkpoint
#endif
