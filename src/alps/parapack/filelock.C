#include <chrono>
#include <thread>
/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "filelock.h"
#include <boost/filesystem/operations.hpp>
#include <boost/throw_exception.hpp>
#include <iostream>
#include <cstdio>
#include <cerrno>
#include <system_error>
#include <stdexcept>

namespace alps {

filelock::filelock() : auto_release_(true), is_locking_(false) {}

filelock::filelock(boost::filesystem::path const& file, bool lock_now, int wait,
  bool auto_release) : auto_release_(auto_release), is_locking_(false) {
  set_file(file);
  if (lock_now) lock(wait);
}

filelock::~filelock() {
  if (is_locking_) {
    if (!auto_release_)
      std::cerr << "Warning: lock for \"" << file_ << "\" is being removed\n";
    release();
  }
}

void filelock::set_file(boost::filesystem::path const& file) {
  if (is_locking_) {
    std::cerr << "Warning: lock for \"" << file_ << "\" is being removed\n";
    release();
  }
  file_ = file.string();
  lock_ = file.parent_path() / (file.filename().string() + ".lck");
}

void filelock::lock(int wait) {
  if (is_locking_) {
    std::cerr << "Error: file \"" << file_ << "\" is already locked.\n";
    boost::throw_exception(std::logic_error("filelock"));
  }
  for (int i = 0; wait < 0 || i < wait+1; ++i) {
    if (i != 0) {
      std::cerr << "Warning: file \"" << file_ << "\" is locked.  Still trying.\n";
      std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    if (auto* stream = std::fopen(lock_.string().c_str(), "wx")) {
      std::fclose(stream);
      is_locking_ = true;
      break;
    }
    if (errno != EEXIST)
      throw std::system_error(errno, std::generic_category(), "create lock file");
  }
  if (!is_locking_) {
    std::cerr << "Error: lock for file \"" << file_ << "\" failed.\n";
    boost::throw_exception(std::logic_error("filelock"));
  }
}

void filelock::release() {
  if (!is_locking_) {
    std::cerr << "Error: file \"" << file_ << "\" is not locked\n";
    boost::throw_exception(std::logic_error("filelock"));
  }
  remove(lock_);
  is_locking_ = false;
}

bool filelock::locking() const { return is_locking_; }

bool filelock::locked() const { return is_locking_ || exists(lock_); }

} // end namespace alps
