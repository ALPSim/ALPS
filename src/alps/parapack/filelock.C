/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#include "filelock.h"
#include <boost/filesystem/operations.hpp>
#include <boost/throw_exception.hpp>
#include <iostream>
#include <fcntl.h> // for open()
#include <sys/stat.h>
#include <stdexcept>

#include <alps/config.h>
#if defined(ALPS_HAVE_UNISTD_H)
# include <unistd.h>
#elif defined(ALPS_HAVE_WINDOWS_H)
# include <windows.h>
#endif
#if defined(ALPS_HAVE_SYS_STAT_H)
# include <sys/stat.h>
#endif
#if defined(ALPS_HAVE_WINDOWS_H)
# include <io.h>
#endif

namespace alps {

filelock::filelock() : is_locking_(false) {}

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
      std::cerr << "Waring: file \"" << file_ << "\" is locked.  Still trying.\n";
#if defined(ALPS_HAVE_UNISTD_H)
      sleep(1);    // sleep 1 Sec
#elif defined(ALPS_HAVE_WINDOWS_H)
      Sleep(1000); // sleep 1000 mSec
#else
# error "sleep not found"
#endif
      //sleep(1);
    }
#if defined(ALPS_HAVE_WINDOWS_H)
    int fd = _open(lock_.string().c_str(), O_WRONLY | O_CREAT | O_EXCL , _S_IWRITE);
#else
    int fd = open(lock_.string().c_str(), O_WRONLY | O_CREAT | O_EXCL , S_IWRITE);
#endif

    if (fd > 0) {
      is_locking_ = true;
#if defined(ALPS_HAVE_WINDOWS_H)
      _close(fd);
#else
      close(fd);
#endif	  
      break;
    }
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
