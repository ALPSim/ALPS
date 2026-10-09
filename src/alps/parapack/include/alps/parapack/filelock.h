/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef PARAPACK_FILELOCK_H
#define PARAPACK_FILELOCK_H

#include <alps/config.h>

#include <boost/filesystem/path.hpp>

namespace alps {

class filelock {
public:
  filelock();
  filelock(boost::filesystem::path const& file, bool lock_now = false, int wait = -1,
    bool auto_release = true);
  ~filelock();

  void set_file(boost::filesystem::path const& file);

  void lock(int wait = -1); // wait = -1 means waiting forever
  void release();

  bool locking() const;
  bool locked() const;

private:
  std::string file_;
  boost::filesystem::path lock_;
  bool auto_release_;
  bool is_locking_;
};

} // end namespace alps

#endif // PARAPACK_FILELOCK_H
