/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_LOCATION_IMPL_H
#define LOOPER_LOCATION_IMPL_H

#include "location.h"
#include <alps/osiris.h>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace looper {

class location {
public:
  location(int pos = 0, bool is_bond = true) : loc_(pos << 1 | (is_bond ? 1 : 0)) {}
  int pos() const { return loc_ >> 1; }
  bool is_bond() const { return loc_ & 1; }
  bool is_site() const { return !is_bond(); }
  bool operator==(const location& rhs) const { return loc_ == rhs.loc_; }
  bool operator!=(const location& rhs) const { return !operator==(rhs); }
  void save(alps::ODump& dump) const { dump << loc_; }
  void load(alps::IDump& dump) { dump >> loc_; }
  static location bond_location(int pos) { return location(pos, true); }
  static location site_location(int pos) { return location(pos, false); }
private:
  int loc_;
};

inline int pos(const location& loc) { return loc.pos(); }
inline bool is_bond(const location& loc) { return loc.is_bond(); }
inline bool is_site(const location& loc) { return loc.is_site(); }

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dump, const looper::location& loc) {
  loc.save(dump); return dump;
}

inline alps::IDump& operator>>(alps::IDump& dump, looper::location& loc) {
  loc.load(dump); return dump;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif


namespace looper {

// optimized version for models with bond terms only

class location_bond {
public:
  location_bond(int pos = 0, bool is_bond = true) : loc_(pos) {
    if (!is_bond) boost::throw_exception(std::invalid_argument("location_bond"));
  }
  int pos() const { return loc_; }
  static bool is_bond() { return true; }
  static bool is_site() { return false; }
  bool operator==(const location_bond& rhs) const { return loc_ == rhs.loc_; }
  bool operator!=(const location_bond& rhs) const { return !operator==(rhs); }
  static location_bond bond_location(int pos) { return location_bond(pos); }
  static location_bond site_location(int) {
    boost::throw_exception(std::logic_error("location_bond"));
    return location_bond();
  }
  void save(alps::ODump& dump) const { dump << loc_; }
  void load(alps::IDump& dump) { dump >> loc_; }
private:
  int loc_;
};

inline int pos(const location_bond& loc) { return loc.pos(); }
inline bool is_bond(const location_bond&) { return true; }
inline bool is_site(const location_bond&) { return false; }

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dump, const looper::location_bond& loc) {
  loc.save(dump); return dump;
}

inline alps::IDump& operator>>(alps::IDump& dump, looper::location_bond& loc) {
  loc.load(dump); return dump;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif


namespace looper {

// for long-range models

class location_longrange {
public:
  location_longrange(int pos = 0) : source_(pos << 1) {}
  location_longrange(int source, int target) : source_(source << 1 | 1), target_(target) {}
  int pos() const { return source_ >> 1; }
  int source() const { return source_ >> 1; }
  int target() const { return target_; }
  bool is_bond() const { return source_ & 1; }
  bool is_site() const { return !is_bond(); }
  bool operator==(const location_longrange& rhs) const {
    return source_ == rhs.source_ && target_ == rhs.target_;
  }
  bool operator!=(const location_longrange& rhs) const { return !operator==(rhs); }
  void save(alps::ODump& dump) const { dump << source_ << target_; }
  void load(alps::IDump& dump) { dump >> source_ >> target_; }
  static location_longrange bond_location(int source, int target) {
    return location_longrange(source, target);
  }
  static location_longrange site_location(int pos) { return location_longrange(pos); }
private:
  int source_;
  int target_;
};

inline int pos(const location_longrange& loc) { return loc.pos(); }
inline int source(const location_longrange& loc) { return loc.source(); }
inline int target(const location_longrange& loc) { return loc.target(); }
inline bool is_bond(const location_longrange&) { return true; }
inline bool is_site(const location_longrange&) { return false; }

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dump, const looper::location_longrange& loc) {
  loc.save(dump); return dump;
}

inline alps::IDump& operator>>(alps::IDump& dump, looper::location_longrange& loc) {
  loc.load(dump); return dump;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_LOCATION_IMPL_H
