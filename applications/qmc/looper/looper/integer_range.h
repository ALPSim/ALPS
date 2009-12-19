/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef LOOPER_INTEGER_RANGE_H
#define LOOPER_INTEGER_RANGE_H

#include <alps/expression.h>
#include <alps/osiris.h>
#include <boost/call_traits.hpp>
#include <boost/spirit/core.hpp>
#include <boost/throw_exception.hpp>
#include <iosfwd>
#include <limits>
#include <stdexcept>
#include <string>

namespace looper {

//
// class tempalte integer_range
//

template<class T>
class integer_range {
public:
  typedef T value_type;
  typedef typename boost::call_traits<value_type>::param_type param_type;

  integer_range() : mi_(0), ma_(0) {}
  explicit integer_range(param_type v) : mi_(v), ma_(v) {}
  explicit integer_range(param_type vmin, param_type vmax) : mi_(vmin), ma_(vmax) {
    if (mi_ > ma_)
      boost::throw_exception(std::runtime_error("integer_range: range error"));
  }
  integer_range(integer_range const& r) : mi_(r.mi_), ma_(r.ma_) {}
  integer_range(std::string const& str) :
    mi_(std::numeric_limits<value_type>::min()), ma_(std::numeric_limits<value_type>::max()) {
    init(str, alps::Parameters());
  }
  integer_range(std::string const& str, alps::Parameters const& p) :
    mi_(std::numeric_limits<value_type>::min()), ma_(std::numeric_limits<value_type>::max()) {
    init(str, p);
  }
  integer_range(std::string const& str, alps::Parameters const& p, param_type def_mi,
    param_type def_ma) : mi_(def_mi), ma_(def_ma) {
    init(str, p);
  }

  integer_range& operator=(param_type v) {
    mi_ = ma_ = v;
    return *this;
  }

  template<typename U>
  integer_range& operator*=(U x) {
    if (x > U(0)) {
      if (x > U(1) && mi_ > static_cast<value_type>(std::numeric_limits<value_type>::max() / x))
        mi_ = std::numeric_limits<value_type>::max();
      else
        mi_ = static_cast<value_type>(x * mi_);
      if (x > U(1) && ma_ > static_cast<value_type>(std::numeric_limits<value_type>::max() / x))
        ma_ = std::numeric_limits<value_type>::max();
      else
        ma_ = static_cast<value_type>(x * ma_);
    } else {
      boost::throw_exception(std::runtime_error("integer_range: multiplier should be positive"));
    }
    return *this;
  }

  void include(param_type v) {
    if (v < mi_) mi_ = v;
    if (v > ma_) ma_ = v;
  }

  value_type min() const { return mi_; }
  value_type max() const { return ma_; }
  value_type size() const { return ma_ - mi_ + 1; }
  bool is_included(param_type v) const { return (v >= min()) && (v <= max()); }

  void save(alps::ODump& dp) const { dp << mi_ << ma_; }
  void load(alps::IDump& dp) { dp >> mi_ >> ma_; }

protected:
  void init(std::string const& str, alps::Parameters const& p) {
    using namespace boost::spirit;
    std::string mi_str, ma_str;
    if (!parse(
      str.c_str(),
          ( ch_p('[')
            >> (*(anychar_p-'['-':'-']'))[assign_a(mi_str)]
            >> ':'
            >> (*(anychar_p-'['-':'-']'))[assign_a(ma_str)]
            >> ']' )
        | ( ch_p('[')
            >> (+(anychar_p-'['-':'-']'))[assign_a(mi_str)][assign_a(ma_str)]
            >> ']' )
        | (+(anychar_p-'['-':'-']'))[assign_a(mi_str)][assign_a(ma_str)]
        >> end_p,
      space_p).full)
      boost::throw_exception(std::runtime_error("integer_range: parse error: " + str));
    if (mi_str.empty()) {
      mi_ = std::numeric_limits<value_type>::min();
    } else {
      double v = alps::evaluate(mi_str, p);
      if (v < std::numeric_limits<value_type>::min() ||
          v > std::numeric_limits<value_type>::max())
        boost::throw_exception(std::runtime_error("integer_range: range error"));
      mi_ = static_cast<value_type>(v);
      // Note: boost::numeric_cast<value_type>(v) does not work for intel C++ for x86_64
    }
    if (ma_str.empty()) {
      ma_ = std::numeric_limits<value_type>::max();
    } else {
      double v = alps::evaluate(ma_str, p);
      if (v < std::numeric_limits<value_type>::min() ||
          v > std::numeric_limits<value_type>::max())
        boost::throw_exception(std::runtime_error("integer_range: range error"));
      ma_ = static_cast<value_type>(v);
      // Note: boost::numeric_cast<value_type>(v) does not work for intel C++ for x86_64
    }
    if (mi_ > ma_)
      boost::throw_exception(std::runtime_error("integer_range: range error"));
  }

private:
  value_type mi_, ma_;
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<typename T, typename U>
looper::integer_range<T> operator*(looper::integer_range<T> const& t, U x) {
  looper::integer_range<T> result = t;
  result *= x;
  return result;
}

template<typename T, typename U>
looper::integer_range<T> operator*(U x, looper::integer_range<T> const& t) {
  return t * x;
}

template<class T>
std::ostream& operator<<(std::ostream& os, integer_range<T> const& ir) {
  os << '[' << ir.min() << ':' << ir.max() << ']';
  return os;
}

template<class T>
alps::ODump& operator<<(alps::ODump& dp, looper::integer_range<T> const& ir) {
  ir.save(dp);
  return dp;
}

template<class T>
alps::IDump& operator>>(alps::IDump& dp, looper::integer_range<T>& ir) {
  ir.load(dp);
  return dp;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_INTEGER_RANGE_H
