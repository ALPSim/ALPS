/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
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

/* $Id$ */

#ifndef ALPS_PARSER_XMLATTRIBUTES_H
#define ALPS_PARSER_XMLATTRIBUTES_H

// for MSVC
#if defined(_MSC_VER)
# pragma warning(disable:4251)
#endif

#include <alps/config.h>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace alps {

class ALPS_DECL XMLAttribute
{
public:
  typedef std::string key_type;
  typedef std::string value_type;

  XMLAttribute(const XMLAttribute& attr)
    : key_(attr.key_), value_(attr.value_) {}
  XMLAttribute(const key_type& k) : key_(k), value_() {}
  XMLAttribute(const key_type& k, const value_type& v) : key_(k), value_(v) {}
  XMLAttribute(const key_type& k, const char * v) : key_(k), value_(v) {}
  template<class T>
  XMLAttribute(const key_type& k, const T& v)
    : key_(k), value_(boost::lexical_cast<std::string, T>(v)) {}

  key_type& key() { return key_; }
  const key_type& key() const { return key_; }
  value_type& value() { return value_; }
  const value_type& value() const { return value_; }

private:
  key_type key_;
  value_type value_;
};

class ALPS_DECL XMLAttributes
{
public:
  typedef XMLAttribute::key_type    key_type;
  typedef XMLAttribute::value_type  value_type;

  typedef std::vector<XMLAttribute> list_type;
  typedef list_type::size_type      size_type;
  typedef list_type::iterator       iterator;
  typedef list_type::const_iterator const_iterator;

public:
  XMLAttributes() {}
  XMLAttributes(const std::string& str);

  void clear() { list_.clear(); map_.clear(); }
  size_type size() const { return list_.size(); }

  bool defined(const key_type& k) const {
    return (map_.find(k) != map_.end());
  }

  // accessing elements by key
  value_type& operator[](const key_type& k) {
    if (defined(k)) {
      return list_[map_.find(k)->second].value();
    } else {
      push_back(k, value_type());
      return list_.rbegin()->value();
    }
  }
  const value_type& operator[](const key_type& k) const {
    if (!defined(k))
      boost::throw_exception(std::runtime_error("attribute not defined"));
    return list_[map_.find(k)->second].value();
  }
  value_type value_or_default(const key_type& k, const value_type& v) const {
    return defined(k) ? (*this)[k] : v;
  }

  iterator begin() { return list_.begin(); }
  const_iterator begin() const { return list_.begin(); }
  iterator end() { return list_.end(); }
  const_iterator end() const { return list_.end(); }

  void push_back(const XMLAttribute& attr);
  void push_back(const key_type& k, const value_type& v) {
    push_back(XMLAttribute(k, v));
  }
  XMLAttributes& operator<<(const XMLAttribute& a) {
    (*this)[a.key()] = a.value();
    return *this;
  }
  XMLAttributes& operator<<(const XMLAttributes& attr) {
    for (const_iterator itr = attr.begin(); itr != attr.end(); ++itr)
      (*this) << *itr;
    return *this;
  }

private:
  typedef std::map<key_type, size_type> map_type;

  list_type list_;
  map_type map_;
};

} // namespace alps

#endif // ALPS_PARSER_XMLATTRIBUTE_H
