/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parser/xmlattributes.h>
#include <boost/classic_spirit.hpp>

namespace alps {

namespace {

struct assign_string {
  assign_string(const std::string& name)
    : ptr_(const_cast<std::string*>(&name)) {}
  void operator()(const char* first, const char* last) const
  { *ptr_ = std::string(first, last); }
  std::string* ptr_;
};

struct append_attr {
  append_attr(const XMLAttributes& attr, const std::string& name,
              const std::string& value)
    : ptr_(const_cast<XMLAttributes*>(&attr)), name_(name), value_(value) {}
  void operator()(const char*, const char*) const
  { ptr_->push_back(XMLAttribute(name_, value_)); }
  XMLAttributes* ptr_;
  const std::string& name_;
  const std::string& value_;
};

struct attr_parser : public boost::spirit::grammar<attr_parser> {
  attr_parser(XMLAttributes& a) : attr(a) {}

  template<typename ScannerT>
  struct definition {
    boost::spirit::rule<ScannerT> name_p, value_p, attribute_p, attributes_p;
    definition(const attr_parser& self)
    {
      using boost::spirit::alpha_p;
      using boost::spirit::alnum_p;
      using boost::spirit::anychar_p;
      using boost::spirit::ch_p;
      name_p = ((alpha_p | ch_p('_')) >> *(alnum_p | ch_p('_')))[assign_string(self.name)];
      value_p = ch_p('\"') >> (*(anychar_p - ch_p('\"')))[assign_string(self.value)] >> ch_p('\"');
      attributes_p = *(name_p >> ch_p('=') >> value_p)[append_attr(self.attr, self.name, self.value)];
    }
    const boost::spirit::rule<ScannerT>& start() const { return attributes_p; }
  };

  XMLAttributes& attr;
  std::string name, value;
};

}

XMLAttributes::XMLAttributes(const std::string& str)
{
  if (!boost::spirit::parse(str.c_str(), attr_parser(*this),
    boost::spirit::space_p).full)
      boost::throw_exception(std::runtime_error("parse failed"));
}

void XMLAttributes::push_back(const XMLAttribute& attr)
{
  if (defined(attr.key()))
    boost::throw_exception(std::runtime_error("duplicated attribute " +
      attr.key()));
  map_[attr.key()] = list_.size();
  list_.push_back(attr);
}

} // namespace alps
