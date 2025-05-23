/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2008 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

#include "parameterlist.h"
#include "parameterlist_p.h"
#include <alps/parser/parser.h>
#include <boost/foreach.hpp>
#include <boost/throw_exception.hpp>
#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace bs = boost::spirit;

namespace alps {

void ParameterList::parse(std::istream& is, bool replace_env) {
  std::deque<char> buff;
  std::copy(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>(),
    std::back_inserter(buff));
  ParameterListParser plist_p(*this);
  bs::parse_info<std::deque<char>::iterator> info = bs::parse(
    buff.begin(), buff.end(),
    plist_p >> bs::end_p,
    bs::blank_p | bs::comment_p("//") | bs::comment_p("/*", "*/"));

  /* // ST 2006.10.06: following in-situ version does not work with Intel C++ on IA64
  typedef bs::multi_pass<std::istreambuf_iterator<char> > iterator_t;
  iterator_t first = bs::make_multi_pass(std::istreambuf_iterator<char>(is));
  iterator_t last = bs::make_multi_pass(std::istreambuf_iterator<char>());
  bs::parse_info<bs::multi_pass<std::istreambuf_iterator<char> > > info = bs::parse(
    first, last,
    ...
  */

  if (!(info.full || plist_p.stop)) {
    std::deque<char>::iterator itr = info.stop;
    std::string err = "parameter parse error at \"";
    for (int i = 0; itr != buff.end() && i < 32; ++itr, ++i)
      err += (*itr != '\n' ? *itr : ' ');
    boost::throw_exception(std::runtime_error(err + "\""));
  }
  if (replace_env) replace_envvar();
}

void ParameterList::replace_envvar() {
  for (iterator itr = this->begin(); itr != this->end(); ++itr) itr->replace_envvar();
  // BOOST_FOREACH(Parameters& p, *this) p.replace_envvar();
}


//
// XML support
//

ParameterListXMLHandler::ParameterListXMLHandler(ParameterList& list)
  : CompositeXMLHandler("PARAMETERLIST"), list_(list),
    parameter_(), default_(), current_(), parameter_handler_(parameter_),
    current_handler_(current_) {
  add_handler(parameter_handler_);
  add_handler(current_handler_);
}

void ParameterListXMLHandler::start_child(const std::string& name,
  const XMLAttributes& /* attributes */, xml::tag_type type) {
  if (type == xml::element) {
    if (name == "PARAMETER") {
      parameter_ = Parameter();
    } else if (name == "PARAMETERS") {
      current_ = default_;
    }
  }
}
void ParameterListXMLHandler::end_child(const std::string& name,
                                        xml::tag_type type) {
  if (type == xml::element) {
    if (name == "PARAMETER") {
      default_[parameter_.key()] = parameter_.value();
    } else if (name == "PARAMETERS") {
      list_.push_back(current_);
    }
  }
}

}
