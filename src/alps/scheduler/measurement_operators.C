/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@comp-phys.org>,
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

#include <alps/scheduler/measurement_operators.h>

// some file (probably a python header) defines a tolower macro ...
#undef tolower
#undef toupper

#include <boost/regex.hpp> 

alps::MeasurementOperators::MeasurementOperators (Parameters const& parms)
{
  boost::regex expression("^MEASURE_AVERAGE\\[(.*)]$");
  boost::smatch what;
  for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
    std::string lhs = it->key();
    if (boost::regex_match(lhs, what, expression))
      average_expressions[what.str(1)]=it->value();
  }

  expression = boost::regex("^MEASURE_LOCAL\\[(.*)]$");
  for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
    std::string lhs = it->key();
    if (boost::regex_match(lhs, what, expression))
      local_expressions[what.str(1)]=it->value();
  }

  expression = boost::regex("^MEASURE_CORRELATIONS\\[(.*)]$");
  for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
    std::string lhs = it->key();
    if (boost::regex_match(lhs, what, expression)) {
      std::string key = what.str(1);
      std::string value = it->value();
      boost::regex expression2("^(.*):(.*)$");
      if (boost::regex_match(value, what, expression2))
        correlation_expressions[key] = std::make_pair(what.str(1), what.str(2));
      else
        correlation_expressions[key] = std::make_pair(value, value);
    }
  }

  expression = boost::regex("^MEASURE_STRUCTURE_FACTOR\\[(.*)]$");
  for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
    std::string lhs = it->key();
    if (boost::regex_match(lhs, what, expression)) {
      std::string key = what.str(1);
      std::string value = it->value();
      boost::regex expression2("^(.*):(.*)$");
      if (boost::regex_match(value, what, expression2))
        structurefactor_expressions[key] = std::make_pair(what.str(1), what.str(2));
      else
        structurefactor_expressions[key]=std::make_pair(value, value);
    }
  }
}
