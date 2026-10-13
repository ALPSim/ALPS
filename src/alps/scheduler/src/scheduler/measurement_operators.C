/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
