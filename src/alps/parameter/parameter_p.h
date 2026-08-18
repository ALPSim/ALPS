/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2006-2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_PARAMETER_PARAMTER_P_H
#define ALPS_PARAMETER_PARAMTER_P_H

#include "parameter.h"
#include <alps/xml.h>
#include <boost/classic_spirit.hpp>

namespace bs = boost::spirit;

namespace alps {

/// \brief Text-form parser for the Parameter class
struct ALPS_DECL ParameterParser : public bs::grammar<ParameterParser> {

  template<typename ScannerT>
  struct definition {

    bs::rule<ScannerT> parameter;
    bs::rule<ScannerT> key;
    bs::rule<ScannerT> value;

    definition(ParameterParser const& self) {
      parameter = key >> '=' >> value;
      key =
        ( bs::alpha_p
          >> *( bs::alnum_p | '_' | '\'' | '#' | bs::confix_p('[', *bs::print_p, ']') )
        )[bs::assign_a(self.param.key())],
      value =
        bs::lexeme_d
          [ bs::confix_p('"', (*bs::print_p)[bs::assign_a(self.param.value())], '"')
          | bs::confix_p('\'', (*bs::print_p)[bs::assign_a(self.param.value())], '\'')
          /* | bs::confix_p('[', (*bs::print_p)[bs::assign_a(self.param.value())], ']') */
          | ( *( bs::alnum_p | '#' | bs::range_p('\'', '+') | bs::range_p('-', '/')
               | bs::range_p('^', '_')
               )
              % ( bs::ch_p('$') >> bs::confix_p('{', *bs::graph_p, '}') )
            )[bs::assign_a(self.param.value())]
          ];
    }

    bs::rule<ScannerT> const& start() const {
      return parameter;
    }
  };

  ParameterParser(Parameter& p) : param(p) {}

  Parameter& param;
};

/// \brief ALPS XML handler for the Parameter class
class ALPS_DECL ParameterXMLHandler : public XMLHandlerBase {

public:
  ParameterXMLHandler(Parameter& p);

  void start_element(const std::string& name, const XMLAttributes& attributes, xml::tag_type type);
  void end_element(const std::string& name, xml::tag_type type);
  void text(const std::string& text);

private:
  Parameter& parameter_;
};

} // end namespace alps

#endif // ALPS_PARAMETER_PARAMETER_P_H
