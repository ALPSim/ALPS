/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2006-2009 by Synge Todo <wistaria@comp-phys.org>
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
