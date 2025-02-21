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

#ifndef ALPS_PARAMETER_PARAMETERS_P_H
#define ALPS_PARAMETER_PARAMETERS_P_H

#include "parameters.h"
#include "parameter_p.h"
#include <boost/classic_spirit.hpp>

namespace bs = boost::spirit;

namespace alps {

/// \brief Text-form parser for the Parameters class
struct ALPS_DECL ParametersParser : public bs::grammar<ParametersParser> {

  template<typename ScannerT>
  struct definition {

    bs::rule<ScannerT> parameters;

    definition(ParametersParser const& self) {
      parameters =
        *bs::eol_p
        >> self.parameter_p[bs::assign_key_a(self.params, self.param.value(), self.param.key())]
           % ( ( bs::ch_p(";") | bs::ch_p(",") | bs::eol_p ) >> *bs::eol_p )
        >> !bs::ch_p(";") >> *bs::eol_p;
    }

    bs::rule<ScannerT> const& start() const {
      return parameters;
    }
  };

  ParametersParser(Parameters& p) : params(p), parameter_p(param) {}

  Parameters& params;
  mutable Parameter param;
  ParameterParser parameter_p;
};

//
// XML support
//

/// \brief ALPS XML handler for the Parameters class
class ALPS_DECL ParametersXMLHandler : public CompositeXMLHandler
{
public:
  ParametersXMLHandler(Parameters& p);

protected:
  void start_child(const std::string& name,
                   const XMLAttributes& attributes,
                   xml::tag_type type);
  void end_child(const std::string& name, xml::tag_type type);

private:
  Parameters& parameters_;
  Parameter parameter_;
  ParameterXMLHandler handler_;
};

} // namespace alps

#endif // ALPS_PARAMETER_PARAMETERS_P_H
