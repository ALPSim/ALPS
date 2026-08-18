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
