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

#ifndef ALPS_PARAMETER_PARAMETERLIST_P_H
#define ALPS_PARAMETER_PARAMETERLIST_P_H

#include "parameterlist.h"
#include "parameters_p.h"
#include <boost/classic_spirit.hpp>

namespace bs = boost::spirit;

namespace alps {

// parser for alps::ParameterList

struct ALPS_DECL ParameterListParser : public bs::grammar<ParameterListParser> {

  template<typename ScannerT>
  struct definition {

    bs::rule<ScannerT> parameterlist;

    definition(ParameterListParser const& self) {
      self.stop = 0;
      parameterlist =
        +( self.global_p
            | ( bs::ch_p('{') >> *bs::eol_p >> bs::ch_p('}') >> *bs::eol_p
              )[bs::push_back_a(self.plist, self.global)]
            | ( bs::ch_p('{')[bs::assign_a(self.local, self.global)] >> *bs::eol_p
                >> self.local_p >> bs::ch_p('}') >> *bs::eol_p
              )[bs::push_back_a(self.plist, self.local)]
            | ( bs::str_p("#clear") >> !bs::ch_p(";") >> *bs::eol_p )[bs::clear_a(self.global)]
         )
        >> !( bs::str_p("#stop") >> !bs::ch_p(";") >> *bs::eol_p )[bs::increment_a(self.stop)];
    }

    bs::rule<ScannerT> const& start() const {
      return parameterlist;
    }
  };

  ParameterListParser(ParameterList& p) :
    plist(p), global_p(global), local_p(local), stop(0) {}

  ParameterList& plist;
  mutable Parameters global, local;
  ParametersParser global_p, local_p;
  mutable unsigned stop;
};

/// \brief Implementation handler of the ALPS XML parser for the ParameterList class
class ALPS_DECL ParameterListXMLHandler : public CompositeXMLHandler
{
public:
  ParameterListXMLHandler(ParameterList& list);

protected:
  void start_child(const std::string& name,
                   const XMLAttributes& attributes,
                   xml::tag_type type);
  void end_child(const std::string& name, xml::tag_type type);

private:
  ParameterList& list_;
  Parameter parameter_;
  Parameters default_, current_;
  ParameterXMLHandler parameter_handler_;
  ParametersXMLHandler current_handler_;
};

} // end namespace alps

#endif // ALPS_PARAMETER_PARAMETERLIST_P_H
