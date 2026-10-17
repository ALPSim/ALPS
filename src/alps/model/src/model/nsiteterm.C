/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// XML reading/writing and operator substitution for n-site terms; mirrors
// bondterm.C.

#include <alps/model/nsiteoperator.h>
#include <alps/model/operatorsubstitution.h>
#include <boost/algorithm/string/split.hpp>
#include <algorithm>
#include <limits>
#include <sstream>
#include <boost/algorithm/string/classification.hpp>

#ifndef ALPS_WITHOUT_XML

void alps::NSiteOperator::read_xml(const XMLTag& intag, std::istream& is)
{
  XMLTag tag(intag);
  name_ = tag.attributes["name"];
  sites_.clear();
  std::string s = tag.attributes["sites"];
  std::vector<std::string> words;
  boost::algorithm::split(words, s, boost::algorithm::is_any_of(" \t\n,"), boost::algorithm::token_compress_on);
  for (std::size_t i = 0; i < words.size(); ++i)
    if (!words[i].empty()) {
      if (std::find(sites_.begin(), sites_.end(), words[i]) != sites_.end())
        boost::throw_exception(std::runtime_error("site name " + words[i] + " appears twice in <" + intag.name + " sites=\"" + s + "\">"));
      sites_.push_back(words[i]);
    }
  if (sites_.empty())
    boost::throw_exception(std::runtime_error("<" + intag.name + "> needs a sites attribute naming its sites, e.g. sites=\"i j k l\""));
  if (tag.type != XMLTag::SINGLE) {
    term_ = parse_content(is);
    while (true) {
      tag = parse_tag(is, false);
      if (tag.name == "/" + intag.name)
        return;
      if (tag.name == "PARAMETER") {
        parms_[tag.attributes["name"]] = tag.attributes["default"];
        if (tag.type != XMLTag::SINGLE) {
          tag = parse_tag(is);
          if (tag.name != "/PARAMETER")
            boost::throw_exception(std::runtime_error("End tag </PARAMETER> missing while parsing " + name() + " Hamiltonian"));
        }
      }
      else if (tag.type != XMLTag::COMMENT)
        boost::throw_exception(std::runtime_error("Illegal tag <" + tag.name + "> in <" + intag.name + "> element"));
      std::string next_part = parse_content(is);
      if (!term_.empty() && !next_part.empty())
        term_ += " ";
      term_ += next_part;
    }
  }
}

void alps::NSiteOperator::substitute_operators(const ModelLibrary& m, const Parameters& p)
{
  OperatorSubstitution<std::complex<double> > subs(m, p, sites());
  Expression e(term());
  e.partial_evaluate(subs);
  e.simplify();
  // Print with full precision: lexical_cast uses the stream default of 6
  // significant digits, which would truncate numeric coefficients in the term
  // (e.g. two-body integrals written as literals).
  std::ostringstream os;
  os.precision(std::numeric_limits<double>::max_digits10);
  os << e;
  term_ = os.str();
}

std::string alps::NSiteOperator::sites_string() const
{
  std::string s;
  for (std::size_t i = 0; i < sites_.size(); ++i)
    s += (i ? " " : "") + sites_[i];
  return s;
}

alps::NSiteTermDescriptor::NSiteTermDescriptor(const XMLTag& intag, std::istream& is)
{
  XMLTag tag(intag);
  type_ = tag.attributes["type"] == "" ? -1 : boost::lexical_cast<int, std::string>(tag.attributes["type"]);
  read_xml(intag, is);
}

void alps::NSiteOperator::write_xml(oxstream& os) const
{
  os << start_tag("NSITEOPERATOR");
  if (!name().empty())
    os << attribute("name", name());
  os << attribute("sites", sites_string());
  for (Parameters::const_iterator it = parms().begin(); it != parms().end(); ++it)
    os << start_tag("PARAMETER") << attribute("name", it->key())
       << attribute("default", it->value()) << end_tag("PARAMETER");
  os << term();
  os << end_tag("NSITEOPERATOR");
}

void alps::NSiteTermDescriptor::write_xml(oxstream& os) const
{
  os << start_tag("NSITETERM");
  if (type_ >= 0)
    os << attribute("type", type_);
  os << attribute("sites", sites_string());
  for (Parameters::const_iterator it = parms().begin(); it != parms().end(); ++it)
    os << start_tag("PARAMETER") << attribute("name", it->key())
       << attribute("default", it->value()) << end_tag("PARAMETER");
  os << term();
  os << end_tag("NSITETERM");
}

#endif
