/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2005 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/model/globaloperator.h>
#include <boost/foreach.hpp>

void alps::GlobalOperator::substitute_operators(const ModelLibrary& m, const Parameters& p)
{
  for (std::vector<SiteTermDescriptor>::iterator it=siteterms_.begin();it!=siteterms_.end();++it)
    it->substitute_operators(m,p);
  for (std::vector<BondTermDescriptor>::iterator it=bondterms_.begin();it!=bondterms_.end();++it)
    it->substitute_operators(m,p);
  for (std::vector<NSiteTermDescriptor>::iterator it=nsiteterms_.begin();it!=nsiteterms_.end();++it)
    it->substitute_operators(m,p);
  default_siteterm_.substitute_operators(m,p);
  default_bondterm_.substitute_operators(m,p);
  if (!default_nsiteterm_.term().empty())
    default_nsiteterm_.substitute_operators(m,p);
}

alps::SiteOperator alps::GlobalOperator::site_term(unsigned int type) const
{
  for (std::vector<SiteTermDescriptor>::const_iterator it =siteterms_.begin();it!=siteterms_.end();++it)
    if (it->match_type(type))
      return it->site_operator();
  return default_siteterm_.get(type).site_operator();
}

alps::BondOperator alps::GlobalOperator::bond_term(unsigned int type) const
{
  for (std::vector<BondTermDescriptor>::const_iterator it =bondterms_.begin();it!=bondterms_.end();++it)
    if (it->match_type(type))
      return it->bond_operator();
  return default_bondterm_.get(type).bond_operator();
}

alps::NSiteOperator alps::GlobalOperator::nsite_term(unsigned int type) const
{
  for (std::vector<NSiteTermDescriptor>::const_iterator it =nsiteterms_.begin();it!=nsiteterms_.end();++it)
    if (it->match_type(type))
      return it->nsite_operator();
  if (default_nsiteterm_.term().empty())
    return NSiteOperator();
  return default_nsiteterm_.get(type).nsite_operator();
}

boost::optional<alps::Parameters> alps::GlobalOperator::create_nsite_term(unsigned int type)
{
  // A typed term's own <PARAMETER> defaults must reach the Hamiltonian
  // parameters as well, or a parameter declared only there cannot be evaluated.
  for (std::vector<NSiteTermDescriptor>::const_iterator it =nsiteterms_.begin();it!=nsiteterms_.end();++it)
    if (it->match_type(type))
      return it->parms().empty() ? boost::optional<Parameters>() : boost::optional<Parameters>(it->parms());
  if (default_nsiteterm_.term().empty())
    return boost::optional<Parameters>();
  nsiteterms_.push_back(default_nsiteterm_.get(type));
  return default_nsiteterm_.parms(type);
}

alps::Parameters alps::GlobalOperator::create_nsite_terms(std::set<unsigned int> const& types)
{
  Parameters p;
  BOOST_FOREACH(unsigned int const& t, types) {
    boost::optional<Parameters> newp = create_nsite_term(t);
    if (newp)
      p << newp.get();
  }
  return p;
}

boost::optional<alps::Parameters> alps::GlobalOperator::create_site_term(unsigned int type)
{
  for (std::vector<SiteTermDescriptor>::const_iterator it =siteterms_.begin();it!=siteterms_.end();++it)
    if (it->match_type(type))
      return boost::optional<Parameters>();
  siteterms_.push_back(default_siteterm_.get(type));
  return default_siteterm_.parms(type);
}

boost::optional<alps::Parameters> alps::GlobalOperator::create_bond_term(unsigned int type)
{
  for (std::vector<BondTermDescriptor>::const_iterator it =bondterms_.begin();it!=bondterms_.end();++it)
    if (it->match_type(type))
      return boost::optional<Parameters>();
  bondterms_.push_back(default_bondterm_.get(type));
  return default_bondterm_.parms(type);
}

alps::Parameters alps::GlobalOperator::create_site_terms(std::set<unsigned int> const& types)
{
  Parameters p;
  BOOST_FOREACH(unsigned int const& t, types) {
    boost::optional<Parameters> newp = create_site_term(t);
    if (newp)
      p << newp.get();
  }
  return p;
}

alps::Parameters alps::GlobalOperator::create_bond_terms(std::set<unsigned int> const& types)
{
  Parameters p;
  BOOST_FOREACH(unsigned int const& t, types) {
    boost::optional<Parameters> newp = create_bond_term(t);
    if (newp)
      p << newp.get();
  }
  return p;
}



#ifndef ALPS_WITHOUT_XML
alps::GlobalOperator::GlobalOperator(const XMLTag& intag, std::istream& is)
{
  XMLTag tag(intag);
  name_=tag.attributes["name"];
  if (tag.type!=XMLTag::SINGLE) {
    tag = parse_tag(is);
    tag = read_xml(tag,is);
    if (tag.name != "/"+intag.name)
      boost::throw_exception(std::runtime_error("Encountered illegal tag <"+tag.name+"> in element <"+intag.name+">"));
  }
}

    
alps::XMLTag alps::GlobalOperator::read_xml(const XMLTag& intag, std::istream& is) {
  XMLTag tag(intag);
  while (true) {
    if (tag.name=="SITETERM") {
      if (tag.attributes["type"]=="")
        default_siteterm_=DefaultSiteTermDescriptor(tag,is);
      else
        siteterms_.push_back(SiteTermDescriptor(tag,is));
    }
    else if (tag.name=="BONDTERM") {
      if (tag.attributes["type"]=="")
        default_bondterm_=DefaultBondTermDescriptor(tag,is);
      else
        bondterms_.push_back(BondTermDescriptor(tag,is));
    }
    else if (tag.name=="NSITETERM") {
      if (tag.attributes["type"]=="")
        default_nsiteterm_=DefaultNSiteTermDescriptor(tag,is);
      else
        nsiteterms_.push_back(NSiteTermDescriptor(tag,is));
    }
    else
      return tag;
    tag=parse_tag(is);
  }
}

void alps::GlobalOperator::write_xml(oxstream& os) const
{
  os << start_tag("GLOBALOPERATOR") << attribute("name", name());
  write_operators_xml(os);
  os << end_tag("GLOBALOPERATOR");
}

void alps::GlobalOperator::write_operators_xml(oxstream& os) const
{
  if (!default_siteterm_.term().empty())
    default_siteterm_.write_xml(os);
  for (std::vector<SiteTermDescriptor>::const_iterator it=siteterms_.begin();it!=siteterms_.end();++it)
    it->write_xml(os);
  if (!default_bondterm_.term().empty())
    default_bondterm_.write_xml(os);
  for (std::vector<BondTermDescriptor>::const_iterator it=bondterms_.begin();it!=bondterms_.end();++it)
    it->write_xml(os);
  if (!default_nsiteterm_.term().empty())
    default_nsiteterm_.write_xml(os);
  for (std::vector<NSiteTermDescriptor>::const_iterator it=nsiteterms_.begin();it!=nsiteterms_.end();++it)
    it->write_xml(os);
}

#endif
