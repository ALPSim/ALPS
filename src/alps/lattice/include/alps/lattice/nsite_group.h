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

#ifndef ALPS_LATTICE_NSITE_GROUP_H
#define ALPS_LATTICE_NSITE_GROUP_H

// n-site groups of a graph.
//
// An n-site group is an ordered tuple of distinct sites plus a type. The order is
// the order of the roles in an <NSITETERM sites="i j k l"> of the model.
// N-site groups are kept in a list next to the Boost graph; the graph itself
// still contains only sites and 2-site bonds.

#include <alps/lattice/graphproperties.h>
#include <alps/parser/parser.h>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace alps {

struct nsite_group
{
  nsite_group() : type(0) {}
  type_type type;
  std::vector<std::size_t> sites;   // 0-based vertex indices, in role order
};

typedef std::vector<nsite_group> nsite_group_list;

// An n-site group as written in a <UNITCELL>: each site is a unit-cell vertex
// (0-based) plus a cell offset, exactly like the SOURCE/TARGET of an EDGE.
struct unitcell_nsite_group
{
  unitcell_nsite_group() : type(0) {}
  type_type type;
  std::vector<std::pair<std::size_t, std::vector<int> > > sites;
};

inline void check_distinct_sites(const nsite_group& p)
{
  std::vector<std::size_t> s(p.sites);
  std::sort(s.begin(), s.end());
  if (std::adjacent_find(s.begin(), s.end()) != s.end())
    boost::throw_exception(std::runtime_error(
      "an <NSITE> of type " + boost::lexical_cast<std::string>(p.type)
      + " contains the same site twice (lattice too small for this n-site group?)"));
}

// Reads <NSITE type="t"> <SITE id="1"/> ... </NSITE> inside a <GRAPH>.
// Site ids are 1-based, as for <EDGE source= target=>.
inline nsite_group read_graph_nsite_group_xml(const XMLTag& intag, std::istream& in)
{
  XMLTag tag(intag);
  nsite_group p;
  if (tag.attributes["type"] != "")
    p.type = boost::lexical_cast<type_type>(tag.attributes["type"]);
  if (tag.type != XMLTag::SINGLE) while (true) {
    tag = parse_tag(in);
    if (tag.name == "/NSITE")
      break;
    if (tag.name != "SITE")
      boost::throw_exception(std::runtime_error("illegal element <" + tag.name + "> in <NSITE>"));
    std::string id = tag.attributes["id"];
    if (id == "")
      boost::throw_exception(std::runtime_error("<SITE> in a <GRAPH> <NSITE> needs an id attribute"));
    int i = boost::lexical_cast<int>(id);
    if (i < 1)
      boost::throw_exception(std::runtime_error("illegal site id " + id + " in <NSITE>"));
    p.sites.push_back(i - 1);
    if (tag.type != XMLTag::SINGLE) {
      tag = parse_tag(in);
      if (tag.name != "/SITE")
        boost::throw_exception(std::runtime_error("</SITE> tag missing in <NSITE>"));
    }
  }
  if (p.sites.empty())
    boost::throw_exception(std::runtime_error("empty <NSITE>"));
  check_distinct_sites(p);
  return p;
}

} // end namespace alps

#endif // ALPS_LATTICE_NSITE_GROUP_H
