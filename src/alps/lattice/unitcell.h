/***************************************************************************
* ALPS++/lattice library
*
* lattice/unitcell.h     the unit cell of a lattice
*
* $Id$
*
* Copyright (C) 2001-2003 by Matthias Troyer <troyer@comp-phys.org>
*                            Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person or organization 
* obtaining a copy of the software covered by this license (the "Software") 
* to use, reproduce, display, distribute, execute, and transmit the Software, 
* and to prepare derivative works of the Software, and to permit others
* to do so for non-commerical academic use, all subject to the following:
*
* The copyright notice in the Software and this entire statement, including 
* the above license grant, this restriction and the following disclaimer, 
* must be included in all copies of the Software, in whole or in part, and 
* all derivative works of the Software, unless such copies or derivative 
* works are solely in the form of machine-executable object code generated by 
* a source language processor.

* In any scientific publication based in part or wholly on the Software, the
* use of the Software has to be acknowledged and the publications quoted
* on the web page http://www.alps.org/license/ have to be referenced.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
**************************************************************************/

#ifndef ALPS_LATTICE_UNITCELL_H
#define ALPS_LATTICE_UNITCELL_H

#include <alps/config.h>
#ifndef ALPS_WITHOUT_XML
# include <alps/parser/parser.h>
#endif
#include <alps/lattice/graph.h>
#include <alps/lattice/graphproperties.h>
#include <alps/lattice/dimensional_traits.h>
#include <boost/graph/adjacency_list.hpp>

namespace alps {

class EmptyUnitCell {
public:
  EmptyUnitCell(std::size_t d=0) : dim_(d) {}
  std::size_t dimension() const {return dim_;}
private:
  std::size_t dim_;	
};

inline dimensional_traits<EmptyUnitCell>::dimension_type
dimension(EmptyUnitCell c)
{
  return c.dimension();
}

class GraphUnitCell
{
public:
  typedef std::vector<int> offset_type;
  typedef detail::coordinate_type coordinate_type;
  typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::directedS,
                                // vertex property
                                boost::property<coordinate_t,detail::coordinate_type,
				  boost::property<vertex_type_t,int> >,
				// edge property
				boost::property<target_offset_t,offset_type,
				  boost::property<source_offset_t,offset_type,
				    boost::property<edge_type_t,int > > >
				> graph_type;

  GraphUnitCell();
  GraphUnitCell(const EmptyUnitCell& e);
#ifndef ALPS_WITHOUT_XML
  GraphUnitCell(const alps::XMLTag&, std::istream&);
#endif

  const GraphUnitCell& operator=(const EmptyUnitCell& e);

#ifndef ALPS_WITHOUT_XML
  void write_xml(std::ostream&, const std::string& = "") const;
#endif

  graph_type& graph() { return graph_;}
  const graph_type& graph() const { return graph_;}
  std::size_t dimension() const { return dim_;}
  const std::string& name() const { return name_;}
  
private:	
  graph_type graph_;
  std::size_t dim_;
  std::string name_;
};

template<>
struct graph_traits<GraphUnitCell> {
  typedef GraphUnitCell::graph_type graph_type;
};

inline dimensional_traits<GraphUnitCell>::dimension_type
dimension(const GraphUnitCell& c)
{
  return c.dimension();
}

typedef std::map<std::string,GraphUnitCell> UnitCellMap;

} // end namespace alps

#ifndef ALPS_WITHOUT_XML

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

inline std::ostream& operator<<(std::ostream& out, const alps::GraphUnitCell& u)
{
  u.write_xml(out);
  return out;	
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace alps
#endif

#endif

#endif // ALPS_LATTICE_UNITCELL_H
