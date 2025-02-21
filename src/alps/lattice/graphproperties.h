/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
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

#ifndef ALPS_LATTICE_GRAPH_PROPERTIES_H
#define ALPS_LATTICE_GRAPH_PROPERTIES_H

#include <vector>
#include <boost/graph/properties.hpp>

namespace alps {

struct vertex_type_t { typedef boost::vertex_property_tag kind; };
typedef vertex_type_t site_type_t;
struct coordinate_t { typedef boost::vertex_property_tag kind; };
struct parity_t { typedef boost::vertex_property_tag kind; };

struct edge_type_t { typedef boost::edge_property_tag kind; };
typedef edge_type_t bond_type_t;
struct source_offset_t { typedef boost::edge_property_tag kind; };
struct target_offset_t { typedef boost::edge_property_tag kind; };
struct boundary_crossing_t { typedef boost::edge_property_tag kind; };
struct edge_vector_t { typedef boost::edge_property_tag kind; };
typedef edge_vector_t bond_vector_t;
struct edge_vector_relative_t { typedef boost::edge_property_tag kind; };
typedef edge_vector_relative_t bond_vector_relative_t;

struct graph_name_t { typedef boost::graph_property_tag kind; };
struct dimension_t { typedef boost::graph_property_tag kind; };

using boost::vertex_index_t;
typedef vertex_index_t site_index_t;
using boost::edge_index_t;
typedef edge_index_t bond_index_t;

typedef std::vector<double> coordinate_type;
typedef std::vector<int> offset_type;
typedef std::vector<int> distance_type;
typedef unsigned int type_type;

} // end namespace alps

#endif // ALPS_LATTICE_GRAPH_PROPERTIES_H
