/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2004 by Ian McCulloch <ianmcc@physik.rwth-aachen.de>
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

#ifndef ALPS_LATTICE_BOND_COMPARE_H
#define ALPS_LATTICE_BOND_COMPARE_H

#include <alps/lattice/graph_traits.h>

#include <algorithm>

namespace alps {

// bond_descriptor_compare - a comparison functor for bond_descriptor's,
// use for example in std::set<bond_descriptor, bond_descriptor_compare<graph_type> >

template <typename G>
struct bond_descriptor_compare
{
  typedef G graph_type;
  typedef graph_traits<graph_type> traits_type;
  typedef typename traits_type::bond_descriptor bond_descriptor;
  typedef typename traits_type::site_descriptor site_descriptor;

  typedef bool            result_type;
  typedef bond_descriptor first_argument_type;
  typedef bond_descriptor second_argument_type;

  bond_descriptor_compare(G const* Gr) : Graph(Gr) {}

  bool operator()(bond_descriptor b1, bond_descriptor b2) const;

  graph_type const* Graph;
};

template <typename G>
inline
bool bond_descriptor_compare<G>::operator()(bond_descriptor b1, bond_descriptor b2) const
{
  using boost::source;
  using boost::target;
  site_descriptor s1 = source(b1, *Graph); 
  site_descriptor s2 = source(b2, *Graph);
  return s1 < s2 || (s1 == s2 && target(b1, *Graph) < target(b2, *Graph)); 
}

// bond_descriptor_compare_undirected - a comparison functor for bond_descriptor's
// that treats bonds (source, target) and (target, source) as equivalent

template <typename G>
struct bond_descriptor_compare_undirected
{
  typedef G graph_type;
  typedef graph_traits<graph_type> traits_type;
  typedef typename traits_type::bond_descriptor bond_descriptor;
  typedef typename traits_type::site_descriptor site_descriptor;

  typedef bool            result_type;
  typedef bond_descriptor first_argument_type;
  typedef bond_descriptor second_argument_type;

  bond_descriptor_compare_undirected(G const* Gr) : Graph(Gr) {}

  bool operator()(bond_descriptor b1, bond_descriptor b2) const;

  graph_type const* Graph;
};

template <typename G>
inline
bool bond_descriptor_compare_undirected<G>::operator()(bond_descriptor b1, bond_descriptor b2) const
{
  using std::swap;
  using boost::source;
  using boost::target;
  site_descriptor s1 = source(b1, *Graph); 
  site_descriptor s2 = source(b2, *Graph);
  site_descriptor t1 = target(b1, *Graph); 
  site_descriptor t2 = target(b2, *Graph);
  if (s1 < t1) swap(s1, t1);
  if (s2 < t2) swap(s2, t2);
  return s1 < s2 || (s1 == s2 && t1 < t2);
}

} // end namespace alps

#endif // ALPS_LATTICE_BOND_COMPARE_H
