/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2002 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

#ifndef ALPS_LATTICE_SIMPLECELL_H
#define ALPS_LATTICE_SIMPLECELL_H

#include <alps/config.h>
#include <alps/lattice/unitcell.h>
#include <alps/lattice/graph_traits.h>
#include <alps/lattice/dimensional_traits.h>
#include <alps/lattice/cell_traits.h>
#include <alps/utility/vectorio.hpp>

namespace alps {

template <class UnitCell=EmptyUnitCell, class Offset=typename std::vector<int> >
class simple_cell  {
public:
  typedef Offset offset_type;
  typedef UnitCell unit_cell_type;
  typedef typename alps::dimensional_traits<UnitCell>::dimension_type dimension_type;

  simple_cell() : dim_(0) {}
  simple_cell(const unit_cell_type& u, const offset_type& o)
   : dim_(alps::dimension(u)), offset_(o) {}
  
  const offset_type& offset() const { return offset_;}
  dimension_type dimension() { return dim_;}
private:
  dimension_type dim_;
  offset_type offset_;
};

template <class UnitCell,class Offset>
inline typename simple_cell<UnitCell,Offset>::dimension_type
dimension(const simple_cell<UnitCell,Offset>& c)
{
  return c.dimension();
}

} // end namespace alps

#endif // ALPS_LATTICE_SIMPLECELL_H
