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

#ifndef ALPS_LATTICE_BOUNDARY_H
#define ALPS_LATTICE_BOUNDARY_H

#include <alps/config.h>
#include <alps/lattice/propertymap.h>
#include <alps/osiris/dump.h>
#include <boost/limits.hpp>
#include <boost/pending/property.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <string>
#include <vector>

namespace alps {

struct boundary_crossing {
  typedef unsigned int dimension_type;
  typedef int direction_type;

  boundary_crossing() : bc(0) {}
  operator bool() const { return bc!=0;}
  
  direction_type crosses(dimension_type d) const 
  { 
    return (bc&(1<<2*d)) ? +1 : ( (bc&(2<<2*d)) ? -1 : 0);
  }
  
  const boundary_crossing& set_crossing(dimension_type d, direction_type dir) 
  { 
    bc &= ~(3<<2*d);
    bc |= (dir>0 ? (1<<2*d) : (dir <0? (2<<2*d) : 0));
    return *this;
  }
  
  const  boundary_crossing& invert() 
  {
    integer_type rest=bc;
    int dim=0;
    while (rest) {
      invert(dim);
      dim++;
      rest >>=2;
    }
    return *this;
  }
  
  void save (ODump& dump) const { dump << bc;}
  void load (IDump& dump) { dump >> bc;}
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  { ar & bc; }
private:  
  typedef uint8_t integer_type;
  integer_type bc;
  const  boundary_crossing& invert(dimension_type d) {
    integer_type mask = 3<<2*d;
    if (bc&mask)
      bc^=mask;
    return *this;
  }
};

} // end namespace alps

#endif // ALPS_LATTICE_GRAPH_PROPERTIES_H
