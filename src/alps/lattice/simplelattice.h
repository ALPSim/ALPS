/***************************************************************************
* ALPS++/lattice library
*
* lattice/simplelattice.h    the lattice class
*
* $Id$
*
* Copyright (C) 2001-2003 by Matthias Troyer <troyer@comp-phys.org>
*                            Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS library, published under the 
* ALPS Library License; you can use, redistribute it and/or modify 
* it under the terms of the License, either version 1 or (at your option) 
* any later version.
*
* You should have received a copy of the ALPS Library License along with 
* the ALPS Library; see the file License.txt. If not, the license is also 
* available from http://alps.comp-phys.org/. 

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

#ifndef ALPS_LATTICE_SIMPLELATTICE_H
#define ALPS_LATTICE_SIMPLELATTICE_H

#include <alps/config.h>
#include <alps/lattice/cell_traits.h>
#include <alps/lattice/dimensional_traits.h>
#include <alps/lattice/coordinate_traits.h>
#include <alps/lattice/simplecell.h>
#include <alps/lattice/lattice.h>

namespace alps {

template <class UnitCell=EmptyUnitCell, class Cell = simple_cell<UnitCell> >
class simple_lattice {
public:
  typedef UnitCell unit_cell_type;
  typedef Cell cell_descriptor;
  typedef typename alps::dimensional_traits<unit_cell_type>::dimension_type dimension_type;
  typedef typename cell_traits< cell_descriptor>::offset_type offset_type;
  
  simple_lattice() {}

  template <class U2, class C2>
  simple_lattice(const simple_lattice<U2,C2>& l)
   : unit_cell_(alps::unit_cell(l))
  {    std::cerr << "simple_lattice\n";}
  
  template <class U2, class C2>
  const simple_lattice& operator=(const simple_lattice<U2,C2>& l)
  {
    unit_cell_ = alps::unit_cell(l);
    return *this;
  }

  const unit_cell_type& unit_cell() const { return unit_cell_;}
  
  cell_descriptor cell(offset_type o) const 
  { return cell_descriptor(unit_cell(),o); }
  
  dimension_type dimension() const { return alps::dimension(unit_cell_);}
protected:
  unit_cell_type unit_cell_;
};

template <class U, class C>
inline typename dimensional_traits<simple_lattice<U,C> >::dimension_type
dimension (const simple_lattice<U,C>& l)
{
  return l.dimension();
}

template <class UnitCell, class Cell>
struct lattice_traits<simple_lattice<UnitCell,Cell> >
{
  typedef typename simple_lattice<UnitCell,Cell>::unit_cell_type unit_cell_type;
  typedef typename simple_lattice<UnitCell,Cell>::cell_descriptor cell_descriptor;
  typedef typename simple_lattice<UnitCell,Cell>::offset_type offset_type;
};

} // end namespace alps

#endif // ALPS_LATTICE_SIMPLELATTICE_H
