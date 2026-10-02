/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2002 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_LATTICE_CELL_TRAITS_H
#define ALPS_LATTICE_CELL_TRAITS_H

namespace alps {

template <class C>
struct cell_traits 
{
  typedef typename C::offset_type offset_type;
};

} // end namespace alps

#endif // ALPS_LATTICE_CELL_TRAITS_H
