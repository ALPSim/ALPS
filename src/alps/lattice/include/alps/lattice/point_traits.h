/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_LATTICE_POINT_TRAITS_H
#define ALPS_LATTICE_POINT_TRAITS_H

namespace alps {

template <class P>
struct point_traits
{
  typedef P vector_type;        
};

} // end namespace alps

#endif // ALPS_LATTICE_POINT_TRAITS_H
