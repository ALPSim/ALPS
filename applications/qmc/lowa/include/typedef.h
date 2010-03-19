/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Applications
*
* Copyright (C) 2006-2010 by Lode Pollet <lpollet@physics.harvard.edu>,
*                            Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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
*****************************************************************************/

/* $Id: typedefs.h 3520 2010-03-03 16:55:00Z tamama $ */


#ifndef ALPS_APPLICATIONS_LOWA_TYPEDEF_H
#define ALPS_APPLICATIONS_LOWA_TYPEDEF_H

#include <iostream>
#include <alps/config.h>

namespace alps {
namespace applications {
namespace lowa {

  typedef int32_t site_type;                  // do not make unsigned!!!  (Lode)
  typedef int32_t index_type;                 // do not make unsigned!!!  (Lode)
  typedef int32_t dist_type;                  // do not make unsigned!!!  (Lode)
  typedef uint8_t dim_type;
  typedef uint8_t fock_basis_type;
  typedef uint32_t particle_number_type;
  typedef uint32_t vertex_type;
  typedef double  time_type;
  typedef double  real_coordinate_type;
  typedef double  parm_type;
  typedef double  obs_type;

  typedef uint8_t kink_assoc_size_type;
  typedef uint8_t kink_assoc_index_type;

} // ending namespace lowa
} // ending namespace applications
} // ending namespace alps

#endif
