/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
*                       Matthias Troyer <troyer@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
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

/* $Id$ */

/// @file hirschfye.C
/// @brief implements the Hirsch-Fye solver
/// @sa HirschFyeSolver

#include "hirschfye.h"
#include "hirschfyesim.h" // your old hirschfye.h
#include <alps/osiris.h>
#include <cstdio>


multiple_vector_type HirschFyeSolver::solve(const multiple_vector_type& G0, const alps::Parameters& parms)
{
  alps::ProcessList where(1,alps::local_process());
  /*
  HirschFyeSimulation sim(where,parms,0); // your old IsingSimulation
  sim.set_GO(G0); // to implement
  
  while (sim.work_done() <1.)
    sim.run();
  */
  multiple_vector_type g;// = sim.get_G(); // to implement

  
  // .... maybe write a checkpoint ....
  return g;
}


