/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
*                       Matthias Troyer <troyer@comp-phys.org>
*                       Emanuel Gull <gullc@comp-phys.org>
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

/// @file selfconsistency.C
/// @brief implements the selfconsistency loop functions

#include "selfconsistency.h"
/// @brief Run the self consistency loop for G and G0 mainly in imaginary time. Perform Fourier transformations if needed.
///
/// @param parms contains the ALPS parameters needed for the simulation.
/// @param solver is the impurity solver (e.g. Hirsch Fye) that creates G out of G0.
/// @param hilbert is the HilbertTransformer that solves the Dyson equation, i.e. generates G0 out of G
/// @param G0 is the bare Green's function in imaginary time. It has to be provided as an initial guess
/// @param G is the Green's function, it does not have to be initialized but reasonable values will be returned upon completion of the loop
void selfconsistency_loop(alps::Parameters& parms, ImpuritySolver& solver, HilbertTransformer& hilbert, itime_multiple_vector_t& G0, itime_multiple_vector_t& G)
{
  
  int N = static_cast<int>(parms["N"]);
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  double h = static_cast<double>(parms["H"]);
  double converged = static_cast<double>(parms["CONVERGED"]);
  int paramagnet = static_cast<int>(parms["PARAMAGNET"]);
  
  multiple_vector_type G0_old(G0);
    
  double max_diff;	
  do {
    G0_old = G0;
    G = solver.solve(G0, parms);
    G0 = hilbert(G, mu, h, beta, paramagnet);
	
	max_diff=0;
    for(int i=0; i<N; i++) {
	  if (fabs(G0.first[i]-G0_old.first[i]) > max_diff)
	    max_diff = fabs(G0.first[i]-G0_old.first[i]);
	  if (fabs(G0.second[i]-G0_old.second[i]) > max_diff)
	    max_diff = fabs(G0.second[i]-G0_old.second[i]);	
	}	
	
  } while (max_diff > converged);
  
}
