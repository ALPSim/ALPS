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

#ifndef ALPS_DMFT_SELFCONSISTENCY_H
#define ALPS_DMFT_SELFCONSISTENCY_H

/// @file selfconsistency.h
/// @brief declares the selfconsistency loop functions

#include "solver.h"
#include "hilberttransformer.h"
#include "green_function.h"

/// performs a DMFT self-conistency loop until convergence
///
/// @param parms The input parameters for the simulation.
/// @param solver   The impurity solver. It is left in the state after the final iteration to 
///                 retrieve additional information.
/// @param hilbert  a Hilbert transformation object. It performs the Hilbert transformation
///                 for the density of states of the given model.
/// @param G0_tau The bare Green's function in imaginary time
/// @param G_tau  The Green's function G in imaginary time.

extern void selfconsistency_loop
(
 alps::Parameters& parms, 
 ImpuritySolver& solver, 
 HilbertTransformer& hilbert,
 itime_green_function_t& G0_tau, 
 itime_green_function_t& G_tau
 );

/// performs a DMFT self-conistency loop until convergence
///
/// @param parms The input parameters for the simulation.
/// @param solver   An impurity solver that takes both the bare GF in imaginary time AND in Matsubara freqencies, but returns
/// @param hilbert  a Hilbert transformation object. It performs the Hilbert transformation, taking its arguments in Frequency space.
/// @param G0_omega The bare Green's function in Matsubara frequencies.
/// @param G_omega  The dressed Green's function G in Matsubara frequencies.
void selfconsistency_loop_omega_tau(const alps::Parameters& parms, FrequencyTimeImpuritySolver& solver, FrequencySpaceHilbertTransformer& hilbert,
									matsubara_green_function_t& G0_omega, matsubara_green_function_t& G_omega);


#endif

