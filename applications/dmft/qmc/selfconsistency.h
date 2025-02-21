 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
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

/* $Id: selfconsistency.h 370 2009-08-05 10:08:34Z fuchs $ */

#ifndef ALPS_DMFT_SELFCONSISTENCY_H
#define ALPS_DMFT_SELFCONSISTENCY_H

/// @file selfconsistency.h
/// @brief declares the selfconsistency loop functions

#include "solver.h"
#include "xml.h"
#include "hilberttransformer.h"
#include "fouriertransform.h"
#include "green_function.h"

/// performs a DMFT self-consistency loop until convergence
///
/// @param parms The input parameters for the simulation.
/// @param solver   The impurity solver. It is left in the state after the final iteration to 
///                 retrieve additional information.
/// @param hilbert  a Hilbert transformation object. It performs the Hilbert transformation
///                 for the density of states of the given model.

extern void selfconsistency_loop(alps::Parameters& parms, ImpuritySolver& solver, HilbertTransformer& hilbert);

extern void F_selfconsistency_loop(alps::Parameters& parms, ImpuritySolver& solver, HilbertTransformer& hilbert);


/// performs a DMFT self-consistency loop until convergence
///
/// @param parms The input parameters for the simulation.
/// @param solver   An impurity solver that takes both the bare GF in imaginary time AND in Matsubara freqencies, but returns
/// @param hilbert  a Hilbert transformation object. It performs the Hilbert transformation, taking its arguments in Frequency space.
void selfconsistency_loop_omega(alps::Parameters& parms, MatsubaraImpuritySolver& solver,  FrequencySpaceHilbertTransformer& hilbert);


//void selfconsistency_loop_DCA(const alps::Parameters& parms, MatsubaraImpuritySolver& solver, DCATransformer& clustertrans);

#endif

