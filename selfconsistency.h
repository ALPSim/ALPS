/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2014 by Emanuel Gull <gullc@comp-phys.org>
 *
 * To be used with permission by the authors
 *
 *****************************************************************************/

/* $Id$ */

#ifndef ALPS_DMFT_SELFCONSISTENCY_H
#define ALPS_DMFT_SELFCONSISTENCY_H

/// @file selfconsistency.h
/// @brief declares the selfconsistency loop functions

#include "solver.h"
#include "hilberttransformer.h"
#include "externalsolver.h"
#include "cluster.h"
#include "fouriertransform.h"
#include "green_function.h"

/// performs a DMFT self-conistency loop until convergence
void selfconsistency_loop_single_site(const alps::Parameters& parms, MatsubaraImpuritySolver& solver,  FrequencySpaceHilbertTransformer& hilbert);


void selfconsistency_loop_DCA(const alps::Parameters& parms, MatsubaraImpuritySolver& solver, DCATransformer& clustertrans);
void selfconsistency_loop_DCA(const alps::Parameters& parms, NambuImpuritySolver& solver, DCATransformer& clustertrans);
void selfconsistency_loop_CDMFT(const alps::Parameters& parms, MatsubaraImpuritySolver& solver, CDMFTTransformer& clustertrans);
void selfconsistency_loop_CDMFT(const alps::Parameters& parms, NambuImpuritySolver& solver, CDMFTTransformer& clustertrans);


void perform_single_site_loop(const alps::Parameters &parms);

void perform_DCA_loop(const alps::Parameters &parms);
void perform_CDMFT_loop(const alps::Parameters &parms);
void perform_cluster_loop(const alps::Parameters &parms);

FrequencySpaceHilbertTransformer *setup_DOS_hilbert_transform(const alps::Parameters &parms);
FrequencySpaceHilbertTransformer *setup_2D_hilbert_transform(const alps::Parameters &parms);
FrequencySpaceHilbertTransformer *setup_Bethe_hilbert_transform(const alps::Parameters &parms);

#endif

