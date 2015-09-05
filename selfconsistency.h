/*****************************************************************************
 *
 * Copyright (C) 2005-2015 by Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#pragma once
namespace alps{ class params;}
class ImpuritySolver;
class FrequencySpaceHilbertTransformer;
class DCATransformer;
class CDMFTTransformer;
/// performs a DMFT self-conistency loop until convergence
void selfconsistency_loop_single_site(const alps::params& parms, ImpuritySolver& solver,  FrequencySpaceHilbertTransformer& hilbert);


void selfconsistency_loop_DCA(const alps::params& parms, ImpuritySolver& solver, DCATransformer& clustertrans);
void selfconsistency_loop_DCA_Nambu(const alps::params& parms, ImpuritySolver& solver, DCATransformer& clustertrans);
void selfconsistency_loop_CDMFT(const alps::params& parms, ImpuritySolver& solver, CDMFTTransformer& clustertrans);
void selfconsistency_loop_CDMFT_Nambu(const alps::params& parms, ImpuritySolver& solver, CDMFTTransformer& clustertrans);


void perform_single_site_loop(const alps::params &parms);

void perform_DCA_loop(const alps::params &parms);
void perform_CDMFT_loop(const alps::params &parms);
void perform_cluster_loop(const alps::params &parms);

FrequencySpaceHilbertTransformer *setup_DOS_hilbert_transform(const alps::params &parms);
FrequencySpaceHilbertTransformer *setup_2D_hilbert_transform(const alps::params &parms);
FrequencySpaceHilbertTransformer *setup_Bethe_hilbert_transform(const alps::params &parms);
