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

#ifndef ___MOVES___
#define ___MOVES___

#include "impurity.h"
#include "update.h"

// length of inserted segment (mu>0) or anti-segment (mu<0)
inline double compute_length(double r, double l_max, double mu) {

  if (mu == 0)
    return r*l_max;  
  else
    return 1/mu*log(r*(exp(mu*l_max)-1)+1);
}


void insert_remove_full_line(prng_t& rng, double mu, double BETA, blas_matrix &u, int& full_line, std::vector<segment_container_t>& other_segments, std::vector<int>& other_full_line, int this_flavor);

void insert_remove_segment(prng_t& rng, double t, double BETA, double mu, blas_matrix &u, hybridization_t& F, segment_container_t& segments, blas_matrix & M, double & sign, std::vector<segment_container_t>& other_segments, std::vector<int> other_full_line, int this_flavor);

void insert_remove_antisegment(prng_t& rng, double t, double BETA, double mu,blas_matrix &u, hybridization_t& F, int& full_line, segment_container_t& segments, blas_matrix & M, double & sign, std::vector<segment_container_t>& other_segments, std::vector<int> other_full_line, int this_flavor);


// shift segment
void shift_segment(prng_t& rng, segment_container_t& segments, double BETA, double mu, blas_matrix &u, hybridization_t& F, blas_matrix & M, double & sign, std::vector<segment_container_t>& other_segments, std::vector<int>& other_full_line, int this_flavor) ;
// swap segment configurations
void swap_segments(prng_t& rng, double BETA, hybridization_t& F_up, hybridization_t& F_down, segment_container_t& segments_up, segment_container_t& segments_down, int& full_line_up, int& full_line_down, double & sign_up, double & sign_down, blas_matrix & M_up, blas_matrix& M_down);

void measure_GF_imp(prng_t& rng, int N, double BETA, double mu, blas_matrix &u, hybridization_t& F, segment_container_t& segments, blas_matrix & M, std::vector<segment_container_t>& other_segments, std::vector<int> other_full_line, int this_flavor, std::valarray<double> & G_meas);

#endif
