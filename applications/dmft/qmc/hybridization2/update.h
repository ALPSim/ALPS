/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2012 by Emanuel Gull <gull@phys.columbia.edu>
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

#ifndef ___UPDATE___
#define ___UPDATE___

#include "impurity.h"


// invert matrix A and calculate its determinant
void invert(blas_matrix & A, double & det);

void construct_matrix(blas_matrix & M, segment_container_t & segments, double BETA,  hybridization_t & F);

double construct_inverse(blas_matrix & M, segment_container_t& segments, double BETA,  hybridization_t &F) ;

// determine F(\tau)
inline double interpolate_F(double t, double BETA, hybridization_t& F) {

  double sign=1;
  if (t<0) {
    t += BETA;
    sign=-1;
  }

  int N = F.size()-1;
  double n = t/BETA*N;
  int n_lower = (int)n; // interpolate linearly between n_lower and n_lower+1
  
  return sign*(F[n_lower] + (n-n_lower)*(F[n_lower+1]-F[n_lower]));

}


// compute distances up/down to the next segment and iterators of these segments
// note: s_down always points to a physical segment, while s_up may point to segments.end() 
void compute_intervals(double t, double BETA, double& t_up, double& t_down, segment_container_t& segments, segment_container_t::iterator& s_up, segment_container_t::iterator& s_down);

// compute overlap between a segment and a list of segments
// requires segment with 0<=t_begin<t_end<=BETA
double segment_overlap(times segment, segment_container_t& other_segments, int other_full_line, double BETA);

double compute_overlap(times segment, segment_container_t& other_segments, int other_full_line, double BETA);

// functions required to compute determinant ratios and perform fast matrix updates 

double det_rat_up(times & new_segment, blas_matrix & M, segment_container_t& segments_old, hybridization_t & F, vector_t& Fs, vector_t& Fe, double BETA, double & det_rat_sign, double & overlap); 

void compute_M_up(int k, blas_matrix & M, vector_t& Fs, vector_t& Fe, double det_rat);

double det_rat_down(std::size_t k, blas_matrix & M, segment_container_t& segments_old, double & det_rat_sign);

void compute_M_down(int k, blas_matrix & M);

// move segment without changin its length
double det_rat_move(times & new_segment, int k, blas_matrix & M, segment_container_t& segments_old, hybridization_t& F, double BETA, double & det_rat_sign, double & overlap);

void compute_M_move(times & new_segment, int k, blas_matrix & M, segment_container_t& segments_old, hybridization_t& F, double BETA, double det_rat);

// shift end point of segment
double det_rat_shift(times & new_segment, std::size_t k, blas_matrix & M, segment_container_t& segments_old, hybridization_t& F, double BETA, double & det_rat_sign, double & overlap);

void compute_M_shift(times & new_segment, std::size_t k, blas_matrix & M, segment_container_t & segments_old, hybridization_t& F, double BETA, double det_rat);

double det_rat_insert_anti(times & anti_segment, blas_matrix & M, segment_container_t& segments_old, hybridization_t& F, double BETA, double & det_rat_sign, double & overlap, vector_t& R);

void compute_M_insert_anti(times & anti_segment, int s, int r, blas_matrix & M, segment_container_t& segments_old, hybridization_t& F, double BETA, double det_rat, vector_t& R);

double det_rat_remove_anti(times anti_segment, int r, std::size_t s, blas_matrix & M, segment_container_t& segments_old, hybridization_t& F, double BETA, double & det_rat_sign);

void compute_M_remove_anti(blas_matrix & M, int s, int r); 

double get_occupation(segment_container_t &segments, int full_line, double tau, double BETA);

#endif
