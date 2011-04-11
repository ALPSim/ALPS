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
