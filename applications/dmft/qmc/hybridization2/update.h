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

#endif
