/***********************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Emanuel Gull <gull@phys.columbia.edu>,
 *                              Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
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
 ***********************************************************************************/

#ifndef ___IMP___
#define ___IMP___

#include <alps/ngs.hpp>

#include <alps/alea.h>
#include <cmath>
#include <alps/numeric/detail/matrix.hpp>
#include <alps/numeric/detail/vector.hpp>
#include "green_function.h"
#include "alps_solver.h"

//typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> blas_matrix;
typedef blas::matrix blas_matrix;
typedef boost::variate_generator<boost::mt19937, boost::uniform_real<> > prng_t;
typedef std::vector<double> vector_t;

class times
{
public:
  times() {t_start_=0; t_end_=0; };
  times(double t_start, double t_end) {t_start_=t_start; t_end_=t_end; };
  double t_start() const {return t_start_;} // begin of segment
  double t_end() const {return t_end_;} // end of segment
  void set_t_start(double t_start) {t_start_ = t_start;}
  void set_t_end(double t_end) {t_end_ = t_end;}
  
private:
  double t_start_, t_end_;
};

inline bool operator<(const times& t1, const times& t2) {
  return t1.t_start() < t2.t_start();
}
inline bool operator<(const times& t1, const double t2) {
  return t1.t_start() < t2;
}

inline bool operator>(times t1, times t2) {
  return t1.t_start() > t2.t_start();
}

inline bool operator==(times t1, times t2) {
  return t1.t_start() == t2.t_start();
}

typedef std::set<times> segment_container_t;
typedef std::vector<double> hybridization_t;
typedef std::vector<std::vector<double> > hybridization_container_t;


class hybridization : public alps::mcbase
{
public:
  hybridization(const alps::params& p, int rank);
  void update();
  void measure();
  //getting data into the ALPS solver
  //...from an external program
  void read_external_input_data(const parameters_type &parms); 
  //from the ALPS self consistency cycle
  void read_alps_framework_input_data_omega(const parameters_type &parms);
  void read_alps_framework_input_data_tau(const parameters_type &parms);
  
  double fraction_completed() const {return std::max((sweeps-thermalization_sweeps)/(double)total_sweeps,0.);}
  double work_done() const;
  bool is_thermalized() const;
  
  void create_measurements();
  void set_measurement_vectors();
  void reset_measurement_vectors();
  void resize_measurement_vectors(int crank);
  
private:
  int sweeps;                                        // sweeps done so far
  int thermalization_sweeps;                         // sweeps to be done for equilibration
  int total_sweeps;                                  // sweeps to be done after equilibration
  //double mu;                                       // chemical potential
  std::vector<double> mu_e;                          // mu-<\epsilon>
  blas_matrix u;                                    // interaction matrix
  double t;                                          // energy unit
  hybridization_container_t F;                       // F_up(\tau) = -G_{0,down}^{-1}(-\tau) + (iw + mu) hybridization function
  hybridization_container_t Fomega_imag;
  std::vector<segment_container_t >  segments;       // stores configurations with 0,1,... segments (but not full line)
  std::vector<int>          full_line;               // if 1 means that particle occupies full time-line
  std::vector<double>        sign;                   // sign of Z_n_up
  std::vector<blas_matrix>      M;                   // inverse matrix for up-spins
  std::valarray<double> G_meas;
  std::valarray<double> order_meas;
  std::valarray<double> n_meas;
  std::valarray<double> matrix_size;

  std::valarray<double> gw_meas_re;
  std::valarray<double> gw_meas_im;
  std::valarray<double> fw_meas_re;
  std::valarray<double> fw_meas_im;
  std::valarray<double> gl_meas;
  std::valarray<double> fl_meas;

  std::valarray<std::complex<double> >  mw_meas;
  std::valarray<std::complex<double> > nmw_meas;
  std::valarray<double> g2w_meas_re;
  std::valarray<double> g2w_meas_im;
  std::valarray<double> hw_meas_re;
  std::valarray<double> hw_meas_im;

  std::vector<std::vector<double> > n_vectors;
  std::valarray<double> nn_meas;
  std::valarray<double> nnt_meas;

  double sign_meas;
 
  const double BETA;
  const int FLAVORS;
  const int N_order;
  const int N_meas;
  const int N_shift;
  const int N_swap;

  const int N;
  const int N_nn;
  const int N_w;
  const int N_w2;
  const int N_W;
  const int N_l;
  int N_w_aux;
  
  const int MEASURE_gw;
  const int MEASURE_fw;
  const int MEASURE_gl;
  const int MEASURE_fl;
  const int MEASURE_g2w;
  const int MEASURE_hw;
  const int MEASURE_nn;
  const int MEASURE_nnt;
  
};

std::ostream &operator<<(std::ostream & os, segment_container_t segments);
#endif

