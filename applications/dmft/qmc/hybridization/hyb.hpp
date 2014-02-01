/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                   
 *  based on an earlier version by Philipp Werner and Emanuel Gull
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
#ifndef HYB_HPP
#define HYB_HPP

#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>
#include <alps/ngs/signal.hpp>
#include "green_function.h"
#include "hybsegment.hpp"
#include "hyblocal.hpp"
#include "hybconfig.hpp"
#include "boost/chrono/chrono.hpp"

#ifdef HYB_SIM_MAIN
boost::uint64_t sweep_count;
std::vector<boost::uint64_t> nacc,nprop;
std::vector<std::string> update_type;
#else
extern boost::uint64_t sweep_count;
extern std::vector<boost::uint64_t> nacc,nprop;
extern std::vector<std::string> update_type;
#endif

class hybridization:public alps::mcbase
{
public:
  //constructor
  hybridization(const alps::params &parms, int crank);
  void sanity_check(const alps::params &parms); //check whether parameters make sense
  void show_info(const alps::params &parms, int crank);
  //Monte Carlo update and measurements functions
  void measure();
  void update();
  bool is_thermalized() const { return (sweeps >= thermalization_sweeps); }
  double fraction_completed() const;
  friend std::ostream &operator<<(std::ostream &os, const hybridization &hyb);

private:
  bool VERBOSE;
  int crank;
  int csize;
  int output_period;
  clock_t start_time;
  clock_t end_time;
  //boost::chrono::steady_clock::time_point lasttime;
  //boost::chrono::steady_clock::duration delay;

  //initialize all measurements and measurement vectors (with 0)
  void create_measurements();
  //measure_* functions perform the actual measurements
  void measure_order();
  void measure_G(std::vector<std::map<double,double> > &F_prefactor);
  void measure_Gw(std::vector<std::map<double,double> > &F_prefactor);
  void measure_Gl(std::vector<std::map<double,double> > &F_prefactor);
  void measure_G2w(std::vector<std::map<double,double> > &F_prefactor);
  void measure_nn();
  void measure_nnt();
  void measure_nnw();
  void measure_sector_statistics();

  //accumulate_* functions pipe the result into ALPS observables and clear measurment vectors
  void accumulate_order();
  void accumulate_G();
  void accumulate_Gw();
  void accumulate_Gl();
  void accumulate_nn();
  void accumulate_nnt();
  void accumulate_nnw();
  void accumulate_sector_statistics();

  //Monte Carlo update routines
  void change_zero_order_state_update();
  void shift_segment_update();
  void insert_remove_segment_update();
  void insert_remove_antisegment_update();
  void insert_remove_spin_flip_update();

  //details of update routines
  void insert_segment_update(int orbital);
  void remove_segment_update(int orbital);
  void insert_antisegment_update(int orbital);
  void remove_antisegment_update(int orbital);
  void spin_flip_update(int orbital);
  
  //programming and debug functions
  double full_weight() const;

  //algorithm parameters
  boost::uint64_t sweeps;
  boost::uint64_t thermalization_sweeps;
  boost::uint64_t total_sweeps;
  std::size_t n_orbitals;
  double sign;

  //physics parameters
  double beta,U_,MU_;
  

  //updates parameters
  std::size_t N_meas;
  double fraction; // Idea by Philipp: Do not use the long time intervals but allow
  // only for a certain fraction when inserting segments; supposedly increases
  // acceptance rate if B*U or B*MU is very large
  
  //measurement parameters
  std::size_t N_w;    //number of Matsubara frequency points
  std::size_t N_l;    //number of Legendre coefficients
  std::size_t N_t;    //number of imag time slices
  std::size_t N_w2;   //number of Matsubara frequency points for two-particle measurements
  std::size_t N_W;    //number of bosonic Matsubara frequency points for two-particle measurements
  std::size_t N_w_aux;//number of Matsubara frequency points for the measurment of M(w1,w2)
  std::size_t N_hist_orders;
  std::size_t N_nn;
  bool spin_flip;
  bool MEASURE_nnt;
  bool MEASURE_nnw;
  bool MEASURE_nn;
  bool MEASURE_g2w;
  bool MEASURE_h2w;
  bool MEASURE_time;
  bool MEASURE_freq;
  bool MEASURE_legendre;
  bool MEASURE_sector_statistics;

  //observable names
  std::vector<std::string> g_names;
  std::vector<std::string> f_names;
  std::vector<std::string> density_names;
  std::vector<std::string> order_names;
  std::vector<std::string> order_histogram_names;
  std::vector<std::string> gwr_names, gwi_names, fwr_names, fwi_names;
  std::vector<std::string> gl_names, fl_names;
  std::vector<std::vector<std::string> > g2wr_names, g2wi_names, h2wr_names, h2wi_names;
  std::vector<std::vector<std::string> > nnt_names, nnw_re_names, nn_names;

  //measurement vectors (initialized first time in create_measurements() )
  double sgn;
  std::vector<std::vector<double> >order_histogram;
  std::vector<double>orders;
  std::vector<double>order_histogram_total;
  std::vector<std::vector<double> >G;
  std::vector<std::vector<double> >F;
  std::vector<double>densities;
  std::vector<std::vector<double> >Gwr;
  std::vector<std::vector<double> >Gwi;
  std::vector<std::vector<double> >Fwr;
  std::vector<std::vector<double> >Fwi;
  std::vector<std::vector<double> >Gl;
  std::vector<std::vector<double> >Fl;
  std::vector<std::vector<std::complex<double> > >G2w;
  std::vector<std::vector<std::complex<double> > >F2w;
  std::vector<std::vector<double> >n_vectors;
  std::vector<std::vector<std::vector<double> > >nnt;
  std::vector<std::vector<std::vector<double> > >nnw_re;
  std::vector<std::vector<double> >nn;
  std::vector<double>sector_statistics;
  std::vector<double>g2wr;
  std::vector<double>g2wi;
  std::vector<double>h2wr;
  std::vector<double>h2wi;

  std::vector<std::map<double,double> > F_prefactor;

  //local impurity operator configuration
  local_configuration local_config;

  //impurity operator configuration
  hybridization_configuration hyb_config;
};

std::ostream &operator<<(std::ostream &os, const hybridization &hyb);

#ifndef COLORS
#define COLORS
#define cblack "\033[22;30m"
#define cred "\033[22;31m"
#define cgreen "\033[22;32m"
#define cbrown "\033[22;33m"
#define cblue "\033[22;34m"
#define cmagenta "\033[22;35m"
#define ccyan "\033[22;36m"
#define cgray "\033[22;37m"
#define cdgray "\033[01;30m"
#define clred "\033[01;31m"
#define clgreen "\033[01;32m"
#define clyellow "\033[01;33m"
#define clblue "\033[01;34m"
#define clmagenta "\033[01;35m"
#define clcyan "\033[01;36m"
#define cwhite "\033[01;37m"
#endif

#endif
