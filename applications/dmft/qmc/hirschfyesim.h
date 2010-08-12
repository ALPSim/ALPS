/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
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

/// @file hirschfyesim.h
/// @brief the actual Hirsch-Fye simulation

#ifndef ALPS_DMFT_HIRSCHFYESIM_H
#define ALPS_DMFT_HIRSCHFYESIM_H

#include "types.h"
#include "hirschfyeaux.h"
#include "alps_solver.h"
#include "green_function.h"

#include <alps/scheduler/montecarlo.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>
#include <alps/alea/detailedbinning.h>
#include <alps/alea/detailedbinning.h>
#include <alps/alea/nobinning.h>

#include <boost/numeric/ublas/io.hpp>

#include <stack>
#include <queue>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>


typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;
typedef alps::SimpleRealVectorObservable vec_obs_t;

class HirschFyeSim : public alps::scheduler::MCSimulation, public alps::MatsubaraImpurityTask
{
public:
  HirschFyeSim(const alps::ProcessList& w, const boost::filesystem::path& p) 
  : alps::scheduler::MCSimulation(w,p)
  { 
  }
  
  HirschFyeSim(const alps::ProcessList& w, const alps::Parameters& p) 
  : alps::scheduler::MCSimulation(w,p) 
  {
  }
  
  std::pair<matsubara_green_function_t, itime_green_function_t>get_result();
  
};


///This class implements the Hirsch Fye algorithm.
class HirschFyeRun : public alps::scheduler::MCRun
{
public:
  ///
  HirschFyeRun(const alps::ProcessList&,const alps::Parameters&,int);
  ///store a dump of this simulation
  void save(alps::ODump&) const;
  ///load a previously saved simulation
  void load(alps::IDump&);
  ///perform one Hirsch Fye step
  void dostep();
  ///return true if the simulation has been thermalized
  bool is_thermalized() const;
  ///return a value between zero and one, one meaning that the simulation is about to end.
  double work_done() const;
  
  ///return the Green's function G
  itime_green_function_t get_result() const;
  
  ///compute vector Green's function from a matrix.
  void green_vector_from_matrix(itime_green_function_t &green_function, const dense_matrix & green_matrix_up, const dense_matrix & green_matrix_down)const;
  ///compute matrix Green's function from the vector.
  void green_matrix_from_vector(const itime_green_function_t & bare_green_function, dense_matrix & green_matrix_up, dense_matrix & green_matrix_down)const;
  
  
private:
  int sweeps;                    // sweeps done
  int thermalization_sweeps;        // sweeps to be done for equilibration
  int total_sweeps;                // sweeps to be done after equilibration
  double beta;                    // inverse temperature
  int N;                            // number of time slices
  int n_site;         // for cluster simulations: number of sites
  double u;                        // on-site interaction
  double lambda;                    // cosh(lambda) = exp(delta_tau*u/2)
  int N_check;                    // number of sweeps with single-spin updates
  int check_counter;                // counts the number of single-spin updates
  int fp_interval;                // counts the interval between measurements of the fourpoint function
  bool measure_fourpoint;
  double tolerance;                // tolerance for deterioration of precision
  dense_matrix Green0_up;           // bath Green's function matrix for up spins
  dense_matrix Green_up;            // Green's function matrix for up spins  
  dense_matrix Green0_down;        // bath Green's function matrix for down spins
  dense_matrix Green_down;          // Green's function matrix for down spins     
  std::vector<int> spins;            // auxiliary Ising spins
  double max_time;
  double start_time;
  int sign;
  
  itime_green_function_t bare_green_tau;
  itime_green_function_t green_tau;
  
  
};

#endif
