/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
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

/* $Id$ */

/// @file hirschfyesim.h
/// @brief the actual Hirsch-Fye simulation
///
/// @todo needs to be implemented as soon as possible by Philipp Werner


#include "hirschfyesim.h"
#include <cmath>
using namespace std;
using namespace alps;

// compute largest element of a matrix
double norm_max(dense_matrix & m) {

  double max(0);
  for (dense_matrix::const_iterator1 it = m.begin1(); it != m.end1(); ++it)
    if (std::fabs((*it)) > max)
	  max = std::fabs((*it));
	  
  return max;

}

// compute Green's function for given spin configuration from Green0
void update_from_zero(dense_matrix & Green, dense_matrix & Green0, std::vector<int> & spins, double l) {

  int N(spins.size());

  // calculate a = 1 + (1-g)*(exp(v')-1) 
  dense_matrix a(N, N);
  
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
	  a(i,j) = -Green0(i,j)*(exp(l*spins[j])-1);
	}
	a(i,i) += exp(l*spins[i]); 
  }
  
  // solve a*Green=Green0
  // gesv solves A*X=B, input for B is Green0, output (=solution X) is Green   
  Green = Green0; 
  boost::numeric::bindings::lapack::gesv(a,Green);
  
}


HirschFyeSimulation::HirschFyeSimulation(const alps::ProcessList& where,const alps::Parameters& p,int node, multiple_vector_type& G0)
  : alps::scheduler::MCRun(where,p,node),  
    sweeps(0),
    thermalization_sweeps(static_cast<int>(parms["THERMALIZATION"])),
    total_sweeps(static_cast<int>(parms["SWEEPS"])),
    beta(static_cast<double>(parms["BETA"])),
    N(static_cast<int>(parms["N"])),	
    u(static_cast<double>(parms["U"])),
	N_check(1000), // value is automatically adjusted during the simulation
	check_counter(0),
	tolerance(static_cast<double>(parms["TOLERANCE"])),	
	G0_up(static_cast<int>(parms["N"])),	
	G_up(static_cast<int>(parms["N"]),0.),
	Green0_up(static_cast<int>(parms["N"]), static_cast<int>(parms["N"])),
	Green_up(static_cast<int>(parms["N"]), static_cast<int>(parms["N"])),
	G0_down(static_cast<int>(parms["N"])),	
	G_down(static_cast<int>(parms["N"]),0.),
	Green0_down(static_cast<int>(parms["N"]), static_cast<int>(parms["N"])),
	Green_down(static_cast<int>(parms["N"]), static_cast<int>(parms["N"])),			
    spins(static_cast<int>(parms["N"]), 0) // at the beginning, spins are formally set to 0

{

  // initialize G0
  for (int i=0; i<G0.fist.size(); i++) {
    G0_up[i] = G0.first[i];
	G0_down[i] = G0.second[i];
  }

  // initialize lambda
  double tmp = exp(beta/N*u/2);
  lambda = log(tmp + sqrt(pow(tmp,2)-1));

  // initialize bath Green's function matrices
  for (int i=0; i<N; i++) {
    for (int j=i; j<N; j++) {	  
      Green0_up(j,i) = G0_up[j-i];
      Green0_down(j,i) = G0_down[j-i];
    }
  }
  for (int i=0; i<N; i++) {
    for (int j=i+1; j<N; j++) {
      Green0_up(i,j) = -G0_up[N-(j-i)];
      Green0_down(i,j) = -G0_down[N-(j-i)];
	}
  }

  // choose random spin configuration
  for (int i=0; i<spins.size(); i++)
    spins[i] = (random_01()<0.5 ? 1 : -1);

  // initialize Green's function matrices
  update_from_zero(Green_up, Green0_up, spins, lambda);
  update_from_zero(Green_down, Green0_down, spins, -lambda);
 
  // create measurement objects
  measurements << RealVectorObservable("G_meas_up");
  measurements << RealVectorObservable("G_meas_down");  
  
}

void HirschFyeSimulation::load(alps::IDump& dump) {
  dump >> sweeps;
  if(!where.empty()) {// skip reading the spins if we are just evaluating
    dump >> G0_up >> G_up >> G0_down >> G_down >> spins; 	 
  }
  else
    measurements.compact(); // and throw away time series
}

void HirschFyeSimulation::save(alps::ODump& dump) const
{
  dump << G0_up << G_up << G0_down << G_down << spins;
}

bool HirschFyeSimulation::change_parameter(const std::string& name, const alps::StringValue& value)
{
  if(name=="SWEEPS")
    total_sweeps=static_cast<uint32_t>(value);
  else if (name=="THERMALIZATION" && !is_thermalized())
    thermalization_sweeps=static_cast<uint32_t>(value);
  else
    return false; // cannot do it
  return true; // could do it
}

bool HirschFyeSimulation::is_thermalized() const
{
  return (sweeps >= thermalization_sweeps);
}

double HirschFyeSimulation::work_done() const
{
  return (is_thermalized() ? (sweeps-thermalization_sweeps)/double(total_sweeps) : 0.);
}

void HirschFyeSimulation::dostep()
{

  // increment sweep count
  sweeps++;
  check_counter++;
 
  // do local update
  for (int i=0; i<N; i++)
    update_single_spin(random_01(), Green_up, Green_down, spins, lambda);
	  
  // calculate Greens function	
  std::valarray<double> G_tmp_up(N), G_tmp_down(N);			
  G_tmp_up = 0.;
  G_tmp_down = 0.;
  
  for (int i=0; i<N; i++) {
    for (int j=i; j<N; j++) {
      G_tmp_up[j-i] += Green_up(j,i);
	  G_tmp_down[j-i] += Green_down(j,i); 
	  
	  if (j>i) {
        G_tmp_up[N-(j-i)] += -Green_up(i,j);
	    G_tmp_down[N-(j-i)] += -Green_down(i,j); 	  
	  }
	  
    }
  }
  
  G_tmp_up /= N;
  G_tmp_down /= N;
  
  // measure Greens function
  measurements.get<RealVectorObservable>("G_meas_up") << G_tmp_up;
  measurements.get<RealVectorObservable>("G_meas_down") << G_tmp_down;

  // check if precision has deteriorated
  if (check_counter%N_check == 0) {
  
    dense_matrix Green_diff_up(Green_up), Green_diff_down(Green_down);

    update_from_zero(Green_up, Green0_up, spins, lambda);
    update_from_zero(Green_down, Green0_down, spins, -lambda);

    Green_diff_up -= Green_up;
	Green_diff_down -= Green_down;

    if (norm_max(Green_diff_up) > tolerance || norm_max(Green_diff_down) > tolerance)
	  N_check /= 2;
	else 
	  N_check *= 2;
	  
	check_counter = 0;

  }

}


multiple_vector_type HirschFyeSimulation::get_G() {

  multiple_vector_type g;

  alps::RealVectorObsevaluator G_meas_up=get_measurements()["G_meas_up"];
  alps::RealVectorObsevaluator G_meas_down=get_measurements()["G_meas_down"];  
  
  // not sure how to implement this correctly
  for (int i=0; i<G_meas_up.size(); i++) {
    g.first[i] = G_meas_up.mean()[i];
	g.second[i] = G_meas_down.mean()[i];
  }
  
  return g;

}

