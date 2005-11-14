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

#ifndef ALPS_DMFT_HIRSCHFYESIM_H

#include <stack>
#include <queue>
#include <vector>
#include <iostream>
#include <algorithm>
#include <alps/scheduler/montecarlo.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>
#include <cmath>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include "types.h"


inline alps::ODump& operator << (alps::ODump& odump, const dense_matrix& m) {
  odump << m.size1() << m.size2();
  for (int i=0; i<m.size1(); i++) {
    for (int j=0; j<m.size2(); j++) {
          odump << m(i,j);
	}
  }
  return odump;
}


inline alps::IDump& operator >> (alps::IDump& idump, dense_matrix& m) {
  int size1, size2;
  idump >> size1;
  idump >> size2;
  m.resize(size1, size2);
  for (int i=0; i<m.size1(); i++) {
    for (int j=0; j<m.size2(); j++) {
          idump >> m(i,j);
    }
  }  
  return idump;
}


// compute Green's function matrix for a given spin configuration, from Greens function for spins set to zero
void update_from_zero(dense_matrix & Green, dense_matrix & Green0, std::vector<int> & spins, double l);


// attempt single spin flip and (if successful) compute Green's function matrix for the new configuration
template<class RNG> 
void update_single_spin(RNG & rng, dense_matrix & Green_up, dense_matrix & Green_down, std::vector<int> & spins, double lambda) 
{

  int N(spins.size());
  
  // choose random site
  int site = (int)(N*rng());

  // calculate ratio of determinants
  double r_up = 1 + (1-Green_up(site,site))*(exp(-2*lambda*spins[site])-1);
  double r_down = 1 + (1-Green_down(site,site))*(exp(2*lambda*spins[site])-1);
  double det_rat = r_up*r_down; 

  // what about sign changes?
  det_rat = fabs(det_rat);
  
  // if update successful ...
  if (rng() < det_rat/(1+det_rat)) {
  
    // update Green's function
    dense_matrix Green_new_up(Green_up), Green_new_down(Green_down);
	
	double tmp1 = exp(-2*lambda*spins[site])-1;
	double tmp2 = exp(+2*lambda*spins[site])-1;
	for (int i=0; i<N; i++) {
	  for (int j=0; j<N; j++) {	
        Green_new_up(i,j)   += (Green_up(i,site)-(i==site))*tmp1/(1+(1-Green_up(site,site))*tmp1)*Green_up(site,j);
        Green_new_down(i,j) += (Green_down(i,site)-(i==site))*tmp2/(1+(1-Green_down(site,site))*tmp2)*Green_down(site,j);			
      }	
	}
  
    Green_up.swap(Green_new_up); //avoid copying again. never mind the contents of Green_new_up now.
	Green_down.swap(Green_new_down);
	
	// update spin
	spins[site] = -spins[site];
  
  }
}


class HirschFyeSimulation : public alps::scheduler::MCRun
{
public:
  HirschFyeSimulation(const alps::ProcessList&,const alps::Parameters&,int, const multiple_vector_type&);
  void save(alps::ODump&) const;
  void load(alps::IDump&);
  void dostep();
  bool is_thermalized() const;
  double work_done() const;
  bool change_parameter(const std::string& name, const alps::StringValue& value);
  multiple_vector_type get_G();

  ///compute matrix Green's function from the vector.
  void green_vector_from_matrix(std::valarray<double> &green_vector, const dense_matrix & green_matrix)const;
  void green_matrix_from_vector(const std::valarray<double> &green_vector, dense_matrix & green_matrix)const;
  
  
private:
  int sweeps;						// sweeps done
  int thermalization_sweeps;		// sweeps to be done for equilibration
  int total_sweeps;					// sweeps to be done after equilibration
  double beta;						// inverse temperature
  int N;							// number of time slices
  double u;							// on-site interaction
  double lambda;					// cosh(lambda) = exp(delta_tau*u/2)
  int N_check;						// number of sweeps with single-spin updates
  int check_counter;				// counts the number of single-spin updates
  double tolerance;					// tolerance for deterioration of precision
                                    // (maximum deviation of Green's function matrix elements)
  std::valarray<double> G_up;		// Green's function for up spins
  dense_matrix Green0_up;       	// bath Green's function matrix for up spins
  dense_matrix Green_up;        	// Green's function matrix for up spins  
  std::valarray<double> G_down;		// Green's function for down spins
  dense_matrix Green0_down;	        // bath Green's function matrix for down spins
  dense_matrix Green_down;      	// Green's function matrix for down spins     
  std::vector<int> spins;		    // auxiliary Ising spins
  
};

#define ALPS_DMFT_HIRSCHFYESIM_H

#endif
