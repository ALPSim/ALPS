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

#include "interaction_expansion.hpp"
#include <complex>
#include <alps/alea.h>
#include <alps/alea/simpleobseval.h>
#include <alps/scheduler/montecarlo.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>



typedef alps::SignedObservable<alps::SimpleObservable<double,alps::DetailedBinning<double> > > signed_obs_t;
typedef alps::SignedObservable<alps::RealVectorObservable> signed_vec_obs_t;
typedef alps::RealVectorObservable vec_obs_t;
typedef alps::SimpleObservable<double,alps::DetailedBinning<double> > simple_obs_t;
typedef const alps::SimpleObservable<double,alps::DetailedBinning<double> > const_simple_obs_t;




///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables. It is also called at the start of every DMFT
///iteration.
void InteractionExpansionRun::initialize_observables(void) 
{
  if(measurements.has("Sign")){
    measurements.clear();
  }
  measurements<<alps::RealObservable("Sign");
  measurements<<alps::RealVectorObservable("PertOrder");  
  if(measurement_method==selfenergy_measurement_itime_rs) {
    for(unsigned int flavor=0;flavor<n_flavors;++flavor){
      for(unsigned int i=0;i<n_site;++i){
        for(unsigned int j=0;j<n_site;++j){
          std::stringstream obs_name;
          obs_name<<"W_"<<flavor<<"_"<<i<<"_"<<j;
          measurements<<signed_vec_obs_t(obs_name.str().c_str());
        }
      }
    }
  }
  else {
    for(unsigned int flavor=0;flavor<n_flavors;++flavor){
      for (unsigned int k=0; k<n_site; k++) {                   
        std::stringstream obs_name_real, obs_name_imag;
        obs_name_real<<"Wk_real_"<<flavor<<"_"<<k << "_" << k;
        obs_name_imag<<"Wk_imag_"<<flavor<<"_"<<k << "_" << k;
        measurements<<signed_vec_obs_t(obs_name_real.str().c_str());
        measurements<<signed_vec_obs_t(obs_name_imag.str().c_str());
      }
    }
  }
  measurements << signed_vec_obs_t("densities");
  for(unsigned int flavor=0;flavor<n_flavors;++flavor)
    measurements << signed_vec_obs_t("densities_"+boost::lexical_cast<std::string>(flavor));
  measurements << signed_obs_t("density_correlation");
  measurements << signed_vec_obs_t("n_i n_j");
  for(unsigned int flavor=0;flavor<n_flavors;++flavor){
    for(unsigned int i=0;i<n_site;++i){
      std::stringstream density_name, sz_name;
      density_name<<"density_"<<flavor<<"_"<<i;
      measurements<<signed_obs_t(density_name.str().c_str());
    }
  }
  for(unsigned int i=0;i<n_site;++i){
    std::stringstream sz_name, sz2_name, sz0_szj_name;
    sz_name<<"Sz_"<<i;
    sz2_name<<"Sz2_"<<i;
    sz0_szj_name<<"Sz0_Sz"<<i;
    measurements<<signed_obs_t(sz_name.str().c_str());
    measurements<<signed_obs_t(sz2_name.str().c_str());
    measurements<<signed_obs_t(sz0_szj_name.str().c_str());
  }
  //acceptance probabilities
  measurements<<alps::RealObservable("VertexInsertion");
  measurements<<alps::RealObservable("VertexRemoval");
  measurements.reset(true);
}




///this function is called whenever measurements should be performed. Depending
///on the value of  measurement_method it will choose one particular
///measurement function. 
void InteractionExpansionRun::measure_observables() 
{
  //measure the fermionic sign
  measurements["Sign"]<<sign;
  //measure the expansion order (perturbation order)
  std::valarray<double> pert_order(n_flavors);
  for(unsigned int i=0;i<n_flavors;++i) { 
      assert(num_rows(M[i].matrix()) == num_cols(M[i].matrix()));
    pert_order[i]=num_rows(M[i].matrix());
  }
  measurements["PertOrder"] << pert_order;

  //measure observables in frequency space
  if (measurement_method == selfenergy_measurement_matsubara){
    compute_W_matsubara();
  }
  //measure observables in imaginary time space
  else if (measurement_method == selfenergy_measurement_itime_rs){
    compute_W_itime();
  }
}


