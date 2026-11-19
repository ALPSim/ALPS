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
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "interaction_expansion.hpp"
#include "sz_results.hpp"
#include <complex>
#include <alps/scheduler/montecarlo.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>





///this function is called at the start of the simulation for allocation of
///memory for the ALPS observables. It is also called at the start of every DMFT
///iteration.
void InteractionExpansion::initialize_observables(void) 
{
  if(measurements.has("Sign")){
    measurements = alps::alea::observable_set();
  }
  measurements << alps::alea::observable("Sign");
  measurements << alps::alea::observable("PertOrder");
  if(measurement_method==selfenergy_measurement_itime_rs) {
    for(unsigned int flavor=0;flavor<n_flavors;++flavor){
      for(unsigned int i=0;i<n_site;++i){
        for(unsigned int j=0;j<n_site;++j){
          std::stringstream obs_name;
          obs_name<<"W_"<<flavor<<"_"<<i<<"_"<<j;
          measurements << alps::alea::observable(obs_name.str().c_str());
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
        measurements << alps::alea::observable(obs_name_real.str().c_str());
        measurements << alps::alea::observable(obs_name_imag.str().c_str());
      }
    }
  }
  measurements << alps::alea::observable("densities");
  for(unsigned int flavor=0;flavor<n_flavors;++flavor)
    measurements << alps::alea::observable("densities_"+boost::lexical_cast<std::string>(flavor));
  measurements << alps::alea::observable("density_correlation");
  measurements << alps::alea::observable("n_i n_j");
  for(unsigned int flavor=0;flavor<n_flavors;++flavor){
    for(unsigned int i=0;i<n_site;++i){
      std::stringstream density_name, sz_name;
      density_name<<"density_"<<flavor;
      if (n_site>1) density_name<<"_"<<i;
      measurements << alps::alea::observable(density_name.str().c_str());
    }
  }
  for(unsigned int i=0;i<n_site;++i){
    std::stringstream sz_name, sz2_name, sz0_szj_name;
    sz_name<<"Sz_"<<i;
    sz2_name<<"Sz2_"<<i;
    sz0_szj_name<<"Sz0_Sz"<<i;
    // Measured as value times sign. NGS's signed observables divided their
    // mean by the mean sign; alps::solvers::ctint divides the results by the
    // "Sign" result with divide_sz_results_by_sign before saving them.
    measurements << alps::alea::observable(sz_name.str().c_str());
    measurements << alps::alea::observable(sz2_name.str().c_str());
    measurements << alps::alea::observable(sz0_szj_name.str().c_str());
  }
  //acceptance probabilities
  measurements << alps::alea::observable("VertexInsertion");
  measurements << alps::alea::observable("VertexRemoval");
  measurements << alps::alea::observable("MeasurementTime");
  measurements << alps::alea::observable("UpdateTime");
  measurements << alps::alea::observable("RecomputeTime");
  measurements.reset();
}


void divide_sz_results_by_sign(alps::alea::result_set & results, unsigned int n_site)
{
  if(!results.has("Sign"))
    return;
  alps::alea::observable_result const sign(results["Sign"]);
  for(unsigned int i=0;i<n_site;++i){
    std::stringstream sz_name, sz2_name, sz0_szj_name;
    sz_name<<"Sz_"<<i;
    sz2_name<<"Sz2_"<<i;
    sz0_szj_name<<"Sz0_Sz"<<i;
    std::string const names[] = {sz_name.str(), sz2_name.str(), sz0_szj_name.str()};
    // Only the imaginary-time measurement fills them, and the MPI adapter
    // leaves out observables that no rank measured.
    for(std::string const & name : names)
      if(results.has(name) && results[name].count())
        results.insert(name, results[name] / sign);
  }
}




///this function is called whenever measurements should be performed. Depending
///on the value of  measurement_method it will choose one particular
///measurement function. 
void InteractionExpansion::measure_observables(void) 
{
  measurements["Sign"]<<sign;
  if (measurement_method == selfenergy_measurement_matsubara)
    compute_W_matsubara();
  else if (measurement_method == selfenergy_measurement_itime_rs)
    compute_W_itime();
  std::valarray<double> pert_order(n_flavors);
  for(unsigned int i=0;i<n_flavors;++i) { 
      assert(num_rows(M[i].matrix()) == num_cols(M[i].matrix()));
    pert_order[i]=num_rows(M[i].matrix());
  }
  measurements["PertOrder"] << alps::numeric::valarray2vector(pert_order);
}
