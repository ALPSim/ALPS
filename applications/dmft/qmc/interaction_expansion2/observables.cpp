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
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
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
void InteractionExpansion::initialize_observables(void) 
{
  if(measurements.has("Sign")){
    measurements.clear();
  }
#ifdef ALPS_NGS_USE_NEW_ALEA
  measurements << alps::accumulator::RealObservable("Sign");
  measurements << alps::accumulator::RealVectorObservable("PertOrder");
#else
  measurements << alps::ngs::RealObservable("Sign");
  measurements << alps::ngs::RealVectorObservable("PertOrder");
#endif
  if(measurement_method==selfenergy_measurement_itime_rs) {
    for(unsigned int flavor=0;flavor<n_flavors;++flavor){
      for(unsigned int i=0;i<n_site;++i){
        for(unsigned int j=0;j<n_site;++j){
          std::stringstream obs_name;
          obs_name<<"W_"<<flavor<<"_"<<i<<"_"<<j;
#ifndef ALPS_NGS_USE_NEW_ALEA
          measurements << alps::ngs::RealVectorObservable(obs_name.str().c_str());
#else
          throw std::runtime_error("alps::ngs::SignedRealVectorObservable is not implemented");
#endif //ALPS_NGS_USE_NEW_ALEA
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
#ifndef ALPS_NGS_USE_NEW_ALEA
        measurements << alps::ngs::RealVectorObservable(obs_name_real.str().c_str());
        measurements << alps::ngs::RealVectorObservable(obs_name_imag.str().c_str());
#else
        throw std::runtime_error("alps::ngs::SignedRealVectorObservable is not implemented");
#endif //ALPS_NGS_USE_NEW_ALEA
      }
    }
  }
#ifndef ALPS_NGS_USE_NEW_ALEA
  measurements << alps::ngs::RealVectorObservable("densities");
  for(unsigned int flavor=0;flavor<n_flavors;++flavor)
    measurements << alps::ngs::RealVectorObservable("densities_"+boost::lexical_cast<std::string>(flavor));
  measurements << alps::ngs::RealObservable("density_correlation");
  measurements << alps::ngs::RealVectorObservable("n_i n_j");
#else
  throw std::runtime_error("alps::ngs::SignedRealVectorObservable is not implemented");
#endif //ALPS_NGS_USE_NEW_ALEA
  for(unsigned int flavor=0;flavor<n_flavors;++flavor){
    for(unsigned int i=0;i<n_site;++i){
      std::stringstream density_name, sz_name;
      density_name<<"density_"<<flavor;
      if (n_site>1) density_name<<"_"<<i;
#ifndef ALPS_NGS_USE_NEW_ALEA
      measurements << alps::ngs::RealObservable(density_name.str().c_str());
#else
  throw std::runtime_error("alps::ngs::SignedRealVectorObservable is not implemented");
#endif //ALPS_NGS_USE_NEW_ALEA
    }
  }
  for(unsigned int i=0;i<n_site;++i){
    std::stringstream sz_name, sz2_name, sz0_szj_name;
    sz_name<<"Sz_"<<i;
    sz2_name<<"Sz2_"<<i;
    sz0_szj_name<<"Sz0_Sz"<<i;
// #ifndef ALPS_NGS_USE_NEW_ALEA
//     measurements << alps::ngs::SignedRealObservable(sz_name.str().c_str());
//     measurements << alps::ngs::SignedRealObservable(sz2_name.str().c_str());
//     measurements << alps::ngs::SignedRealObservable(sz0_szj_name.str().c_str());
// #else
  throw std::runtime_error("alps::ngs::SignedRealVectorObservable is not implemented");
// #endif //ALPS_NGS_USE_NEW_ALEA
  }
  //acceptance probabilities
#ifdef ALPS_NGS_USE_NEW_ALEA
  measurements << alps::accumulator::RealObservable("VertexInsertion");
  measurements << alps::accumulator::RealObservable("VertexRemoval");
  measurements << alps::accumulator::RealObservable("MeasurementTime");
  measurements << alps::accumulator::RealObservable("UpdateTime");
  measurements << alps::accumulator::RealObservable("RecomputeTime");
#else
  measurements << alps::ngs::RealObservable("VertexInsertion");
  measurements << alps::ngs::RealObservable("VertexRemoval");
  measurements << alps::ngs::RealObservable("MeasurementTime");
  measurements << alps::ngs::RealObservable("UpdateTime");
  measurements << alps::ngs::RealObservable("RecomputeTime");
#endif
  measurements.reset(true);
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
  measurements["PertOrder"] << pert_order;
}
