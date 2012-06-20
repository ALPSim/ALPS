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
#include "xml.h"
#include <boost/filesystem.hpp>

typedef alps::SignedObservable<alps::SimpleObservable<double,alps::DetailedBinning<double> > > signed_obs_t;
typedef alps::SignedObservable<alps::RealVectorObservable> signed_vec_obs_t;
typedef alps::RealVectorObservable vec_obs_t;
typedef alps::SimpleObservable<double,alps::DetailedBinning<double> > simple_obs_t;
typedef const alps::SimpleObservable<double,alps::DetailedBinning<double> > const_simple_obs_t;



///This function is used in connection with reduce matsubara frequencies to
///store the mean and the error of a function in Matsubara frequencies.
///Call this to get the result of a InteractionExpansion simulation on all nodes. 
std::pair<matsubara_green_function_t,itime_green_function_t> InteractionExpansionSim::get_result()
{
  std::istringstream in_omega(parms["G0(omega)"]);
  clock_t time0=clock();
  if(parms.defined("CHECKPOINT")){ //write checkpoint
    alps::Parameters *p=const_cast<alps::Parameters*>(&parms); //remove G0 from parameters as our xml parser cannot handle nested tags.
    p->erase(std::string("G0(omega)"));
    alps::Parameters::iterator it=p->begin();
    while(it!=p->end()){
      if(strncmp(it->key().c_str(),"G0(omega)",9)==0){
        it->value()="...ignored, nested...";
      }
      it++;
    }
    std::cout<<"end of sim, checkpointing"<<std::endl;
    std::string fns=parms["CHECKPOINT"];
    fns+=".xml";
    boost::filesystem::path fn(fns);
    checkpoint(fn);
  }
  clock_t time1=clock();
  std::cout<<"time for writing checkpoint was: "<<(time1-time0)/(double)CLOCKS_PER_SEC<<std::endl;
  ///end checkpoint
  std::cout<<"collecting results."<<std::endl;
  unsigned int n_matsubara=p_["NMATSUBARA"];
  unsigned int n_matsubara_measurements=p_.value_or_default("NMATSUBARA_MEASUREMENTS", n_matsubara);
  unsigned int n_tau=p_["N"];
  unsigned int n_self=p_.value_or_default("NSELF", 10*n_tau);
  spin_t n_flavors(p_.value_or_default("FLAVORS",2));
  unsigned int n_site(p_.value_or_default("SITES",1));
  double beta(p_["BETA"]);
  itime_green_function_t green_itime_measured(n_tau+1, n_site, n_flavors);
  matsubara_green_function_t green_matsubara_measured(n_matsubara, n_site, n_flavors);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  boost::shared_ptr<FourierTransformer> fourier_ptr_g0;
  FourierTransformer::generate_transformer(p_, fourier_ptr_g0);
  //find whether our data is in imaginary time or frequency:
  bool measure_in_matsubara=true;
  if(p_.value_or_default("HISTOGRAM_MEASUREMENT", false)) 
    measure_in_matsubara=false;
  
  const bool compactit = parms.value_or_default("GET_COMPACTED_MEASUREMENTS", true);
  if (!compactit) 
    std::cout << "getting non-compacted measurements\n";
  alps::ObservableSet gathered_measurements=get_measurements(compactit);
  alps::RealVectorObsevaluator pert_obseval=gathered_measurements["PertOrder"];
  alps::RealObsevaluator vertex_insertion=gathered_measurements["VertexInsertion"];
  alps::RealObsevaluator vertex_removal=gathered_measurements["VertexRemoval"];
  if(vertex_insertion.count()!=0) 
    std::cout<<"vertex insertion: attempts: "<<vertex_insertion.count()<<" acceptance: "<<vertex_insertion.mean()<<std::endl;
  if(vertex_removal  .count()!=0) 
    std::cout<<"vertex removal  : attempts: "<<vertex_removal  .count()<<" acceptance: "<<vertex_removal  .mean()<<std::endl;
  alps::RealObsevaluator sign_obseval=gathered_measurements["Sign"];
  std::valarray<double> mean_order(pert_obseval.mean());
  std::cout<<"average matrix size was: "<<std::endl;
  std::ofstream matrix_size("matrix_size", std::ios::app);
  for(unsigned int i=0;i<n_flavors;++i){
    std::cout<<mean_order[i]<<"\t";
    matrix_size<<mean_order[i]<<"\t";
  }
  std::cout<<std::endl;
  matrix_size<<std::endl;
  std::cout<<"average sign was: "<<sign_obseval.mean()<<"+-"<<sign_obseval.error()<<std::endl;
  //single particle Green function measurements
  matsubara_green_function_t bare_green_matsubara(n_matsubara, n_site, n_flavors);
  std::vector<double> densities(n_flavors);
  read_freq(in_omega, bare_green_matsubara);
  if(measure_in_matsubara) {
    std::cout<<"evaluating Matsubara data"<<std::endl;
    evaluate_selfenergy_measurement_matsubara(gathered_measurements, green_matsubara_measured, 
                                              bare_green_matsubara, densities, 
                                              beta, n_site, n_flavors, n_matsubara_measurements);
  } 
  else {
    std::cout<<"evaluating imaginary time data"<<std::endl;
    itime_green_function_t bare_green_itime(n_tau+1, n_site, n_flavors);
    fourier_ptr_g0->backward_ft(bare_green_itime, bare_green_matsubara);
    evaluate_selfenergy_measurement_itime_rs(gathered_measurements, green_itime_measured, bare_green_itime, 
                                             beta, n_site, n_flavors, n_tau, n_self);
  }
  //Fourier transformations
  if (!measure_in_matsubara) {
    for (unsigned int z=0; z<n_flavors; ++z) {
      densities[z] = 0;
      for (unsigned int i=0; i<n_site; ++i)
        densities[z] -= green_itime_measured(n_tau,i,i,z);
      densities[z] /= n_site;
    }
  }
  FourierTransformer::generate_transformer_U(p_, fourier_ptr, densities);
  if (measure_in_matsubara) {
    fourier_ptr->append_tail(green_matsubara_measured, bare_green_matsubara, n_matsubara_measurements);
    fourier_ptr->backward_ft(green_itime_measured, green_matsubara_measured);
  }
  else 
    fourier_ptr->forward_ft(green_itime_measured, green_matsubara_measured);
  std::cout<<"returning to sc loop"<<std::endl;
  return std::make_pair(green_matsubara_measured, green_itime_measured);
} 
