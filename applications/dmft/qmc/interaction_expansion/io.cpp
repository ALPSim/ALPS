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
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/

#include "interaction_expansion.hpp"

///read in Green's function file. Format:
/// Frequency \t val \t val.... 
void InteractionExpansionRun::read_bare_green(std::ifstream &G0_omega, std::ifstream &G0_tau)
{
  assert(G0_omega.is_open() && G0_tau.is_open());
  double ignored;
  for(frequency_t o=0;o<n_matsubara;++o){
    G0_omega >>ignored;
    for(spin_t flavor=0;flavor<n_flavors;++flavor){
      for(site_t j=0;j<n_site;++j){
        for(site_t k=0;k<n_site;++k){
          //initialize bare Green's function
          G0_omega >> bare_green_matsubara(o,j,k, flavor);
        }
      }
    }
  }
  green_matsubara=bare_green_matsubara; //starting values for the dressed GF
  for(itime_index_t tau=0;tau<=n_tau;++tau){
    G0_tau>>ignored;
    for(spin_t flavor=0;flavor<n_flavors;++flavor){
      for(site_t j=0;j<n_site;++j){
        for(site_t k=0;k<n_site;++k){
          G0_tau>>bare_green_itime(tau, j, k, flavor);
        }
      }
    }
  }
  green_itime=bare_green_itime; //starting values for the dressed green's function
  std::cout<<"G0_omega: "<<G0_omega<<std::endl;
  std::cout<<"G0_tau: "<<G0_tau<<std::endl;
}


void InteractionExpansionRun::print(std::ostream &os){
  os<<"InteractionExpansion solver, parameters: "<<std::endl;
  os<<"max order\t"<<max_order<<"\tn_flavors: "
    <<n_flavors<<"\tn_site: "<<n_site
    <<"\tn_matsubara"<<n_matsubara<<std::endl;
  os<<"n_tau: "<<n_tau<<"\tmc steps: "<<mc_steps
    <<"\ttherm steps: "<<therm_steps<<std::endl;
  
  os<<"beta: "<<beta<<"\talpha: "<<alpha<<"\tU: "<<onsite_U<<std::endl;
  
  os<<"recalc period: "<<recalc_period<<"\tmeasurement period: "<< measurement_period
    <<"\tconvergence period: "<< convergence_check_period<<std::endl;
  os<<"almost zero: "<<almost_zero<<std::endl;
}
