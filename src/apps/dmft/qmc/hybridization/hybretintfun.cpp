/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                   
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include"hybretintfun.hpp"
#include "../retarded_interaction.h"

#include <algorithm>
//construct a retarded interaction function. ntime: number of time slices.
//noffidag_orbitals: number of offdiagonal orbitals. ndiag_orbitals: number
//of diagonal orbitals. 
ret_int_fun::ret_int_fun(const alps::input::hybridization_parameters &input):
green_function<double>(input.n_tau+1, 1, 2)//2 "flavors": K and K'
{
  bool use_retarded_interaction=!input.ret_int_k.empty();
  if(!use_retarded_interaction) return;

  beta_=input.beta;
  
  //read in Green's function from a file
  read_interaction_K_function(input);
  interaction_K_function_sanity_check();

  extern int global_mpi_rank;
  if(global_mpi_rank==0){
    std::cout<<*this<<std::endl;
  }
}

void ret_int_fun::interaction_K_function_sanity_check(void){
  // read_retarded_interaction checks the length, K >= 0 and K(0) = K(beta) = 0.
}

//this routine reads in the retarded interaction function, either from a text file or from an hdf5 file (for easy passing of binary data).
//In  case of text files the file format is index - hyb_1 - hyb2 - hyb3 - ... in columns that go from time=0 to time=beta. Note that
//the retarded interaction function is in imaginary time and always positive both for negative and positive times. It is also symmetric.
void ret_int_fun::read_interaction_K_function(const alps::input::hybridization_parameters &input){
  if(input.ret_int_k.empty()) throw(std::invalid_argument(std::string("Parameter RET_INT_K missing, filename for retarded interaction function not specified.")));
  const retarded_interaction r = read_retarded_interaction(input.ret_int_k, input.k_in_hdf5, ntime()-1);
  for(std::size_t i=0; i<ntime(); i++){
    operator()(i,0)=r.K[i];
    operator()(i,1)=r.Kp[i];
  }
}

std::ostream &operator<<(std::ostream &os, const ret_int_fun &K){
  os<<"the retarded interaction function and derivative are: "<<std::endl;
  for(int i=0;i<10;++i){ std::cout<<i<<" "; for(std::size_t j=0;j<K.nflavor();++j){ std::cout<<K(i,j)<<" ";} std::cout<<std::endl; }
  os<<"... *** etc *** ...\n";
  os<<K.ntime()-1<<" "; for(std::size_t j=0;j<K.nflavor();++j){ std::cout<<K(K.ntime()-1,j)<<" ";} std::cout<<std::endl;
  return os;
}
//linear interpolation of the retarded interaction function 
double ret_int_fun::interpolate(double time) const{
  time=std::abs(time);

  double n = time/beta_*(ntime()-1);
  // at time = beta, interpolate on the last interval instead of reading past the end
  int n_lower = std::min((int)n, (int)ntime()-2); // interpolate linearly between n_lower and n_lower+1

  return operator()(n_lower,0) + (n-n_lower)*(operator()(n_lower+1,0)-operator()(n_lower,0));
}

//linear interpolation of the retarded interaction function.
double ret_int_fun::interpolate_deriv(double time) const{
  //if(time<-beta_ || time > beta_) std::cout<< "!!" << std::endl;
  if(time<0.) return -1.*interpolate_deriv(-time);

  double n = time/beta_*(ntime()-1);
  // at time = beta, interpolate on the last interval instead of reading past the end
  int n_lower = std::min((int)n, (int)ntime()-2); // interpolate linearly between n_lower and n_lower+1

  return operator()(n_lower,1) + (n-n_lower)*(operator()(n_lower+1,1)-operator()(n_lower,1));
}
