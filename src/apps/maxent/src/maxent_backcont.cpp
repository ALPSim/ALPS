/*
 * Copyright (C) 1998-2018 ALPS Collaboration.
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

// Adapted from CQMP/Maxent src/maxent_backcont.cpp at
// a3992534228b0de732e033ae23154969724814be; see ../README.md.
#include "maxent_backcont.hpp"

Backcont::Backcont(const MaxEntParameters *param_in) :
 param(param_in),k_type(param->getKernelType()){
} 

vector_type Backcont::backcontinue(const vector_type &A){
  const int ndat = param->ndat();
  const int nfreq = param->nfreq();

  G = vector_type(ndat);
  // G(X) =\int K(X,\omega)A(\omega)d\omega
  for(int n=0;n<ndat;n++){
    //recall that in MaxEntParameters we scaled the kernel by the 1/error
    double err = param->sigma(n);
    G(n)=0;
    for(int i=1;i<nfreq-1;i++){
      //G(n) += -1./M_PI*(param->K(n,i))*A(i)*(param->delta_omega(i)); 
      double delta = (param->omega_coord(i+1)) - (param->omega_coord(i-1));
      G(n) += err*(param->K(n,i))*A(i)*delta/2; 
    }

    G(n) += err*(param->K(n,0))*A(0)*(param->omega_coord(1)-param->omega_coord(0))/2;
    G(n) += err*(param->K(n,nfreq-1))*A(nfreq-1)*(param->omega_coord(nfreq-1)-param->omega_coord(nfreq-2))/2;
  } 

  return G;
}
double Backcont::max_error(const vector_type &y1, const vector_type &y2){
  double max_err = 0.0;
#ifdef NDEBUG
  if(y1.size() != y2.size()){
    throw std::runtime_error("backcont vector size mismatch!");
  }
#endif
  for(int i=0;i<y1.size();i++){
    double delta = std::abs(y1(i)-y2(i));
    if(delta>max_err){
      max_err = delta;
    }
  }
  return max_err;
}
