/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 *                       Emanuel Gull <gull@pks.mpg.de>,
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/
#ifndef HYB_EVALUATE
#define HYB_EVALUATE
#include <boost/math/special_functions/bessel.hpp>

void evaluate_basics(const alps::results_type<hybridization>::type &results, const alps::parameters_type<hybridization>::type &parms, alps::hdf5::archive &solver_output);
void evaluate_time(const alps::results_type<hybridization>::type &results, const alps::parameters_type<hybridization>::type &parms, alps::hdf5::archive &solver_output);
void evaluate_freq(const alps::results_type<hybridization>::type &results, const alps::parameters_type<hybridization>::type &parms, alps::hdf5::archive &solver_output);
void evaluate_legendre(const alps::results_type<hybridization>::type &results, const alps::parameters_type<hybridization>::type &parms, alps::hdf5::archive &solver_output);
void evaluate_nnt(const alps::results_type<hybridization>::type &results, const alps::parameters_type<hybridization>::type &parms, alps::hdf5::archive &solver_output);
void evaluate_nnw(const alps::results_type<hybridization>::type &results, const alps::parameters_type<hybridization>::type &parms, alps::hdf5::archive &solver_output);
void evaluate_sector_statistics(const alps::results_type<hybridization>::type &results, const alps::parameters_type<hybridization>::type &parms, alps::hdf5::archive &solver_output);
void evaluate_2p(const alps::results_type<hybridization>::type &results, const alps::parameters_type<hybridization>::type &parms, alps::hdf5::archive &solver_output);

inline std::complex<double> t(int n, int l){//transformation matrix from Legendre to Matsubara basis
  std::complex<double> i_c(0., 1.);
  return (std::sqrt(static_cast<double>(2*l+1))/std::sqrt(static_cast<double>(2*n+1))) * std::exp(i_c*(n+0.5)*M_PI) * std::pow(i_c,l) * boost::math::cyl_bessel_j(l+0.5,(n+0.5)*M_PI);
}

struct jackson{//the Jackson kernel
  jackson(int N_l_):N_l(N_l_){};
  double operator()(int n){ return (1./(N_l+1.))*((N_l-n+1)*std::cos(M_PI*n/(N_l+1.))+std::sin(M_PI*n/(N_l+1.))/std::tan(M_PI/(N_l+1.))); }
  int N_l;
};

#endif

