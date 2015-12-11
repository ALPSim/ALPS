/*****************************************************************************
 *
 * Copyright (C) 2005 - 2015 Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#include "boost/lexical_cast.hpp"

#include "fourier_transform.h"

void FourierTransformer::generate_transformer(const alps::params &parms,
                                              boost::shared_ptr<FourierTransformer> &fourier_ptr)
{
  int n_flavors = parms["FLAVORS"]|2;
  int n_site = parms["SITES"]|1;
  bool paramagnet = parms["PARAMAGNET"]|true;
  if (parms.defined("GENERAL_FOURIER_TRANSFORMER")) {
    std::cout << "using general fourier transformer" << "\n";
    std::vector<double> eps(n_flavors);
    std::vector<double> epssq(n_flavors);
    for (int f=0; f<n_flavors; ++f) {
      eps[f] = parms["EPS_"+boost::lexical_cast<std::string>(f)];
      epssq[f] = parms["EPSSQ_"+boost::lexical_cast<std::string>(f)];
    }
    fourier_ptr.reset(new SimpleG0FourierTransformer((double)parms["BETA"], (double)parms["MU"], (double)parms["H"],
                                                     n_flavors, eps, epssq));
  }
  else if (parms.defined("CLUSTER_LOOP")){
    throw std::logic_error("please use a cluster fourier transformer, see libcluster");
  }
  else {
    std::cout << "using Bethe lattice fourier transform" << "\n";
    std::vector<double> eps(n_flavors);
    std::vector<double> epssq(n_flavors);
    double t = parms["t"];
    for (int f=0; f<n_flavors; ++f) {
      eps[f] = 0.;
      epssq[f] = t*t;
    }
    fourier_ptr.reset(new SimpleG0FourierTransformer((double)parms["BETA"], (double)parms["MU"], (double)parms["H"],
                                                     n_flavors, eps, epssq));
  }
}




void FourierTransformer::generate_transformer_U(const alps::params &parms,
                                                boost::shared_ptr<FourierTransformer> &fourier_ptr,
                                                const std::vector<double> &densities)
{
  int n_flavors = parms["FLAVORS"]|2;
  int n_site = parms["SITES"]|1;
  double U = parms["U"];
  if (parms.defined("GENERAL_FOURIER_TRANSFORMER")) {
    std::cout << "using general fourier transformer" << "\n";
    std::vector<std::vector<double> > eps(n_flavors);
    std::vector<std::vector<double> >epssq(n_flavors);
    for (int f=0; f<n_flavors; ++f) {
      eps[f].resize(1);
      epssq[f].resize(1);
      eps[f][0] = parms["EPS_"+boost::lexical_cast<std::string>(f)];
      epssq[f][0] = parms["EPSSQ_"+boost::lexical_cast<std::string>(f)];
    }
    fourier_ptr.reset(new GFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., U,
                                              n_flavors, 1, densities, eps, epssq));
  }
  else if (parms.defined("CLUSTER_LOOP")){
    throw std::logic_error("please use a cluster fourier transformer, see libcluster");
  }
  else {
    std::cout << "using Bethe lattice fourier transform" << "\n";
    std::vector<std::vector<double> > eps(n_flavors);
    std::vector<std::vector<double> >epssq(n_flavors);
    double t = parms["t"];
    for (int f=0; f<n_flavors; ++f) {
      eps[f].resize(1);
      epssq[f].resize(1);
      eps[f][0] = 0.;
      epssq[f][0] = t*t;
    }
    fourier_ptr.reset(new GFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., U,
                                              n_flavors, 1, densities, eps, epssq));
  }
}



void FourierTransformer::generate_transformer_U(const alps::params &parms,
                                                boost::shared_ptr<FourierTransformer> &fourier_ptr,
                                                const std::vector<double> &densities, const std::vector<double> &magnetization)
{
  int n_flavors = parms["FLAVORS"]|2;
  int n_site = parms["SITES"]|1;
  double U = parms["U"];
  if (parms.defined("GENERAL_FOURIER_TRANSFORMER")) {
    std::cout << "using general fourier transformer" << "\n";
    std::vector<double> eps(n_flavors);
    std::vector<double> epssq(n_flavors);
    for (int f=0; f<n_flavors; ++f) {
      eps[f] = parms["EPS_"+boost::lexical_cast<std::string>(f)];
      epssq[f] = parms["EPSSQ_"+boost::lexical_cast<std::string>(f)];
    }
    fourier_ptr.reset(new SimpleAntiFerroGFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., (double)parms["H"], U,
                                                             n_flavors, densities, magnetization,  eps, epssq));
  }
  else if (parms.defined("CLUSTER_LOOP")){
    throw std::logic_error("please use a cluster fourier transformer, see libcluster");
  } else {
    std::cout << "using AFM Bethe lattice fourier transform" << "\n";
    std::vector<double> eps(n_flavors);
    std::vector<double> epssq(n_flavors);
    double t = parms["t"];
    for (int f=0; f<n_flavors; ++f) {
      eps[f] = 0.;
      epssq[f] = t*t;
    }
    fourier_ptr.reset(new SimpleAntiFerroGFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., (double)parms["H"], U,
                                                             n_flavors, densities, magnetization, eps, epssq));
  }
}
