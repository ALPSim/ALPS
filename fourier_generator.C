/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
 *                       Matthias Troyer <troyer@comp-phys.org>
 *                       Emanuel Gull <gullc@comp-phys.org>
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

#include "fouriertransform.h"
#include "cluster.h"
#include "boost/lexical_cast.hpp"

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
    std::string clusterloop=parms["CLUSTER_LOOP"];
    if(clusterloop==std::string("DCA")){
      ClusterTransformer* clusterhandler=ClusterTransformer::generate_transformer(parms);
      std::vector<std::vector<double> > eps_f(n_flavors, clusterhandler->epsav());
      std::vector<std::vector<double> > epssq_f(n_flavors, clusterhandler->epssqav());
      fourier_ptr.reset(new ClusterG0FourierTransformer((double)parms["BETA"], (double)parms["MU"], (double)parms["H"],
                                                        n_flavors, n_site, eps_f, epssq_f));
      delete clusterhandler;
    }else if(clusterloop==std::string("CDMFT")){
      std::cout<<"using CDMFT cluster Fourier transformer"<<std::endl;
      std::cout<<"warning: G0 moments are not implemented for CDMFT"<<std::endl;
      std::vector<std::vector<double> > eps_f(n_flavors, std::vector<double>(n_site*n_site,0.00));
      std::vector<std::vector<double> > epssq_f(n_flavors, std::vector<double>(n_site*n_site,0.0001));
      fourier_ptr.reset(new ClusterG0FourierTransformer((double)parms["BETA"], (double)parms["MU"], (double)parms["H"],
                                                        n_flavors, n_site, eps_f, epssq_f));
    }
    else throw std::logic_error("please figure out how to do fourier transforms in your cluster framework if it is not DCA: "+clusterloop);
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
    std::string clusterloop=parms["CLUSTER_LOOP"];
    if(clusterloop==std::string("DCA")){
      ClusterTransformer* clusterhandler=ClusterTransformer::generate_transformer(parms);
      std::vector<std::vector<double> > eps_f(n_flavors, clusterhandler->epsav());
      std::vector<std::vector<double> > epssq_f(n_flavors, clusterhandler->epssqav());
      fourier_ptr.reset(new GFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., U,
                                                n_flavors, n_site, densities, eps_f, epssq_f));
      delete clusterhandler;
    }else if(clusterloop==std::string("CDMFT")){
      std::cout<<"using CDMFT cluster Fourier transformer"<<std::endl;
      std::cout<<"warning: G0 moments are not implemented for CDMFT"<<std::endl;
      std::vector<std::vector<double> > eps_f(n_flavors, std::vector<double>(n_site*n_site,0));
      std::vector<std::vector<double> > epssq_f(n_flavors, std::vector<double>(n_site*n_site,0));
      fourier_ptr.reset(new ClusterG0FourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2, U,
                                                        n_flavors, n_site, eps_f, epssq_f));
    }
    else throw std::logic_error("please figure out how to do fourier transforms in your cluster framwork if it is not DCA: "+clusterloop);
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
    std::string clusterloop=parms["CLUSTER_LOOP"];
    if(clusterloop==std::string("DCA")){
      std::cout << "using cluster fourier transformer" << "\n";
      ClusterTransformer* clusterhandler=ClusterTransformer::generate_transformer(parms);
      std::vector<std::vector<double> > eps_f(n_flavors, clusterhandler->epsav());
      std::vector<std::vector<double> > epssq_f(n_flavors, clusterhandler->epssqav());
      fourier_ptr.reset(new ClusterAntiFerroGFourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2., (double)parms["H"], U,
                                                                n_flavors, n_site, densities, magnetization, eps_f, epssq_f));
      delete clusterhandler;
    }else if(clusterloop==std::string("CDMFT")){
      std::cout<<"using CDMFT cluster Fourier transformer"<<std::endl;
      std::cout<<"warning: G0 moments are not implemented for CDMFT"<<std::endl;
      std::vector<std::vector<double> > eps_f(n_flavors, std::vector<double>(n_site*n_site,0));
      std::vector<std::vector<double> > epssq_f(n_flavors, std::vector<double>(n_site*n_site,0));
      fourier_ptr.reset(new ClusterG0FourierTransformer((double)parms["BETA"], (double)parms["MU"]+U/2, U,
                                                        n_flavors, n_site, eps_f, epssq_f));
    }
    else throw std::logic_error("please figure out how to do fourier transforms in your cluster framwork if it is not DCA (1)");
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
