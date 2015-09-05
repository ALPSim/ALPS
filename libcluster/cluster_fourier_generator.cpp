/*****************************************************************************
 *
 * Copyright (C) 2005 - 2015 Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#include "boost/lexical_cast.hpp"

#include "cluster_fourier_transform.h"
#include "cluster.h"

void FourierTransformer::generate_cluster_transformer(const alps::params &parms,
                                              boost::shared_ptr<FourierTransformer> &fourier_ptr)
{
  int n_flavors = parms["FLAVORS"]|2;
  int n_site = parms["SITES"]|1;
  bool paramagnet = parms["PARAMAGNET"]|true;
  if (parms.defined("CLUSTER_LOOP")){
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
  else throw std::invalid_argument("cluster transformers should have CLUSTER_LOOP defined");
}




void FourierTransformer::generate_cluster_transformer_U(const alps::params &parms,
                                                boost::shared_ptr<FourierTransformer> &fourier_ptr,
                                                const std::vector<double> &densities)
{
  int n_flavors = parms["FLAVORS"]|2;
  int n_site = parms["SITES"]|1;
  double U = parms["U"];
  if (parms.defined("CLUSTER_LOOP")){
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
  else throw std::invalid_argument("cluster transformers should have CLUSTER_LOOP defined.");
}



void FourierTransformer::generate_cluster_transformer_U(const alps::params &parms,
                                                boost::shared_ptr<FourierTransformer> &fourier_ptr,
                                                const std::vector<double> &densities, const std::vector<double> &magnetization)
{
  int n_flavors = parms["FLAVORS"]|2;
  int n_site = parms["SITES"]|1;
  double U = parms["U"];
  if (parms.defined("CLUSTER_LOOP")){
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
  }
  else throw std::invalid_argument("cluster transformers should have CLUSTER_LOOP defined");
}
