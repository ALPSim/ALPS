/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2007-2010 by Emanuel Gull <gull@phys.columbia.edu>
 * Copyright (C) 2007-2008 by Sebastian Fuchs <fuchs@comp-phys.org>
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


#include "cluster.h"
#include "alps/alea.h"
#include <boost/numeric/ublas/io.hpp>

///self-consistency transform operator for the Antiferromagnetic DCA case
matsubara_green_function_t AntiFerroDCATransformer::operator()(const matsubara_green_function_t &G_omega,
                                                               const matsubara_green_function_t &G0_omega,
                                                               double &mu, const double h, const double beta)
{
  //std::cout<<"runnning operator() with AFM"<<std::endl;
  unsigned int nflavor=G0_omega.nflavor();
  unsigned int nfreq=G0_omega.nfreq();
  
  std::cout<<"using a magnetic field of: "<<h<<std::endl;
  
  //first: compute the self energy using G and G0.
  matsubara_green_function_t G0_omega_new(nfreq,n_cluster_sites(),nflavor, blockdiagonal);
  G0_omega_new.clear();
  matsubara_green_function_t G0_omega_inv = green_function_inverse (G0_omega);
  matsubara_green_function_t G_omega_inv  = green_function_inverse (G_omega );
  matsubara_green_function_t selfenergy   = green_function_subtract(G0_omega_inv, G_omega_inv);
  
  //buffer the dispersion and the weight – no need to evaluate this for every frequency
  //first dimension: sign
  //second dimension: k0 or k1
  //third dimension: cluster momentum kc
  //fourth dimension: lattice momentum kl
  std::vector<std::vector<double> >epsilon_kckl(n_cluster_sites(), std::vector<double>(n_lattice_momenta()));
  for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
    for (site_t kl=0; kl<n_lattice_momenta(); ++kl) {
      vector_type k(dimension());
      for (unsigned int d=0; d<dimension(); ++d){
        k[d] = lattice_momentum(kl,d) + cluster_momentum(kc  ,d);
      }
      epsilon_kckl[kc][kl]=epsilon(k);
    }
  }
  
  for (site_t kl=0; kl<n_cluster_cells(); ++kl) {
    matsubara_green_function_t G_kl(G0_omega);
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      for (site_t kc=0; kc<n_cluster_sites(); kc+=2) {
        G_kl(freq,kc  ,kc  ,0) = iomega_n+mu-epsilon_kckl[kc  ][kl];
        G_kl(freq,kc  ,kc+1,0) = -h/2;
        G_kl(freq,kc+1,kc  ,0) = -h/2;
        G_kl(freq,kc+1,kc+1,0) = iomega_n+mu-epsilon_kckl[kc+1][kl];
        
        G_kl(freq,kc  ,kc  ,1) = iomega_n+mu-epsilon_kckl[kc  ][kl];
        G_kl(freq,kc  ,kc+1,1) = h/2;
        G_kl(freq,kc+1,kc  ,1) = h/2;
        G_kl(freq,kc+1,kc+1,1) = iomega_n+mu-epsilon_kckl[kc+1][kl];
      }
    }
    G_kl = green_function_subtract(G_kl, selfenergy);
    G_kl = green_function_inverse(G_kl);
    G0_omega_new = green_function_add(G0_omega_new, G_kl);
  }
  std::cout<<"done running loop over cells."<<std::endl;
  
  G0_omega_new = green_function_divide(G0_omega_new, std::complex<double>(n_cluster_cells()));
  G0_omega_new = green_function_inverse(G0_omega_new);
  G0_omega_new = green_function_add(G0_omega_new, selfenergy);
  G0_omega_new = green_function_inverse(G0_omega_new);
  return G0_omega_new;
}




matsubara_green_function_t AntiFerroDCATransformer::initial_G0(const alps::Parameters& parms, shape_t shape)
{
  if(shape != blockdiagonal) throw std::logic_error("the shape really should be blockdiagonal.");
  if (parms.defined("G0OMEGA_INPUT")){
    std::ifstream in_file(parms["G0OMEGA_INPUT"].c_str());
    if(in_file.good()){
      return ClusterTransformer::initial_G0(parms, shape);
    }
    else
      std::cout<<"WARNING: Input file could not be found, proceeding with free solution."<<std::endl;
  }
  std::cout<<"calculating initial G0_omega"<<std::endl;
  unsigned int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  unsigned int n_tau=boost::lexical_cast<int>(parms["N"]);
  unsigned int n_flavor=parms.value_or_default("FLAVORS", 2);
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  double h = static_cast<double>(parms["H_INIT"]);
  matsubara_green_function_t G_omega(n_matsubara,n_cluster_sites(),n_flavor, blockdiagonal);
  matsubara_green_function_t G0_omega(n_matsubara,n_cluster_sites(),n_flavor, blockdiagonal);
  G0_omega.clear();
  for (spin_t flavor=0; flavor<n_flavor; ++flavor) {
    for (site_t i=0; i<n_cluster_sites(); ++i) {
      for (frequency_t freq=0; freq<n_matsubara; ++freq) {
        G0_omega(freq, i, i, flavor) = 1.;
      }
    }
  }
  G_omega = G0_omega;
  //invoke transformation with G_omega==G0_omega (selfenergy=0)
  matsubara_green_function_t G0_omega_new = operator()(G_omega, G0_omega, mu, h, beta);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  G0_omega_new.write("G0_omega_input");
  transform_into_real_space(G0_omega_new).write("G0_omega_input_realspace");
  itime_green_function_t G0_tau(n_tau, n_cluster_sites(), n_flavor, blockdiagonal);
  fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega_new);
  G0_tau.write("G0_tau_input");
  transform_into_real_space(G0_tau).write("G0_tau_input_realspace");
  return G0_omega_new;
}





std::map<double, ClusterHelper::Matrix<std::complex<double> >::Type>
AntiFerroDCATransformer::local_green_function_from_selfenergies(const std::map<double,
                                                                std::vector<ClusterHelper::Matrix<std::complex<double> >::Type> >& /*selfenergies*/,
                                                                const double /*mu*/, const double /*delta*/)
{
  boost::throw_exception(std::invalid_argument("to be implemented"));
  std::map<double, ClusterHelper::Matrix<std::complex<double> >::Type> dummy;
  return dummy;
}


alps::RealObsevaluator AntiFerroDCATransformer::kinetic_energy(const std::vector<std::vector<std::vector<alps::RealObsevaluator> > >& /*selfenergies*/,
                                                               const double /*beta*/, const double /*mu*/,
                                                               const bool /*error_propagation*/, const bool /*paramagnet*/)
{
  boost::throw_exception(std::invalid_argument("to be implemented"));
  return alps::RealObsevaluator();
}





///compute the spin symmetries of the Green's function
matsubara_green_function_t AntiFerroDCATransformer::spin_symmetrizer(const matsubara_green_function_t& G, const bool paramagnet) const
{
  ///DEBUG: this is implemented strictly for the half filled AFM case
  std::cerr<<"in overloaded symmetrizer. Assuming half filled mu=0"<<std::endl;
  //this spin symmetrizer does not do anything.
  if (paramagnet){
    boost::throw_exception(std::invalid_argument("AFM and PM? redundant parameters! set paramagnet to false.\n"));
  }
  matsubara_green_function_t G_new(G);
  for(int n=0;n<G.nfreq();++n){
    for(int k=0;k<G.nsite();k+=2){
      int q=k+1;
      //spin-symmetrize the diagonal components:
      std::complex<double> GK=(G(n,k,k,0)+G(n,k,k,1)-std::conj(G(n,q,q,0))-std::conj(G(n,q,q,1)))/4.;
      std::complex<double> GKQ=(G(n,k,q,0)-G(n,k,q,1)+G(n,q,k,0)-G(n,q,k,1))/4.;
      if(n==0){
        //std::cout<<" orig form: "<<G(n,k,k,0)<<" "<<G(n,k,k,1)<<" "<<G(n,q,q,0)<<" "<<G(n,q,q,1)<<" symmetrized to: "<<GK<<std::endl;
        std::cout<<" symmetrized form diag: "<<G(n,k,k,0)<<" "<<G(n,k,k,1)<<" "<<-std::conj(G(n,q,q,0))<<" "<<-std::conj(G(n,q,q,1))<<" symmetrized to: "<<GK<<std::endl;
        std::cout<<" symmetrized form offd: "<<G(n,k,q,0)<<" "<<-G(n,k,q,1)<<" "<<G(n,q,k,0)<<" "<<-G(n,q,k,1)<<" symmetrized to: "<<GKQ<<std::endl;
      }
      G_new(n,k,k,0)=GK;
      G_new(n,k,k,1)=GK;
      G_new(n,q,q,0)=-std::conj(GK);
      G_new(n,q,q,1)=-std::conj(GK);
      
      G_new(n,k,q,0)=GKQ;
      G_new(n,k,q,1)=-GKQ;
      G_new(n,q,k,0)=GKQ;
      G_new(n,q,k,1)=-GKQ;
      
    }
  }
  return G_new;
}

///GF operation: compute the inverse of a Matsubara GF
matsubara_green_function_t AntiFerroDCATransformer::green_function_inverse(const matsubara_green_function_t& G)
{
  matsubara_green_function_t G_inv(G);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (std::size_t k=0; k<G.nsite(); k+=2) {
        std::complex<double> detG = G(t, k, k, f)*G(t, k+1, k+1, f) - G(t, k+1, k, f)*G(t, k, k+1, f);
        if(std::abs(detG) < 1.e-15) throw std::runtime_error("there's a problem with your determinant in the inverse: it is near singular!");
        G_inv(t, k  , k  , f) =  G(t, k+1, k+1, f)/detG;
        G_inv(t, k  , k+1, f) = -G(t, k  , k+1, f)/detG;
        G_inv(t, k+1, k  , f) = -G(t, k+1, k  , f)/detG;
        G_inv(t, k+1, k+1, f) =  G(t, k  , k  , f)/detG;
      }
    }
  }
  return G_inv;
}

///GF operation: add a Matsubara GF to another Matsubara GF and return the result
matsubara_green_function_t AntiFerroDCATransformer::green_function_add(const matsubara_green_function_t& G1,
                                                                       const matsubara_green_function_t& G2)
{
  matsubara_green_function_t G(G1);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (unsigned int k=0; k<G.nsite(); k+=2) {
        G(t, k  , k  , f) += G2(t, k  , k  , f);
        G(t, k  , k+1, f) += G2(t, k  , k+1, f);
        G(t, k+1, k  , f) += G2(t, k+1, k  , f);
        G(t, k+1, k+1, f) += G2(t, k+1, k+1, f);
      }
    }
  }
  return G;
}


///GF operation: subtract a Matsubara GF from another Matsubara GF and return the result
matsubara_green_function_t AntiFerroDCATransformer::green_function_subtract(const matsubara_green_function_t& G1,
                                                                            const matsubara_green_function_t& G2)
{
  green_function<std::complex<double> > G(G1);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (std::size_t k=0; k<G.nsite(); k+=2) {
        G(t, k  , k  , f) -= G2(t, k  , k  , f);
        G(t, k  , k+1, f) -= G2(t, k  , k+1, f);
        G(t, k+1, k  , f) -= G2(t, k+1, k  , f);
        G(t, k+1, k+1, f) -= G2(t, k+1, k+1, f);
      }
    }
  }
  return G;
}


///GF operation: divide a Matsubara GF by a number
matsubara_green_function_t AntiFerroDCATransformer::green_function_divide(const matsubara_green_function_t& G1,
                                                                          const std::complex<double> d)
{
  green_function<std::complex<double> > G(G1);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (unsigned int k=0; k<G.nsite(); k+=2) {
        G(t, k, k, f) /= d;
        G(t, k+1, k+1, f) /= d;
        G(t, k+1, k, f) /= d;
        G(t, k, k+1, f) /= d;
      }
    }
  }
  return G;
}

