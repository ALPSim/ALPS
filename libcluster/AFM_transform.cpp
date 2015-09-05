/*****************************************************************************
 *
 * Copyright (C) 2007-2015 by Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#include "alps/numeric/matrix.hpp"

#include "cluster.h"
#include "AFMDCA_transform.hxx"
#include "fourier_transform.h"

matsubara_green_function_t AntiFerroDCATransformer::operator()(const matsubara_green_function_t &G_omega,
                                                            const matsubara_green_function_t &G0_omega,
                                                            double &mu, const double h, const double beta)
{
  matsubara_green_function_t selfenergy(G_omega);
  return this->operator()(G_omega, G0_omega, selfenergy, mu, h, beta, false);
}

///self-consistency transform operator for the Antiferromagnetic DCA case
matsubara_green_function_t AntiFerroDCATransformer::operator()(const matsubara_green_function_t &G_omega,
                                                            const matsubara_green_function_t &G0_omega,
                                                            matsubara_green_function_t &selfenergy,
                                                            double &mu, const double h, const double beta, const bool use_selfenergy)
{
  //std::cout<<"runnning operator() with AFM"<<std::endl;
  unsigned int nflavor=G0_omega.nflavor();
  unsigned int nfreq=G0_omega.nfreq();
  
  std::cout<<"using a magnetic field of: "<<h<<std::endl;
  
  //first: compute the self energy using G and G0.
  matsubara_green_function_t G0_omega_new(nfreq,n_cluster_sites(),nflavor, blockdiagonal);
  G0_omega_new.clear();
  if(!use_selfenergy){ //compute selfenergy out of G and G0 first
    matsubara_green_function_t G0_omega_inv = green_function_inverse (G0_omega);
    matsubara_green_function_t G_omega_inv  = green_function_inverse (G_omega );
    selfenergy   = green_function_subtract(G0_omega_inv, G_omega_inv);
  }
  //buffer the dispersion and the weight – no need to evaluate this for every frequency
  //first dimension: sign
  //second dimension: k0 or k1
  //third dimension: cluster momentum kc
  //fourth dimension: lattice momentum kl
  std::vector<std::vector<double> >epsilon_kckl(n_cluster_sites(), std::vector<double>(n_lattice_momenta()));
  for (int kc=0; kc<n_cluster_sites(); ++kc) {
    for (int kl=0; kl<n_lattice_momenta(); ++kl) {
      std::vector<double> k(dimension());
      for (unsigned int d=0; d<dimension(); ++d){
        k[d] = lattice_momentum(kl,d) + cluster_momentum(kc  ,d);
      }
      epsilon_kckl[kc][kl]=epsilon(k);
    }
  }
  std::vector<double> weight_kl(n_lattice_momenta());
  for (int kl=0; kl<n_lattice_momenta(); ++kl) {
    weight_kl[kl]= lattice_momentum(kl,dimension());
  }
  
  //this here implements an AFM two-site formalism and then rotates to k-space.
  //k-space rotation matrix
  alps::numeric::matrix<std::complex<double> > U(2,2); U(0,0)=1.; U(1,0)=1.; U(0,1)=1.; U(1,1)=-1.; U/=sqrt(2.);
  //some auxiliary matrices
  alps::numeric::matrix<std::complex<double> > Ginv_2site_up(2,2),Ginv_2site_dn(2,2),Ginv_kq_up(2,2),Ginv_kq_dn(2,2);
  matsubara_green_function_t G_kl(G0_omega);
  for (int kl=0; kl<n_lattice_momenta(); ++kl) {
    for (int freq=0; freq<nfreq; ++freq) {
      std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      for (int kc=0; kc<n_cluster_sites(); kc++) {
        std::vector<double> k(dimension());
        for (unsigned int d=0; d<dimension(); ++d){
          k[d] = lattice_momentum(kl,d) + cluster_momentum(kc  ,d);
        }
        
        Ginv_2site_up(0,0)=iomega_n+mu+h;
        Ginv_2site_up(0,1)=epsilon_kckl[kc][kl];
        Ginv_2site_up(1,0)=epsilon_kckl[kc][kl];
        Ginv_2site_up(1,1)=iomega_n+mu-h;
        
        Ginv_2site_dn(0,0)=iomega_n+mu-h;
        Ginv_2site_dn(0,1)=epsilon_kckl[kc][kl];
        Ginv_2site_dn(1,0)=epsilon_kckl[kc][kl];
        Ginv_2site_dn(1,1)=iomega_n+mu+h;
        
        Ginv_kq_up=U*Ginv_2site_up*U;
        Ginv_kq_dn=U*Ginv_2site_dn*U;
        
        //std::cout<<U<<" "<<Ginv_2site_up<<" "<<Ginv_kq_up<<std::endl;
        
        //std::cout<<kl<<" "<<kc<<" "<<freq<<Ginv_kq_up<<std::endl;
        
        int K=(kc%2==0)?kc:kc-1; int Q=(kc%2==0)?kc+1:kc;
        //std::cout<<"K: "<<K<<" "<<Q<<" "<<std::endl;
        G_kl(freq,K,K,0) = (Ginv_kq_up(0,0)-selfenergy(freq,K,K,0))/weight_kl[kl];
        G_kl(freq,K,Q,0) = (Ginv_kq_up(0,1)-selfenergy(freq,K,Q,0))/weight_kl[kl];
        G_kl(freq,Q,K,0) = (Ginv_kq_up(1,0)-selfenergy(freq,Q,K,0))/weight_kl[kl];
        G_kl(freq,Q,Q,0) = (Ginv_kq_up(1,1)-selfenergy(freq,Q,Q,0))/weight_kl[kl];
        
        G_kl(freq,K,K,1) = (Ginv_kq_dn(0,0)-selfenergy(freq,K,K,1))/weight_kl[kl];
        G_kl(freq,K,Q,1) = (Ginv_kq_dn(0,1)-selfenergy(freq,K,Q,1))/weight_kl[kl];
        G_kl(freq,Q,K,1) = (Ginv_kq_dn(1,0)-selfenergy(freq,Q,K,1))/weight_kl[kl];
        G_kl(freq,Q,Q,1) = (Ginv_kq_dn(1,1)-selfenergy(freq,Q,Q,1))/weight_kl[kl];
        //std::cout<<"end "<<kl<<" "<<freq<<" "<<kc<<std::endl;
      }
    }
    //G_kl = green_function_subtract(G_kl, selfenergy);
    G_kl = green_function_inverse(G_kl);
    G0_omega_new = green_function_add(G0_omega_new, G_kl);
  }
  
  G0_omega_new = green_function_divide(G0_omega_new, std::complex<double>(n_cluster_cells()));
  G0_omega_new = green_function_inverse(G0_omega_new);
  G0_omega_new = green_function_add(G0_omega_new, selfenergy);
  G0_omega_new = green_function_inverse(G0_omega_new);
  return G0_omega_new;
}

/*matsubara_green_function_t AntiFerroDCATransformer::operator()(const matsubara_green_function_t &G_omega,
 const matsubara_green_function_t &G0_omega,
 double &mu, const double h, const double beta)
 {
 unsigned int nflavor=G0_omega.nflavor();
 unsigned int nfreq=G0_omega.nfreq();
 assert(G0_omega.nsite()==n_cluster_sites());
 assert(nflavor>1 && nflavor%2 == 0);
 matsubara_green_function_t G0_omega_new(nfreq,n_cluster_sites(),nflavor, blockdiagonal);
 G0_omega_new.clear();
 matsubara_green_function_t G0_omega_inv = green_function_inverse(G0_omega);
 matsubara_green_function_t G_omega_inv = green_function_inverse(G_omega);
 matsubara_green_function_t selfenergy = green_function_subtract(G0_omega_inv, G_omega_inv);
 for (int kl=0; kl<n_cluster_cells(); ++kl) {
 for (int sign=-1; sign<=1; sign+=2) {
 matsubara_green_function_t G_kl(G0_omega);
 for (int freq=0; freq<nfreq; ++freq) {
 std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
 for (int kc=0; kc<n_cluster_sites(); kc+=2) {
 std::vector<double> k0(dimension());
 for (unsigned int d=0; d<dimension(); ++d)
 k0[d] = sign*lattice_momentum(kl,d) + cluster_momentum(kc,d);
 std::vector<double> k1(dimension());
 for (unsigned int d=0; d<dimension(); ++d)
 k1[d] = sign*lattice_momentum(kl,d) + cluster_momentum(kc+1,d);
 for (int flavor=0; flavor<nflavor; flavor+=2) {
 G_kl(freq,kc,kc,flavor) = iomega_n+mu-epsilon(k0);
 G_kl(freq,kc,kc,flavor+1) = iomega_n+mu-epsilon(k0);
 G_kl(freq,kc+1,kc+1,flavor) = iomega_n+mu-epsilon(k1);
 G_kl(freq,kc+1,kc+1,flavor+1) = iomega_n+mu-epsilon(k1);
 G_kl(freq,kc+1,kc,flavor) = -h/2;
 G_kl(freq,kc+1,kc,flavor+1) = h/2;
 G_kl(freq,kc,kc+1,flavor) = -h/2;
 G_kl(freq,kc,kc+1,flavor+1) = h/2;
 }
 }
 }
 G_kl = green_function_subtract(G_kl, selfenergy);
 G_kl = green_function_inverse(G_kl);
 G0_omega_new = green_function_add(G0_omega_new, G_kl);
 }
 }
 G0_omega_new = green_function_divide(G0_omega_new, std::complex<double>(2*n_cluster_cells()));
 G0_omega_new = green_function_inverse(G0_omega_new);
 G0_omega_new = green_function_add(G0_omega_new, selfenergy);
 G0_omega_new = green_function_inverse(G0_omega_new);
 return G0_omega_new;
 }*/



matsubara_green_function_t AntiFerroDCATransformer::initial_G0(const alps::params& parms, shape_t shape)
{
  if(shape != blockdiagonal) throw std::logic_error("the shape really should be blockdiagonal.");
  if (parms.defined("G0OMEGA_INPUT")){
    std::string g0path=parms["G0OMEGA_INPUT"];
    std::ifstream in_file(g0path.c_str());
    if(in_file.good()){
      return ClusterTransformer::initial_G0(parms, shape);
    }
    else
      std::cout<<"WARNING: Input file could not be found, proceeding with free solution."<<std::endl;
  }
  std::cout<<"calculating initial G0_omega"<<std::endl;
  unsigned int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  unsigned int n_tau=boost::lexical_cast<int>(parms["N"]);
  unsigned int n_flavor=parms["FLAVORS"]|2;
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  double h = static_cast<double>(parms["H_INIT"]);
  matsubara_green_function_t G_omega(n_matsubara,n_cluster_sites(),n_flavor, blockdiagonal);
  matsubara_green_function_t G0_omega(n_matsubara,n_cluster_sites(),n_flavor, blockdiagonal);
  G0_omega.clear();
  for (int flavor=0; flavor<n_flavor; ++flavor) {
    for (int i=0; i<n_cluster_sites(); ++i) {
      for (int freq=0; freq<n_matsubara; ++freq) {
        G0_omega(freq, i, i, flavor) = 1.;
      }
    }
  }
  G_omega = G0_omega;
  //invoke transformation with G_omega==G0_omega (selfenergy=0)
  matsubara_green_function_t G0_omega_new = operator()(G_omega, G0_omega, mu, h, beta);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_cluster_transformer(parms, fourier_ptr);
  G0_omega_new.write("G0_omega_input");
  std::ofstream G0_omega_re_file("G0_omega_input_re.dat");
  std::ofstream G0_omega_im_file("G0_omega_input_im.dat");
  print_imag_green_matsubara(G0_omega_im_file, G0_omega_new, beta);
  print_real_green_matsubara(G0_omega_re_file,G0_omega_new,beta);
  
  transform_into_real_space(G0_omega_new).write("G0_omega_input_realspace");
  itime_green_function_t G0_tau(n_tau, n_cluster_sites(), n_flavor, blockdiagonal);
  fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega_new);
  G0_tau.write("G0_tau_input");
  transform_into_real_space(G0_tau).write("G0_tau_input_realspace");
  return G0_omega_new;
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
  return G_new;
  ///DEBUG
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
    for (int f=0; f<G.nflavor(); ++f) {
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
    for (int f=0; f<G.nflavor(); ++f) {
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
    for (int f=0; f<G.nflavor(); ++f) {
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
    for (int f=0; f<G.nflavor(); ++f) {
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

