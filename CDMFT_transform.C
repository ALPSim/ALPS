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


#include "CDMFT_transform.h"
#include <boost/numeric/ublas/io.hpp>
#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/lapack/driver/gesv.hpp>
#include "CDMFT_transformer.hxx"
#include "fouriertransform.h"

typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;

/*
 *
 *
 *
 These are all foursite CDMFT things. relatively easily extendable to other four-site things...
 *
 *
 *
 */

//Auxiliary functions
//invert a matrix
void dense_matrix_inverse(complex_dense_matrix &M /*input*/, complex_dense_matrix &M_inv/*output*/)
{
  M_inv = boost::numeric::ublas::identity_matrix<std::complex<double> >(M.size1());
  boost::numeric::ublas::vector<fortran_int_t> ipivot(M.size1());
  //memset(&(ipivot[0]), 0, sizeof(fortran_int_t)*M.size1());
  boost::numeric::bindings::lapack::gesv(M,ipivot, M_inv);
}

//write a gf to a matrix
template<class T> void gf_to_matrix(int omega, int flavor, const green_function<T> &gf, boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> &M){
  for(unsigned int i=0;i<gf.nsite();++i){
    for(unsigned int j=0;j<gf.nsite();++j){
      M(i,j)=gf(omega, i,j,flavor);
    }
  }
}
//write a matrix to a gf
template<class T> void matrix_to_gf(int omega, int flavor,  green_function<T> &gf, const boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> &M){
  for(unsigned int i=0;i<gf.nsite();++i){
    for(unsigned int j=0;j<gf.nsite();++j){
      gf(omega, i,j,flavor)=M(i,j);
    }
  }
}

CDMFTTransformer::CDMFTTransformer(const alps::params& params)
: ClusterTransformer(params),
t_(params.defined("W") ? 1. : params["t"]),
t2_(params["t'"]|0.),
t3_(params["t3"]|0.),
t4_(params["t4"]|0.),
bandwidth_(params["W"]|0.),
bandwidth_scale_factor_(1.),
L_(params["L"]),
half_filling_(0.),
dispersion_moments_(2, 0.)
{
  if (params.defined("t") && params.defined("W"))
    boost::throw_exception(std::invalid_argument("both t and W are defined"));
  if (params.defined("W"))
    boost::throw_exception(std::invalid_argument("W is defined, check that it works."));

  //calculating dispersion moments, bandwidth and half filling
  std::cout<<"warning: the CDMFT transformer has NO bandwidths, dispersions, etc implemented.";
}
matsubara_green_function_t FoursiteCDMFTTransformer::operator()(const matsubara_green_function_t &G_omega,
                                                                const matsubara_green_function_t &G0_omega,
                                                                double &mu, const double h, const double beta)
{
  if(h!=0) { // && G0_omega.nsite()!=1 <-this not understood - looks like a bug to me.
    std::cerr<<"H is: "<<h<<std::endl;
    throw std::invalid_argument("don't know how to handle magn. fields (use the AFM transformer for that!), exiting");
  }
  unsigned int nflavor=G0_omega.nflavor();
  unsigned int nfreq=G0_omega.nfreq();
  if(G0_omega.nsite()!=n_cluster_sites()) throw std::runtime_error("problem with cluster sites");
    matsubara_green_function_t selfenergy(G0_omega);
    matsubara_green_function_t G0_omega_new(G0_omega);
    matsubara_green_function_t G_omega_lattice(G0_omega);
    
    std::cout<<"computing self energy."<<std::endl;
  //compute the Matsubara self energy
  complex_dense_matrix G0_mat(n_cluster_sites(),n_cluster_sites()), G_mat(n_cluster_sites(),n_cluster_sites()), Sigma_mat(n_cluster_sites(),n_cluster_sites());
  complex_dense_matrix G0_inv_mat(n_cluster_sites(),n_cluster_sites()), G_inv_mat(n_cluster_sites(),n_cluster_sites());
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (frequency_t omega=0; omega<nfreq; ++omega) {
      gf_to_matrix(omega, flavor, G0_omega, G0_mat);
      gf_to_matrix(omega, flavor, G_omega, G_mat);
      dense_matrix_inverse(G0_mat, G0_inv_mat);
      dense_matrix_inverse(G_mat, G_inv_mat);
      Sigma_mat=G0_inv_mat-G_inv_mat;
      matrix_to_gf(omega, flavor, selfenergy, Sigma_mat);
    }
  }
  std::cout<<"computing k summation."<<std::endl;
  
  selfenergy = transform_into_real_space(selfenergy); //the integration is done in real space
  
  //the k-summation
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      complex_dense_matrix diag=boost::numeric::ublas::identity_matrix<std::complex<double> >(n_cluster_sites())*(iomega_n+mu);
      complex_dense_matrix H0(n_cluster_sites(), n_cluster_sites());
      //complex_dense_matrix H0K(n_cluster_sites(), n_cluster_sites());
      complex_dense_matrix tmp(n_cluster_sites(), n_cluster_sites());
      complex_dense_matrix Ginv(n_cluster_sites(), n_cluster_sites());
      G_mat=boost::numeric::ublas::zero_matrix<std::complex<double> >(n_cluster_sites(),n_cluster_sites());
      H0   =boost::numeric::ublas::zero_matrix<std::complex<double> >(n_cluster_sites(),n_cluster_sites());
      tmp  =boost::numeric::ublas::zero_matrix<std::complex<double> >(n_cluster_sites(),n_cluster_sites());
      gf_to_matrix(freq, flavor, selfenergy, Sigma_mat);
      for(int nx=0;nx<L_;++nx){
        double kx=-M_PI+2*M_PI*nx/L_;
        for(int ny=0;ny<L_;++ny){
          //transformation of the H0 - matrix to k-space
          double ky=-M_PI+2*M_PI*ny/L_;
          set_Tij_matrix(kx,ky,H0);
          //Ginv = iomega +mu - epsilon - Sigma
          Ginv=diag+H0+Sigma_mat;
          dense_matrix_inverse(Ginv,tmp);
          G_mat+=tmp/(L_*L_);
        }
      }
      matrix_to_gf(freq, flavor, G_omega_lattice, G_mat);
    }
  }
  std::cout<<"lattice GF is: "<<std::endl;
  std::ofstream lgf_file("lgf.dat");
  std::cout<<"in real space: "<<std::endl;
  lgf_file<<std::make_pair(G_omega_lattice,beta)<<std::endl;
  std::ofstream lgfks_file("lgfks.dat");
  std::cout<<"in k space: "<<std::endl;
  G_omega_lattice = transform_into_k_space(G_omega_lattice);
  lgfks_file<<std::make_pair(G_omega_lattice,beta)<<std::endl;
  
  std::cout<<"computing G0 from sigma"<<std::endl;
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (frequency_t omega=0; omega<nfreq; ++omega) {
      gf_to_matrix(omega, flavor, selfenergy, Sigma_mat);
      gf_to_matrix(omega, flavor, G_omega_lattice, G_mat);
      dense_matrix_inverse(G_mat, G_inv_mat);
      G0_inv_mat=Sigma_mat+G_inv_mat;
      dense_matrix_inverse(G0_inv_mat, G0_mat);
      matrix_to_gf(omega, flavor, G0_omega_new, G0_mat);
    }
  }
  std::cout<<"returning from CDMFT self consistency"<<std::endl;
  std::ofstream g0newks_file("g0newks.dat");
  g0newks_file<<std::make_pair(G0_omega_new,beta)<<std::endl;
  return G0_omega_new;
}



matsubara_green_function_t FoursiteCDMFTTransformer::initial_G0(const alps::params& parms, shape_t shape)
{
  if(shape != nondiagonal) throw std::logic_error("shape should probably be nondiagonal");
  if (parms.defined("G0OMEGA_INPUT")){
    std::string g0path=parms["G0OMEGA_INPUT"];
    std::ifstream in_file(g0path.c_str());
    if(in_file.good()){
      return ClusterTransformer::initial_G0(parms, shape);
    }
    else
      std::cout<<"WARNING: Input file could not be found, proceeding with free solution."<<std::endl;
  }
  //proceed with computation of new G0 if we can't find input G0
  std::cout<<"calculating initial G0_omega"<<std::endl;
  unsigned int n_matsubara=parms["NMATSUBARA"];
  unsigned int n_tau=parms["N"];
  unsigned int n_flavor=parms["FLAVORS"]|2;
  double beta = parms["BETA"];
  double mu =parms["MU"];
  matsubara_green_function_t G0_omega_new(n_matsubara,n_cluster_sites(),n_flavor, nondiagonal);
  G0_omega_new.clear();
  if (parms["INSULATING"]|false) {
    std::cout <<"calculating insulating G0_omega";
    throw std::runtime_error("insulating g0 in foursite transformer not implemented.");
  }
  else if (parms["ISOLATED"]|false) {
    std::cout << "calculating G0_omega for isolated cluster"<<std::endl;
    throw std::runtime_error("isolated g0 in foursite transformer not implemented.");
  }
  else {
    matsubara_green_function_t G0_omega(n_matsubara,n_cluster_sites(),n_flavor, nondiagonal);
    matsubara_green_function_t G_omega(n_matsubara,n_cluster_sites(),n_flavor, nondiagonal);
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
    G0_omega_new = operator()(G_omega, G0_omega, mu, 0, beta);
  }
  itime_green_function_t G0_tau(n_tau+1, n_cluster_sites(), n_flavor, nondiagonal);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  G0_omega_new.write("G0_omega_input");
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega_new);
  G0_tau.write("G0_tau_input");
  matsubara_green_function_t G0_omega_real = transform_into_real_space(G0_omega_new);
  G0_omega_real.write("G0_omega_input_realspace");
  itime_green_function_t G0_tau_real = transform_into_real_space(G0_tau);
  G0_tau_real.write("G0_tau_input_realspace");
  std::ofstream g0_omega_input_hr("g0_omega_input_hr.dat"); g0_omega_input_hr<<std::make_pair(G0_omega_new,beta);
  std::ofstream g0_omega_input_rs_hr("g0_omega_input_rs_hr.dat"); g0_omega_input_rs_hr<<std::make_pair(G0_omega_real,beta);
  return G0_omega_new;
}


nambu_matsubara_green_function_t FoursiteCDMFTTransformer::initial_G0_nambu(alps::params const&){
  throw std::runtime_error("initial nambu_matsubara_green_function_t CDMFT transformer not implemented.");
}
nambu_matsubara_green_function_t FoursiteCDMFTTransformer::operator()(const nambu_matsubara_green_function_t &G_omega,
                                                                      const nambu_matsubara_green_function_t &G0_omega,
                                                                      nambu_matsubara_green_function_t &selfenergy,
                                                                      double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy){
  throw std::runtime_error("nambu operator() in CDMFT transformer not implemented.");
}
/*nambu_green_function<double> FoursiteCDMFTTransformer::transform_into_k_space(nambu_green_function<double> const& gf) const{
 throw std::runtime_error("transform_into_k_space CDMFT transformer not implemented.");
 }
 nambu_green_function<std::complex<double> > FoursiteCDMFTTransformer::transform_into_k_space(nambu_green_function<std::complex<double> > const& gf) const{
 throw std::runtime_error("transform_into_k_space CDMFT transformer not implemented.");
 }*/
matsubara_green_function_t FoursiteCDMFTTransformer::transform_into_real_space(const matsubara_green_function_t& G) const{
  spin_t nflavor=G.nflavor();
  frequency_t nfreq=G.nfreq();
  matsubara_green_function_t G_transformed(G);
  complex_dense_matrix U=U_;
  complex_dense_matrix tmp=U_;
  complex_dense_matrix G_mat(n_cluster_sites(),n_cluster_sites()), G_mat_new(n_cluster_sites(),n_cluster_sites());
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (frequency_t omega=0; omega<nfreq; ++omega) {
      gf_to_matrix(omega, flavor, G, G_mat);
      tmp=boost::numeric::ublas::prod(G_mat,U);
      G_mat_new=boost::numeric::ublas::prod(U,tmp);
      matrix_to_gf(omega, flavor, G_transformed, G_mat_new);
    }
  }
  return G_transformed;
}
matsubara_green_function_t FoursiteCDMFTTransformer::transform_into_k_space(const matsubara_green_function_t& G) const{
  return transform_into_real_space(G); //this is a special foursite case.
}
itime_green_function_t FoursiteCDMFTTransformer::transform_into_real_space(const itime_green_function_t& G) const{
  spin_t nflavor=G.nflavor();
  itime_index_t ntime=G.ntime();
  itime_green_function_t G_transformed(G);
  dense_matrix U(U_.size1(), U_.size1());
  for(unsigned int i=0;i<U_.size1();++i){ for (unsigned int j=0;j<U_.size1();++j){ U(i,j)=U_(i,j).real();} }
  dense_matrix tmp=U;
  dense_matrix G_mat(n_cluster_sites(),n_cluster_sites()), G_mat_new(n_cluster_sites(),n_cluster_sites());
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (itime_index_t tau=0; tau<ntime; ++tau) {
      gf_to_matrix(tau, flavor, G, G_mat);
      tmp=boost::numeric::ublas::prod(G_mat,U);
      G_mat_new=boost::numeric::ublas::prod(U,tmp);
      matrix_to_gf(tau, flavor, G_transformed, G_mat_new);
    }
  }
  return G_transformed;
}
itime_green_function_t FoursiteCDMFTTransformer::transform_into_k_space(const itime_green_function_t& G) const{
  return transform_into_real_space(G);
}
nambu_matsubara_green_function_t FoursiteCDMFTTransformer::transform_into_real_space(const nambu_matsubara_green_function_t& G) const{
  throw std::runtime_error("this has not yet been implemented.");
}
nambu_matsubara_green_function_t FoursiteCDMFTTransformer::transform_into_k_space(const nambu_matsubara_green_function_t& G) const{
  throw std::runtime_error("this has not yet been implemented.");
}
nambu_itime_green_function_t FoursiteCDMFTTransformer::transform_into_real_space(const nambu_itime_green_function_t& G) const{
  throw std::runtime_error("this has not yet been implemented.");
}
nambu_itime_green_function_t FoursiteCDMFTTransformer::transform_into_k_space(const nambu_itime_green_function_t& G) const{
  throw std::runtime_error("this has not yet been implemented.");
}
void FoursiteCDMFTTransformer::set_Tij_matrix(double kx, double ky, complex_dense_matrix &H0){
  H0(0,0)=0;  H0(1,1)=0; H0(2,2)=0; H0(3,3)=0;
  double coskx=cos(kx);
  double cosky=cos(ky);
  
  std::complex<double> exp_ikx=std::exp(std::complex<double>(0, kx));
  std::complex<double> exp_mikx=std::conj(exp_ikx);
  std::complex<double> exp_iky=std::exp(std::complex<double>(0, ky));
  std::complex<double> exp_miky=std::conj(exp_iky);
  //these are the factors of 't'
  H0(0,1)=-2.*t_*coskx*exp_ikx;
  H0(0,3)=-2.*t_*cosky*exp_miky;
  H0(2,1)=-2.*t_*cosky*exp_iky;
  H0(2,3)=-2.*t_*coskx*exp_mikx;
  H0(1,0)=-2.*t_*coskx*exp_mikx;
  H0(1,2)=-2.*t_*cosky*exp_miky;
  H0(3,0)=-2.*t_*cosky*exp_iky;
  H0(3,2)=-2.*t_*coskx*exp_ikx;
  //these are the factors of 'tprime'
  throw std::logic_error("go to xi and ask for help.");
  /*H0(0,2)=-4.*std::exp(std::complex<double>(0, ( kx - ky)))*tprime_*coskx*cosky;
   H0(1,3)=-4.*std::exp(std::complex<double>(0, (-kx - ky)))*tprime_*coskx*cosky;
   H0(0,2)=-4.*std::exp(std::complex<double>(0, (-kx + ky)))*tprime_*coskx*cosky;
   H0(0,2)=-4.*std::exp(std::complex<double>(0, ( kx + ky)))*tprime_*coskx*cosky;*/
}

void FoursiteCDMFTTransformer::set_U_(unsigned int nsite, complex_dense_matrix &U_){
  
}

std::vector<double> FoursiteCDMFTTransformer::epsav() const{
  return epsav_;
}
std::vector<double> FoursiteCDMFTTransformer::epssqav() const{
  return epssqav_;
}
