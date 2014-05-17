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
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
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

//write a gf to a matrix normal state
template<class T> void gf_to_matrix(int omega, int flavor, const green_function<T> &gf, boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> &M){
  for(unsigned int i=0;i<gf.nsite();++i){
    for(unsigned int j=0;j<gf.nsite();++j){
      M(i,j)=gf(omega, i,j,flavor);
    }
  }
}
bool InvertMatrix(complex_dense_matrix &input, complex_dense_matrix &inverse)
{
  using namespace boost::numeric::ublas;
	typedef permutation_matrix<std::size_t> pmatrix;
  
	// create a working copy of the input
	complex_dense_matrix A(input);
  
	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());
  
	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;
  
	// create identity matrix of "inverse"
	inverse.assign(identity_matrix< std::complex<double> > (A.size1()));
  
	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);
  
	return true;
}
//write a gf to a matrix Nambu
template<class T> void gf_to_matrix(int omega, const nambu_green_function<T> &gf, boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> &M){
  for(unsigned int i=0;i<gf.nsite();++i){
    for(unsigned int j=0;j<gf.nsite();++j){
      M(2*i  ,2*j  )=gf(omega,i,j,0,0);
      M(2*i  ,2*j+1)=gf(omega,i,j,0,1);
      M(2*i+1,2*j  )=gf(omega,i,j,1,0);
      M(2*i+1,2*j+1)=gf(omega,i,j,1,1);
    }
  }
}
//write a matrix to a gf normal state
template<class T> void matrix_to_gf(int omega, int flavor,  green_function<T> &gf, const boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> &M){
  for(unsigned int i=0;i<gf.nsite();++i){
    for(unsigned int j=0;j<gf.nsite();++j){
      gf(omega, i,j,flavor)=M(i,j);
    }
  }
}
//write a matrix to a gf Nambu
template<class T> void matrix_to_gf(int omega, nambu_green_function<T> &gf, boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> &M){
  for(unsigned int i=0;i<gf.nsite();++i){
    for(unsigned int j=0;j<gf.nsite();++j){
      gf(omega,i,j,0,0)=M(2*i  ,2*j  );
      gf(omega,i,j,0,1)=M(2*i  ,2*j+1);
      gf(omega,i,j,1,0)=M(2*i+1,2*j  );
      gf(omega,i,j,1,1)=M(2*i+1,2*j+1);
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


matsubara_green_function_t LadderCDMFTTransformer::operator()(const matsubara_green_function_t &G_omega,
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



matsubara_green_function_t LadderCDMFTTransformer::initial_G0(const alps::params& params, shape_t shape)
{
  if(shape != nondiagonal) throw std::logic_error("shape should probably be nondiagonal");
  if (params.defined("G0OMEGA_INPUT")){
    std::string g0path=params["G0OMEGA_INPUT"];
    std::ifstream in_file(g0path.c_str());
    if(in_file.good()){
      return ClusterTransformer::initial_G0(params, shape);
    }
    else
      std::cout<<"WARNING: Input file could not be found, proceeding with free solution."<<std::endl;
  }
  //proceed with computation of new G0 if we can't find input G0
  std::cout<<"calculating initial G0_omega"<<std::endl;
  unsigned int n_matsubara=params["NMATSUBARA"];
  unsigned int n_tau=params["N"];
  unsigned int n_flavor=params["FLAVORS"]|2;
  double beta = params["BETA"];
  double mu =params["MU"];
  matsubara_green_function_t G0_omega_new(n_matsubara,n_cluster_sites(),n_flavor, nondiagonal);
  G0_omega_new.clear();
  if (params["INSULATING"]|false) {
    std::cout <<"calculating insulating G0_omega";
    throw std::runtime_error("insulating g0 in foursite transformer not implemented.");
  }
  else if (params["ISOLATED"]|false) {
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
  FourierTransformer::generate_transformer(params, fourier_ptr);
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


nambu_matsubara_green_function_t LadderCDMFTTransformer::initial_G0_nambu(const alps::params& params){
  if (params.defined("G0_OMEGA_INPUT")){
    std::string g0path=params["G0OMEGA_INPUT"];
    std::ifstream in_file(g0path.c_str());
    if(in_file.good()){
      return ClusterTransformer::initial_G0_nambu(params);
    }
    else
      std::cout<<"WARNING: Input file "<<params["G0OMEGA_INPUT"]<<" defined but could not be found, proceeding with free solution."<<std::endl;
  }
  //proceed with computation of new G0 if we can't find input G0
  std::cout<<"calculating initial G0_omega"<<std::endl;
  unsigned int n_matsubara=params["NMATSUBARA"];
  unsigned int n_tau=params["N"];
  unsigned int n_flavor=params["FLAVORS"]|2;
  double beta = params["BETA"];
  double mu =params["MU"];
  double eta1_init = params["ETA1_INIT"];
  double eta2_init = params["ETA2_INIT"];
  double h_init = params["H_INIT"];
  
  nambu_matsubara_green_function_t G0_omega_new(n_matsubara,n_cluster_sites(),n_flavor, nondiagonal);
  nambu_matsubara_green_function_t selfenergy(n_matsubara,n_cluster_sites(),n_flavor, nondiagonal);
  G0_omega_new.clear();
  if (params.defined("INSULATING")) {
    std::cout <<"calculating insulating G0_omega starting solution";
    for (frequency_t freq=0; freq<n_matsubara; ++freq) {
      const std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      const std::complex<double> G0 = 1./(iomega_n+mu);//-dispersion_moments(0));// for four site dispersion_moments(0)=0;
      for (spin_t flavor=0; flavor<n_flavor; ++flavor) {
        for (site_t s=0; s<n_cluster_sites(); ++s) {
          G0_omega_new(freq, s, s, 0,0) = G0;
          G0_omega_new(freq, s, s, 1,1) = G0;
          G0_omega_new(freq, s, s, 0,1) = 0.;
          G0_omega_new(freq, s, s, 1,0) = 0.;
        }
      }
    }
    
  }else if (params["ISOLATED"]|false){
    throw std::runtime_error("Isolated G0 not implemented in CDMFT");
  }
  else {
    nambu_matsubara_green_function_t G0_omega(n_matsubara,n_cluster_sites(),n_flavor, nondiagonal);
    nambu_matsubara_green_function_t G_omega(n_matsubara,n_cluster_sites(),n_flavor, nondiagonal);
    G0_omega.clear();
    for (site_t i=0; i<n_cluster_sites(); ++i) {
      for (frequency_t freq=0; freq<n_matsubara; ++freq) {
        G0_omega(freq, i, i, 0,0) = 1.;
        G0_omega(freq, i, i, 0,1) = 0.;
        G0_omega(freq, i, i, 1,0) = 0.;
        G0_omega(freq, i, i, 1,1) = 1.;
      }
    }
    //above a G0 function in k-space is initialized, so it needs to be transformed to real-space for CDMFT
    G0_omega_new=transform_into_real_space(G0_omega_new);
    
    G_omega = G0_omega;
    //invoke transformation with G_omega==G0_omega (selfenergy=0)
    std::cout<<"invoking free transform in initial."<<std::endl;
    G0_omega_new = operator()(G_omega, G0_omega, selfenergy, mu, h_init, eta1_init, eta2_init, beta, false);
  }
  
  nambu_itime_green_function_t G0_tau(n_tau+1, n_cluster_sites(), n_flavor, nondiagonal);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  G0_omega_new.write("G0_omega_input");
  FourierTransformer::generate_transformer(params, fourier_ptr);
  fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega_new);
  std::cout<<"write G0_tau_input"<<std::endl;
  G0_tau.write("G0_tau_input");

  return G0_omega_new;
  //throw std::runtime_error("initial nambu_matsubara_green_function_t CDMFT transformer not implemented.");
}
nambu_matsubara_green_function_t LadderCDMFTTransformer::operator()(const nambu_matsubara_green_function_t &G_omega,
                                                                      const nambu_matsubara_green_function_t &G0_omega,
                                                                      nambu_matsubara_green_function_t &selfenergy,
                                                                      double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy){
  std::cout<<"clustertrans starts"<<std::endl;

  unsigned int nflavor=G0_omega.nflavor();
  unsigned int nfreq=G0_omega.nfreq();
  unsigned int nsite=G0_omega.nsite();
  // if(G0_omega.nsite()!=n_cluster_sites()) throw{ std::runtime_error("problem with cluster sites");}
  // if(G0_omega.nflavor()!=2) throw {std::runtime_error("Nambu CDMFT can only handle flavor equals 2");}
  //nambu_matsubara_green_function_t selfenergy(G0_omega);
  nambu_matsubara_green_function_t G0_omega_new(G0_omega);
  nambu_matsubara_green_function_t G_omega_lattice(G0_omega);
  
  std::cout<<"computing self energy."<<std::endl;
  //compute the Matsubara self energy
  complex_dense_matrix G0_mat(n_cluster_sites()*2,n_cluster_sites()*2), G_mat(n_cluster_sites()*2,n_cluster_sites()*2), Sigma_mat(n_cluster_sites()*2,n_cluster_sites()*2);
  complex_dense_matrix G0_inv_mat(n_cluster_sites()*2,n_cluster_sites()*2), G_inv_mat(n_cluster_sites()*2,n_cluster_sites()*2);
  complex_dense_matrix G_mat2(n_cluster_sites()*2,n_cluster_sites()*2),G_inv_mat2(n_cluster_sites()*2,n_cluster_sites()*2);
  for (frequency_t omega=0; omega<nfreq; ++omega) {
    gf_to_matrix(omega, G0_omega, G0_mat);
    gf_to_matrix(omega, G_omega, G_mat);
    InvertMatrix(G0_mat, G0_inv_mat);
    InvertMatrix(G_mat, G_inv_mat);
    Sigma_mat=G0_inv_mat-G_inv_mat;
    
    matrix_to_gf(omega,selfenergy, Sigma_mat);
  }
  std::cout<<"computing k summation."<<std::endl;
  
  //selfenergy = transform_into_real_space(selfenergy); //the integration is done in real space
  
  //the k-summation
  
  for (frequency_t freq=0; freq<nfreq; ++freq) {
    std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
    complex_dense_matrix diag=boost::numeric::ublas::identity_matrix<std::complex<double> >(n_cluster_sites()*2)*(iomega_n);
    complex_dense_matrix H0(n_cluster_sites()*2, n_cluster_sites()*2);
    //complex_dense_matrix H0K(n_cluster_sites(), n_cluster_sites());
    complex_dense_matrix tmp(n_cluster_sites()*2, n_cluster_sites()*2);
    complex_dense_matrix Ginv(n_cluster_sites()*2, n_cluster_sites()*2);
    G_mat=boost::numeric::ublas::zero_matrix<std::complex<double> >(n_cluster_sites()*2,n_cluster_sites()*2);
    H0   =boost::numeric::ublas::zero_matrix<std::complex<double> >(n_cluster_sites()*2,n_cluster_sites()*2);
    tmp  =boost::numeric::ublas::zero_matrix<std::complex<double> >(n_cluster_sites()*2,n_cluster_sites()*2);
    gf_to_matrix(freq, selfenergy, Sigma_mat);
    
    for(site_t kl=0;kl<n_lattice_momenta();kl++){
      double kx=lattice_momentum(kl,0);
      double ky=lattice_momentum(kl,1);
      set_nambu_Tij_dwave_matrix(kx,ky,eta1,eta2,mu,h,H0);
    /*  if(kl==10&&freq==20){
        for(int i=0;i<2*n_cluster_sites();i++){
          for(int j=0;j<2*n_cluster_sites();j++){
            std::cout<<H0(i,j)<<"\t";
          }
          std::cout<<std::endl;
        }
      }*/
      Ginv=diag+H0-Sigma_mat;
      
      InvertMatrix(Ginv,tmp);
      
      G_mat+=tmp/(L_*L_)*lattice_momentum(kl,dimension());
      //G_mat+=tmp/(L_*L_);
    }
    matrix_to_gf(freq, G_omega_lattice, G_mat);
  }
  
  std::cout<<"lattice GF is: "<<std::endl;
  std::ofstream lgf_file("lgf.dat");
  std::cout<<"in real space: "<<std::endl;
  lgf_file<<std::make_pair(G_omega_lattice,beta)<<std::endl;
  std::ofstream lgfks_file("lgfks.dat");
  std::cout<<"in k space: "<<std::endl;
  //G_omega_lattice = transform_into_k_space(G_omega_lattice);
  lgfks_file<<std::make_pair(G_omega_lattice,beta)<<std::endl;
  
  std::cout<<"computing G0 from sigma"<<std::endl;
  for (frequency_t omega=0; omega<nfreq; ++omega) {
    gf_to_matrix(omega, selfenergy, Sigma_mat);
    gf_to_matrix(omega, G_omega_lattice, G_mat);
    InvertMatrix(G_mat, G_inv_mat);
    G0_inv_mat=Sigma_mat+G_inv_mat;
    
    
    //gf_to_matrix(omega,G_omega,G_mat2);
    //InvertMatrix(G_mat2,G_inv_mat2);
    //G0_inv_mat=2*G_inv_mat-G_inv_mat2;
    
    InvertMatrix(G0_inv_mat, G0_mat);
    matrix_to_gf(omega, G0_omega_new, G0_mat);
    
  }
  
  
  std::cout<<"returning from CDMFT self consistency"<<std::endl;
  std::ofstream g0newks_file("g0newks.dat");
  g0newks_file<<std::make_pair(G0_omega_new,beta)<<std::endl;
  return G0_omega_new;
  

  //throw std::runtime_error("nambu operator() in CDMFT transformer not implemented.");
}
/*nambu_green_function<double> LadderCDMFTTransformer::transform_into_k_space(nambu_green_function<double> const& gf) const{
 throw std::runtime_error("transform_into_k_space CDMFT transformer not implemented.");
 }
 nambu_green_function<std::complex<double> > LadderCDMFTTransformer::transform_into_k_space(nambu_green_function<std::complex<double> > const& gf) const{
 throw std::runtime_error("transform_into_k_space CDMFT transformer not implemented.");
 }*/
matsubara_green_function_t LadderCDMFTTransformer::transform_into_real_space(const matsubara_green_function_t& G) const{
  spin_t nflavor=G.nflavor();
  frequency_t nfreq=G.nfreq();
  matsubara_green_function_t G_transformed(G);
  complex_dense_matrix U=U_;
  complex_dense_matrix tmp=U_;
  complex_dense_matrix U_t=boost::numeric::ublas::trans(U_);
  complex_dense_matrix G_mat(n_cluster_sites(),n_cluster_sites()), G_mat_new(n_cluster_sites(),n_cluster_sites());
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (frequency_t omega=0; omega<nfreq; ++omega) {
      gf_to_matrix(omega, flavor, G, G_mat);
      tmp=boost::numeric::ublas::prod(G_mat,U_t);
      G_mat_new=boost::numeric::ublas::prod(U,tmp);
      matrix_to_gf(omega, flavor, G_transformed, G_mat_new);
    }
  }
  return G_transformed;
}

matsubara_green_function_t LadderCDMFTTransformer::transform_into_k_space(const matsubara_green_function_t& G) const{
  spin_t nflavor=G.nflavor();
  frequency_t nfreq=G.nfreq();
  matsubara_green_function_t G_transformed(G);
  complex_dense_matrix U=U_;
  complex_dense_matrix tmp=U_;
  complex_dense_matrix U_t=boost::numeric::ublas::trans(U_);
  complex_dense_matrix G_mat(n_cluster_sites(),n_cluster_sites()), G_mat_new(n_cluster_sites(),n_cluster_sites());
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (frequency_t omega=0; omega<nfreq; ++omega) {
      gf_to_matrix(omega, flavor, G, G_mat);
      tmp=boost::numeric::ublas::prod(G_mat,U);
      G_mat_new=boost::numeric::ublas::prod(U_t,tmp);
      matrix_to_gf(omega, flavor, G_transformed, G_mat_new);
    }
  }
  return G_transformed;
}

itime_green_function_t LadderCDMFTTransformer::transform_into_real_space(const itime_green_function_t& G) const{
  spin_t nflavor=G.nflavor();
  itime_index_t ntime=G.ntime();
  itime_green_function_t G_transformed(G);
  dense_matrix U(U_.size1(), U_.size1());
  for(unsigned int i=0;i<U_.size1();++i){ for (unsigned int j=0;j<U_.size1();++j){ U(i,j)=U_(i,j).real();} }
  dense_matrix tmp=U;
  dense_matrix Ut=boost::numeric::ublas::trans(U);
  dense_matrix G_mat(n_cluster_sites(),n_cluster_sites()), G_mat_new(n_cluster_sites(),n_cluster_sites());
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (itime_index_t tau=0; tau<ntime; ++tau) {
      gf_to_matrix(tau, flavor, G, G_mat);
      tmp=boost::numeric::ublas::prod(G_mat,Ut);
      G_mat_new=boost::numeric::ublas::prod(U,tmp);
      matrix_to_gf(tau, flavor, G_transformed, G_mat_new);
    }
  }
  return G_transformed;
}
itime_green_function_t LadderCDMFTTransformer::transform_into_k_space(const itime_green_function_t& G) const{
  spin_t nflavor=G.nflavor();
  itime_index_t ntime=G.ntime();
  itime_green_function_t G_transformed(G);
  dense_matrix U(U_.size1(), U_.size1());
  for(unsigned int i=0;i<U_.size1();++i){ for (unsigned int j=0;j<U_.size1();++j){ U(i,j)=U_(i,j).real();} }
  dense_matrix tmp=U;
  dense_matrix Ut=boost::numeric::ublas::trans(U);
  dense_matrix G_mat(n_cluster_sites(),n_cluster_sites()), G_mat_new(n_cluster_sites(),n_cluster_sites());
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (itime_index_t tau=0; tau<ntime; ++tau) {
      gf_to_matrix(tau, flavor, G, G_mat);
      tmp=boost::numeric::ublas::prod(G_mat,U);
      G_mat_new=boost::numeric::ublas::prod(Ut,tmp);
      matrix_to_gf(tau, flavor, G_transformed, G_mat_new);
    }
  }
  return G_transformed;
}
nambu_matsubara_green_function_t LadderCDMFTTransformer::transform_into_real_space(const nambu_matsubara_green_function_t& G) const{
  //spin_t nflavor=G.nflavor();
  frequency_t nfreq=G.nfreq();
  nambu_matsubara_green_function_t G_transformed(G);
  complex_dense_matrix U(n_cluster_sites()*2,n_cluster_sites()*2), tmp(n_cluster_sites()*2,n_cluster_sites()*2);
  U.clear();tmp.clear();
  
  for(unsigned int i=0;i<n_cluster_sites();++i){
    for(unsigned int j=0;j<n_cluster_sites();++j){
      U(2*i,2*j)=U(2*i+1,2*j+1)=U_(i,j);
    }
  }
  complex_dense_matrix U_t=boost::numeric::ublas::trans(U);
  complex_dense_matrix G_mat(n_cluster_sites()*2,n_cluster_sites()*2), G_mat_new(n_cluster_sites()*2,n_cluster_sites()*2);
  for (frequency_t omega=0; omega<nfreq; ++omega) {
    gf_to_matrix(omega, G, G_mat);
    tmp=boost::numeric::ublas::prod(G_mat,U_t);
    G_mat_new=boost::numeric::ublas::prod(U,tmp);
    matrix_to_gf(omega,G_transformed, G_mat_new);
  }
  return G_transformed;
  
  //  throw std::runtime_error("this has not yet been implemented.");
}
nambu_matsubara_green_function_t LadderCDMFTTransformer::transform_into_k_space(const nambu_matsubara_green_function_t& G) const{
  //spin_t nflavor=G.nflavor();
  frequency_t nfreq=G.nfreq();
  nambu_matsubara_green_function_t G_transformed(G);
  complex_dense_matrix U(n_cluster_sites()*2,n_cluster_sites()*2), tmp(n_cluster_sites()*2,n_cluster_sites()*2);
  U.clear();tmp.clear();
  for(unsigned int i=0;i<n_cluster_sites();++i){
    for(unsigned int j=0;j<n_cluster_sites();++j){
      U(2*i,2*j)=U(2*i+1,2*j+1)=U_(i,j);
    }
  }
  complex_dense_matrix U_t=boost::numeric::ublas::trans(U);
  complex_dense_matrix G_mat(n_cluster_sites()*2,n_cluster_sites()*2), G_mat_new(n_cluster_sites()*2,n_cluster_sites()*2);
  for (frequency_t omega=0; omega<nfreq; ++omega) {
    gf_to_matrix(omega, G, G_mat);
    tmp=boost::numeric::ublas::prod(G_mat,U);
    G_mat_new=boost::numeric::ublas::prod(U_t,tmp);
    matrix_to_gf(omega,G_transformed, G_mat_new);
  }
  return G_transformed;
  //throw std::runtime_error("this has not yet been implemented.");
}
nambu_itime_green_function_t LadderCDMFTTransformer::transform_into_real_space(const nambu_itime_green_function_t& G) const{
  //spin_t nflavor=G.nflavor();
  itime_index_t ntime=G.ntime();
  nambu_itime_green_function_t G_transformed(G);
  dense_matrix U(n_cluster_sites()*2,n_cluster_sites()*2), tmp(n_cluster_sites()*2,n_cluster_sites()*2);
  U.clear();tmp.clear();
  
  for(unsigned int i=0;i<n_cluster_sites();++i){
    for(unsigned int j=0;j<n_cluster_sites();++j){
      U(2*i,2*j)=U(2*i+1,2*j+1)=real(U_(i,j));
    }
  }
  dense_matrix U_t=boost::numeric::ublas::trans(U);
  dense_matrix G_mat(n_cluster_sites()*2,n_cluster_sites()*2), G_mat_new(n_cluster_sites()*2,n_cluster_sites()*2);
  for (itime_index_t tau=0; tau<ntime; ++tau) {
    gf_to_matrix(tau, G, G_mat);
    tmp=boost::numeric::ublas::prod(G_mat,U_t);
    G_mat_new=boost::numeric::ublas::prod(U,tmp);
    matrix_to_gf(tau,G_transformed, G_mat_new);
  }
  return G_transformed;
  
  //  throw std::runtime_error("this has not yet been implemented.");
}
nambu_itime_green_function_t LadderCDMFTTransformer::transform_into_k_space(const nambu_itime_green_function_t& G) const{
  //spin_t nflavor=G.nflavor();
  itime_index_t ntime=G.ntime();
  nambu_itime_green_function_t G_transformed(G);
  dense_matrix U(n_cluster_sites()*2,n_cluster_sites()*2), tmp(n_cluster_sites()*2,n_cluster_sites()*2);
  U.clear();tmp.clear();
  
  for(unsigned int i=0;i<n_cluster_sites();++i){
    for(unsigned int j=0;j<n_cluster_sites();++j){
      U(2*i,2*j)=U(2*i+1,2*j+1)=real(U_(i,j));
    }
  }
  dense_matrix U_t=boost::numeric::ublas::trans(U);
  dense_matrix G_mat(n_cluster_sites()*2,n_cluster_sites()*2), G_mat_new(n_cluster_sites()*2,n_cluster_sites()*2);
  for (itime_index_t tau=0; tau<ntime; ++tau) {
    gf_to_matrix(tau, G, G_mat);
    tmp=boost::numeric::ublas::prod(G_mat,U);
    G_mat_new=boost::numeric::ublas::prod(U_t,tmp);
    matrix_to_gf(tau,G_transformed, G_mat_new);
  }
  return G_transformed;
  
  //  throw std::runtime_error("this has not yet been implemented.");
}
void LadderCDMFTTransformer::set_Tij_matrix(double kx, double ky, complex_dense_matrix &H0){
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
void LadderCDMFTTransformer::set_nambu_Tij_dwave_matrix(double kx, double ky, double eta1, double eta2, double mu, double h, complex_dense_matrix &H0){
  int i,j,nsite;
  nsite=n_cluster_sites();
  
  double coskx=cos(kx);
  double cosky=cos(ky);
  
  std::complex<double> exp_ikx=std::exp(std::complex<double>(0, kx));
  std::complex<double> exp_mikx=std::conj(exp_ikx);
  std::complex<double> exp_iky=std::exp(std::complex<double>(0, ky));
  std::complex<double> exp_miky=std::conj(exp_iky);
  
  std::complex<double> exp_2ikx=std::exp(std::complex<double>(0, 2*kx));
  std::complex<double> exp_m2ikx=std::conj(exp_2ikx);
  std::complex<double> exp_2iky=std::exp(std::complex<double>(0, 2*ky));
  std::complex<double> exp_m2iky=std::conj(exp_2iky);
  
  
  std::complex<double> eta=std::complex<double>(-eta1,eta2);
  complex_dense_matrix G(n_cluster_sites(),n_cluster_sites()),mconj_G(n_cluster_sites(),n_cluster_sites()),F(n_cluster_sites(),n_cluster_sites()),F_conj(n_cluster_sites(),n_cluster_sites());
  G.clear();F.clear();mconj_G.clear();F_conj.clear();
  // G(0,0)=mu;  G(1,1)=mu; G(2,2)=mu; G(3,3)=mu; ///epsilon(k)->-2t(coskx+cosky) so mu-epsilon=mu+2t(coskx+cosky)  +h -h is the applied antiferromagnetic field
  //If the second term is deleted, it will be back to the homogeneous problem t=tstar
  G(0,1)=2.*t_*cosky*exp_miky+(tstar_-t_)*exp_m2iky;
  G(0,3)=2.*t_*coskx*exp_mikx+(tstar_-t_)*exp_m2ikx;
  G(1,0)=2.*t_*cosky*exp_iky+(tstar_-t_)*exp_2iky;
  G(1,2)=2.*t_*coskx*exp_mikx+(tstar_-t_)*exp_m2ikx;
  G(2,1)=2.*t_*coskx*exp_ikx+(tstar_-t_)*exp_2ikx;
  G(2,3)=2.*t_*cosky*exp_iky+(tstar_-t_)*exp_2iky;
  G(3,0)=2.*t_*coskx*exp_ikx+(tstar_-t_)*exp_2ikx;
  G(3,2)=2.*t_*cosky*exp_miky+(tstar_-t_)*exp_m2iky;
  
  mconj_G=-boost::numeric::ublas::conj(G);
  //mconj_G=-G;
  F(0,1)=-cosky*exp_miky*eta;
  F(0,3)= coskx*exp_mikx*eta;
  F(1,0)=-cosky*exp_iky*eta;
  F(1,2)= coskx*exp_mikx*eta;
  F(2,1)= coskx*exp_ikx*eta;
  F(2,3)=-cosky*exp_iky*eta;
  F(3,0)= coskx*exp_ikx*eta;
  F(3,2)=-cosky*exp_miky*eta;
  
  F_conj=boost::numeric::ublas::conj(F);
  
  for(i=0;i<nsite;i++){
    for(j=0;j<nsite;j++){
      H0(2*i  ,2*j  )=G(i,j);
      H0(2*i+1,2*j+1)=mconj_G(i,j);
      H0(2*i  ,2*j+1)=F(i,j);
      H0(2*i+1,2*j  )=F_conj(i,j);
    }
  }
  H0(0,0)+=mu+h/2;H0(1,1)+=-mu+h/2;H0(4,4)+=mu+h/2;H0(5,5)+=-mu+h/2;
  H0(2,2)+=mu-h/2;H0(3,3)+=-mu-h/2;H0(6,6)+=mu-h/2;H0(7,7)+=-mu-h/2;
  
}

void LadderCDMFTTransformer::set_U_(unsigned int nsite, complex_dense_matrix &U_){
  double phase=0;
  U_.resize(nsite,nsite);
  for(int i=0;i<nsite;i++){
    for(int j=0;j<nsite;j++){
      phase=0;
      for (unsigned int d=0; d<dimension(); ++d){
        phase += cluster_coordinate(i,d)*cluster_momentum(j,d);
        // std::cout<<i<<j<<" "<<phase<<std::endl;
      }
      U_(i,j)=std::exp(std::complex<double>(0, phase));
      //std::cout<<i<<j<<" "<<U_(i,j)<<std::endl;
    }
  }
  U_/=sqrt(n_cluster_sites());
  for(int i=0;i<nsite;i++){
    for(int j=0;j<nsite;j++){
      std::cout<<U_(i,j)<<" ";
    }
    std::cout<<std::endl;
  }

}

std::vector<double> LadderCDMFTTransformer::epsav() const{
  return epsav_;
}
std::vector<double> LadderCDMFTTransformer::epssqav() const{
  return epssqav_;
}

template<class T>
green_function<T> LadderCDMFTTransformer::spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const
{
  green_function<T> G_new(G);
  if (paramagnet) {
    if(G_new.nflavor()>1 && G_new.nflavor()%2==0) {
      for(spin_t flavor=0; flavor<G_new.nflavor(); flavor+=2) {
        for (site_t k1=0; k1<n_cluster_sites(); ++k1) {
          for (site_t k2=0; k2<n_cluster_sites(); ++k2) {
            for(frequency_t freq=0; freq<G_new.nfreq(); ++freq) {
              G_new(freq, k1, k2, flavor)= 0.5*(G_new(freq, k1, k2, flavor)+G_new(freq, k1, k2, flavor+1));
              G_new(freq, k1, k2, flavor+1)=G_new(freq, k1, k2, flavor);
            }
          }
        }
      }
    }
  }
  return G_new;
}
nambu_green_function<std::complex<double> > LadderCDMFTTransformer::spin_symmetrizer(const nambu_green_function<std::complex<double> >& G, const bool paramagnet) const
{
  nambu_green_function<std::complex<double> > G_new(G);
  //std::ofstream nambu_before_file("beforespin.dat");
  //nambu_before_file<<std::make_pair(G, 10.)<<std::endl;;
  if (paramagnet) {
    for (site_t k1=0; k1<n_cluster_sites(); ++k1) {
      for (site_t k2=0; k2<n_cluster_sites(); ++k2) {
        for(frequency_t freq=0; freq<G_new.nfreq(); ++freq) {
          G_new(freq, k1, k2, 0,0).imag()= 0.5*(G(freq, k1, k2, 0,0).imag()+G_new(freq, k1, k2, 1,1).imag());
          G_new(freq, k1, k2, 0,0).real()= 0.5*(G(freq, k1, k2, 0,0).real()-G_new(freq, k1, k2, 1,1).real());
          G_new(freq, k1, k2, 1,1).real()= -G_new(freq, k1, k2, 0,0).real();
          G_new(freq, k1, k2, 1,1).imag()=  G_new(freq, k1, k2, 0,0).imag();
          //symmetrize the  F.
          G_new(freq, k1, k2, 0,1)= (G(freq, k1, k2, 0,1)+std::conj(G(freq, k1, k2, 1,0)))/2.;
          G_new(freq, k1, k2, 1,0)= std::conj(G_new(freq, k1, k2, 0,1));
        }
        //std::cout<<"symmetrizing: "<<G(0, k1, k2, 0,0)<<" "<<G(0, k1, k2, 1,1)<<" and: "<<G(0, k1, k2, 0,1)<<" "<<G(0, k1, k2, 1,0)<<std::endl;
        //std::cout<<"new values  : "<<G_new(0, k1, k2, 0,0)<<" "<<G_new(0, k1, k2, 1,1)<<" and: "<<G_new(0, k1, k2, 0,1)<<" "<<G_new(0, k1, k2, 1,0)<<std::endl;
      }
    }
  }
  else{
    for (site_t k1=0; k1<n_cluster_sites(); ++k1) {
      for (site_t k2=0; k2<n_cluster_sites(); ++k2) {
        int k1_p,k2_p;
        
        if(abs(k1-k2)%2!=0){//neighbouring site
          k1_p=k2;k2_p=k1;
        }
        else{
          k1_p=k1+1;k2_p=k2+1;
          if(k1_p==n_cluster_sites())
            k1_p-=n_cluster_sites();
          if(k2_p==n_cluster_sites())
            k2_p-=n_cluster_sites();
        }
        for(frequency_t freq=0; freq<G_new.nfreq(); ++freq) {
          G_new(freq, k1, k2, 0,0).imag()= 0.5*(G(freq, k1, k2, 0,0).imag()+G_new(freq, k1_p, k2_p, 1,1).imag());
          G_new(freq, k1, k2, 0,0).real()= 0.5*(G(freq, k1, k2, 0,0).real()-G_new(freq, k1_p, k2_p, 1,1).real());
          G_new(freq, k1, k2, 1,1).imag()= 0.5*(G(freq, k1, k2, 1,1).imag()+G_new(freq, k1_p, k2_p, 0,0).imag());
          G_new(freq, k1, k2, 1,1).real()= 0.5*(G(freq, k1, k2, 1,1).real()-G_new(freq, k1_p, k2_p, 0,0).real());
          //symmetrize the  F.
          G_new(freq, k1, k2, 0,1)= (G(freq, k1, k2, 0,1)+std::conj(G(freq, k1_p, k2_p, 1,0)))/2.;
          G_new(freq, k1, k2, 1,0)= (G(freq, k1, k2, 1,0)+std::conj(G(freq, k1_p, k2_p, 0,1)))/2.;
        }
      }
    }
  }
  //std::ofstream nambu_after_file("afterspin.dat");
  //nambu_after_file<<std::make_pair(G_new, 10.);
  return G_new;
}

