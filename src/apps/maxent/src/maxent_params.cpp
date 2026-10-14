/*
 * Copyright (C) 1998-2018 ALPS Collaboration.
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

// Legacy ALPS data and covariance handling from src/apps/maxent/src/maxent_parms.cpp
// at d4533fd30c65656cf5ab7fa448dd6b921b39aa6f; original notice:
/*****************************************************************************
 *
 * ALPS Project Applications
 *
 * Copyright (C) 2010 by Sebastian Fuchs <fuchs@comp-phys.org>
 *                       Thomas Pruschke <pruschke@comp-phys.org>
 *                       Matthias Troyer <troyer@comp-phys.org>
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

// Adapted from CQMP/Maxent src/maxent_params.cpp, revision
// a3992534228b0de732e033ae23154969724814be, on 2026-10-14.
// Uses ALPS parameter APIs and retains ALPS input/output conventions;
// see ../README.md for provenance and licensing.

#include "maxent.hpp"
#include "eigen_lapack.hpp"
#include <alps/config.h> // needed to set up correct bindings
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>
#include <alps/hdf5/vector.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/legendre.hpp> //needed for Legendre transform
namespace bmth = boost::math;

// ALPS text data consists of indexed rows: index value error. Complex
// frequency data uses alternating real and imaginary rows. Imaginary times
// come from TAU_i parameters or the uniform grid, independently of the index.
void ContiParameters::read_data_from_text_file(const alps::params& p) {
  const std::string fname = p["DATA"].cast<std::string>();
  std::ifstream datstream(fname.c_str());
  if (!datstream)
    throw std::invalid_argument("could not open data text file: " + fname);

  std::vector<bool> loaded(ndat(), false);
  const double norm = p["NORM"];
  int index;
  double value, error;
  while (datstream >> index >> value >> error) {
    if (index < 0)
      throw std::invalid_argument("negative data index in: " + fname);
    if (index < ndat()) {
      y_(index) = value / norm;
      sigma_(index) = error / norm;
      loaded[index] = true;
  }
      }
  if (!datstream.eof())
    throw std::invalid_argument("expected data rows: index value error in " + fname);
  if (std::find(loaded.begin(), loaded.end(), false) != loaded.end())
    throw std::invalid_argument("data file has fewer than NDAT indexed values: " + fname);

  if (!p["COVARIANCE_MATRIX"].cast<std::string>().empty())
    read_covariance_matrix_from_text_file(p["COVARIANCE_MATRIX"].cast<std::string>());
}

///Read data from a hdf5 file, with filename given by p["DATA"] in the parameters.
///The data is stored at /Data.
///The error is stored at /Error
///if the parameter COVARIANCE_MATRIX is specified, then the covariance matrix is read in instead of the error.
///The covariance matrix is expected to be stored at /Covariance

void ContiParameters::read_data_from_hdf5_file(const alps::params& p) {
  std::string fname = p["DATA"];
  //attempt to read from h5 archive
  alps::hdf5::archive ar(fname, "r");
  std::vector<double> tmp(ndat());
  std::stringstream path;
  path << "/Data";
  ar >> alps::make_pvp(path.str(), tmp);
  for (std::size_t i = 0; i < ndat(); i++)
    y_(i) = tmp[i] / static_cast<double>(p["NORM"]);
  path.str("");
  if (p["COVARIANCE_MATRIX"].cast<std::string>()=="") {
    path << "/Error";
    ar >> alps::make_pvp(path.str(), tmp);
    for (std::size_t i = 0; i < ndat(); i++)
      sigma_(i) = tmp[i] / static_cast<double>(p["NORM"]);
  } else {
    path << "/Covariance";
    cov_.resize(ndat(), ndat());
    tmp.clear();
    tmp.resize(ndat() * ndat());
    ar >> alps::make_pvp(path.str(), tmp);
    for (std::size_t i = 0; i < ndat(); i++)
      for (std::size_t j = 0; j < ndat(); j++)
        cov_(i, j) = tmp[i * ndat() + j];
  }
}

void ContiParameters::read_data_from_param_file(const alps::params& p) {
  if (!p.defined("NORM")) {
    throw std::runtime_error("parameter NORM missing!");
  } else
    std::cerr << "Data normalized to: " << static_cast<double>(p["NORM"])
        << std::endl;

  for (int i = 0; i < ndat(); ++i) {
    if (!p.defined("X_" + boost::lexical_cast<std::string>(i))) {
      throw std::runtime_error("parameter X_"+ boost::lexical_cast<std::string>(i)+ " missing!");
    }
    y_(i) = static_cast<double>(p["X_" + boost::lexical_cast<std::string>(i)])
        / static_cast<double>(p["NORM"]);
    if (p["COVARIANCE_MATRIX"].cast<std::string>()=="") {
      if (!p.defined("SIGMA_" + boost::lexical_cast<std::string>(i))) {
        throw std::runtime_error(
            std::string("parameter SIGMA_"+boost::lexical_cast<std::string>(i)+ " missing! "));
      }
      sigma_(i) = static_cast<double>(p["SIGMA_"
          + boost::lexical_cast<std::string>(i)])
          / static_cast<double>(p["NORM"]);
    }
  }
}

ContiParameters::ContiParameters(alps::params& p) :
T_(p["T"] | (1/p["BETA"].cast<double>())),ndat_(p["NDAT"]), nfreq_(p["NFREQ"]),
y_(ndat_),sigma_(ndat_),K_(),grid_(p),inputGrid_(ndat_)
{
  //note: T_=1/beta now taken care of elsewhere
  if (ndat_<4) 
    boost::throw_exception(std::invalid_argument("NDAT too small"));


  if (p.defined("DATA") && p["DATA"].cast<std::string>() != "") {
    if(p.defined("DATA_IN_HDF5") && (p["DATA_IN_HDF5"])) {
      //attempt to read from h5 archive
      read_data_from_hdf5_file(p);
    } else {
      read_data_from_text_file(p);
    }
  } else {
    // ALPS parameters already contain the typed X_i and SIGMA_i values.
    read_data_from_param_file(p);
    if (p["COVARIANCE_MATRIX"].cast<std::string>() != "")
      read_covariance_matrix_from_text_file(p["COVARIANCE_MATRIX"].cast<std::string>());
  }
}

void ContiParameters::scale_data_with_error(const int ntab) {
  //Look around Eq. D.5 in Sebastian's thesis. We have sigma_ = sqrt(eigenvalues of covariance matrix) or, 
  //in case of a diagonal covariance matrix, we have sigma_=SIGMA_X. 
  //Then define y := \bar{G}/sigma_ and K := (1/sigma_)\tilde{K}
  for (int i = 0; i < ndat(); i++) {
    y_[i] /= sigma_[i];
    for (int j = 0; j < ntab; ++j) {
      K_(i, j) /= sigma_[i];
    }
  }
}

void ContiParameters::read_covariance_matrix_from_text_file(
    const std::string& fname) {
  cov_ = matrix_type::Zero(ndat(), ndat());
  std::cerr << "Reading covariance matrix\n";
  std::ifstream covstream(fname.c_str());
  if (!covstream)
    boost::throw_exception(
        std::invalid_argument(
            "could not open covariance matrix file: " + fname));

  int i, j;
  double covariance;
  while (covstream >> i >> j >> covariance) {
    if (i < 0 || j < 0)
      throw std::invalid_argument("negative covariance index in: " + fname);
    if (i < ndat() && j < ndat())
      cov_(i, j) = covariance;
  }
  if (!covstream.eof())
    throw std::invalid_argument("expected covariance rows: i j value in " + fname);
}

void ContiParameters::decompose_covariance_matrix(const alps::params& p){
    vector_type var(ndat());
    //bindings::lapack::syev('V', bindings::upper(cov_) , var, bindings::lapack::optimal_workspace()); 
    //TODO: check if this truly implements lapack's expected overwrite of cov_
    // ALPS's LAPACK covariance decomposition uses the upper triangle.
    Eigen::SelfAdjointEigenSolver<matrix_type> es(cov_.selfadjointView<Eigen::Upper>());
    var=es.eigenvalues();
    cov_=es.eigenvectors();
    matrix_type cov_trans = cov_.transpose();
    matrix_type K_loc = maxent_prec_prod(cov_trans, K_);
    vector_type y_loc = maxent_prec_prod(cov_trans, y_);
    if (p["VERBOSE"])
      std::cout << "# Eigenvalues of the covariance matrix:\n";
    // We drop eigenvalues of the covariance matrix which are smaller than 1e-10
    // as they represent bad data directions (usually there is a steep drop
    // below that value)
    double cutoff = p["CM_EIGENVALUE_CUTOFF"];
    std::cout << "Cutoff for singular eigenvalues of the covariance: " << cutoff << "\n";
    int new_ndat_,old_ndat_=ndat();
    for (new_ndat_ =0;new_ndat_<ndat();new_ndat_++)
      if (var[new_ndat_]>cutoff) break; 
    // This is the number of good data
    ndat_ = old_ndat_ - new_ndat_;
    if (new_ndat_ ==0)
      std::cout << " Ignoring no singular eigenvalues of covariance" << std::endl;
    else
      std::cout << "# Ignoring singular eigenvalues (0-" << new_ndat_-1 << " out of " << old_ndat_ << ")\n";
    // Now resize kernel and data matrix and fill it with the values for the
    // good data directions
    K_.resize(ndat_,nfreq_);
    y_.resize(ndat_);
    sigma_.resize(ndat_);
    for (int i=0; i<ndat(); i++) {
      y_(i) = y_loc(new_ndat_+i);
      for (int j=0; j<nfreq_; j++) {
        K_(i,j)=K_loc(new_ndat_+i,j);
      }
    }
    for (int i=0; i<ndat(); ++i) {
      sigma_[i] = std::sqrt(std::abs(var(new_ndat_+i)))/static_cast<double>(p["NORM"]);
      if (p["VERBOSE"])
        std::cout << "# " << var(new_ndat_+i) << "\n";
    }
}

void MaxEntParameters::compute_minimal_chi2()const {
  matrix_type Ut = U_.transpose(); //U^T has dimension ns_*ndat()
  vector_type t = maxent_prec_prod(Ut, y_); //t has dimension ns_
  vector_type y2 = maxent_prec_prod(U_, t); //y2 has dimension ndat(), which is dimension of y
  double chi = (y_ - y2).norm(); //this measures the loss of precision when transforming to singular space and back.
  std::cout << "minimal chi2: " << chi * chi / y_.size() << std::endl;
}

void MaxEntParameters::truncate_to_singular_space(const vector_type& S) {
  //(truncated) U has dimension ndat() * ns_; ndat() is # of input (matsubara frequency/imag time) points
  //(truncated) Sigma has dimension ns_*ns_ (number of singular eigenvalues)
  //(truncated) V^T has dimensions ns_* nfreq(); nfreq() is number of output (real) frequencies
  //ns_ is the dimension of the singular space.
  U_.conservativeResize(ndat(), ns_);
  Vt_.conservativeResize(ns_, nfreq());
  Sigma_= Eigen::MatrixXd::Zero(ns_,ns_);
  for (int s = 0; s < ns_; ++s) {
    std::cout << "# " << s << "\t" << S[s] << "\n";
    Sigma_(s, s) = S[s];
  }
}

void MaxEntParameters::singular_value_decompose_kernel(bool verbose,
    vector_type& S) {
  /*boost::numeric::bindings::lapack::gesvd('S', 'S', Kt, S, U_, Vt_);
  */

#ifdef HAVE_LAPACK

  matrix_type Kt = K_; // gesvd destroys K!
  lapack_svd(Kt,S, Vt_,U_);

#else
  const double threshold = std::sqrt(std::numeric_limits<double>::epsilon())
      * nfreq();
  Eigen::JacobiSVD<matrix_type> svd(K_,Eigen::ComputeThinU | Eigen::ComputeThinV);
  //svd.setThreshold(threshold);
  S=svd.singularValues();
  U_=svd.matrixU();
  Vt_=svd.matrixV().transpose(); 
#endif
  if (verbose)
    std::cout << "# Singular values of the Kernel:\n";

  const double prec = std::sqrt(std::numeric_limits<double>::epsilon())
      * nfreq() * S[0];
  if (verbose)
    std::cout << "# eps = " << sqrt(std::numeric_limits<double>::epsilon())
        << std::endl << "# prec = " << prec << std::endl;

  for (unsigned int s = 0; s < S.size(); ++s) {
    if (verbose)
      std::cout << "# " << s << "\t" << S[s] << "\n";

    ns_ = (S[s] >= prec) ? s + 1 : ns_;
  }
  if (ns() == 0)
    boost::throw_exception(
        std::logic_error("all singular values smaller than the precision"));
  
}
void MaxEntParameters::check_high_frequency_limit(const vector_type& y,const kernel_type kt){
    //TODO: expand check into time/bosonic domain
    //we know that the limit of a green's function is ~ -1/iw_{n}
    //this checks that we have a good high frequency limit
    //within the tolerance of the error bar
    
    int n=0;
    if (kt==frequency_fermionic_kernel)
        n = (ndat()-2)/2; //this is the real n; due to weird vector structure
    if (kt==frequency_fermionic_ph_kernel)
        n=ndat()-1;
    
    if(n!=0){
        //limit = G(iwn)*iwn
        std::complex<double> iwn(0,(2*n+1)*M_PI*T());
        double limit = y(ndat()-1)*iwn.imag();
        if(std::abs(1+limit)>sigma(ndat()-1)){
            std::cerr<<"The high frequency limit is not 1!: " << std::abs(limit)
            <<" Check norm?"<< std:: endl;
        }
        //now backcontinue default model and check high frequency limit
        // G(iw_{n})=\sum_{m}K_{nm}A_{m}
        std::complex<double> G;
        iwn = std::complex<double>(iwn.real(), (2*2000+1)*M_PI*T());
        for(int j=0;j<nfreq();j++){
            G+= 1.0/(iwn-omega_coord(j))*Default().D(omega_coord(j)) * delta_omega(j);
        }
        limit = G.imag()*iwn.imag();
        if(std::abs(1+limit)>.01){
            std::cerr<<"The high frequency limit of the default model is not 1!: "
            << std::abs(limit) <<" Check norm?"<< std:: endl;
        }
    }
    //time space=> tail_1 = -G(0)-G(beta) = 1
    if(kt==time_fermionic_kernel){
        double err = sqrt(sigma(0)*sigma(0)+sigma(ndat()-1)*sigma(ndat()-1));
        double limit = -y(0)-y(ndat()-1);
        if(std::abs(limit-1)>err)
            std::cerr<<"The high frequency limit is not +1!: " << limit
            <<" Check norm?"<< std:: endl;
    }
}

MaxEntParameters::MaxEntParameters(alps::params& p) :
    ContiParameters(p),
    Default_(make_default_model(p, "DEFAULT_MODEL")),
    U_(ndat(), ndat()), Vt_(ndat(), nfreq()), Sigma_(ndat(), ndat()),
    omega_coord_(nfreq()), delta_omega_(nfreq()), ns_(0)
{
  for (int i=0; i<nfreq(); ++i) {
    omega_coord_[i] = (Default().omega_of_t(grid_(i)) + Default().omega_of_t(grid_(i+1)))/2.;
    delta_omega_[i] = Default().omega_of_t(grid_(i+1)) - Default().omega_of_t(grid_(i));
  }
  //build a kernel matrix
  kernel ker(p,omega_coord_,inputGrid_);
  K_=ker();
  k_type = ker.getKernelType();

  //scale lhs and rhs according to errors, etc.
  if (p["COVARIANCE_MATRIX"].cast<std::string>()!="")
    decompose_covariance_matrix(p);
    
    check_high_frequency_limit(y(),k_type);

  //Look around Eq. D.5 in Sebastian's thesis. We have sigma_ = sqrt(eigenvalues of covariance matrix) or, 
  //in case of a diagonal covariance matrix, we have sigma_=SIGMA_X. 
  //The then define y := \bar{G}/sigma_ and K := (1/sigma_)\tilde{K}
  scale_data_with_error(nfreq());

  // Preserve ALPS's optional normalization constraint on the final data row.
  if (p["ENFORCE_NORMALIZATION"] | false) {
    const double sigma = p["SIGMA_NORMALIZATION"].cast<double>() / p["NORM"].cast<double>();
    for (int j=0; j<nfreq(); ++j)
      K_(ndat()-1,j) = 1. / sigma;
    y_[ndat()-1] = 1. / sigma;
  }

  std::cerr << "Kernel is set up\n";

  vector_type S(ndat());
  bool verbose=p["VERBOSE"];
  //perform the SVD decomposition K = U Sigma V^T
  singular_value_decompose_kernel(verbose, S);

  //truncate all matrix to singular space
  truncate_to_singular_space(S);

  //compute Ut and 
  compute_minimal_chi2();
}
