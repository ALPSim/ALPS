/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2010 by Sebastian Fuchs <fuchs@comp-phys.org>
*                       Thomas Pruschke <pruschke@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
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

#include "maxent_parms.hpp"
#include <boost/numeric/bindings/traits/ublas_symmetric.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/lapack/gesvd.hpp>


ContiParameters::ContiParameters(const alps::Parameters& p) : 
  Default_(make_default_model(p, "DEFAULT_MODEL")), 
  T_(p.value_or_default("T", 1./static_cast<double>(p["BETA"]))),
  ndat_(p["NDAT"]), nfreq_(p["NFREQ"]),
  y_(ndat_), K_(),   t_array_(nfreq_+1)
{
  if (ndat_<4) 
    boost::throw_exception(std::invalid_argument("NDAT too small"));
  if (p["FREQUENCY_GRID"]=="Lorentzian") {
    double cut = p.value_or_default("CUT", 0.01);
    std::vector<double> temp(nfreq_+1);
    for (int i=0; i<nfreq_+1; ++i) 
      temp[i] = tan(M_PI * (double(i)/(nfreq_)*(1.-2*cut)+cut - 0.5));
    for (int i=0; i<nfreq_+1; ++i) 
      t_array_[i] = (temp[i] - temp[0])/(temp[temp.size()-1] - temp[0]);
  }
  else if (p["FREQUENCY_GRID"]=="half Lorentzian") {
    double cut = p.value_or_default("CUT", 0.01);
    std::vector<double> temp(nfreq_+1);
    for (int i=0; i<nfreq_+1; ++i) 
      temp[i] = tan(M_PI * (double(i+nfreq_)/(2*nfreq_)*(1.-2*cut)+cut - 0.5));
    for (int i=0; i<nfreq_+1; ++i) 
      t_array_[i] = (temp[i] - temp[0])/(temp[temp.size()-1] - temp[0]);
  }
  else if (p["FREQUENCY_GRID"]=="quadratic") {
    double s = p.value_or_default("SPREAD", 4);
    if (s<1) 
      boost::throw_exception(std::invalid_argument("the parameter SPREAD must be greater than 1"));
    std::vector<double> temp(nfreq_);
    double t=0;
    for (int i=0; i<nfreq_; ++i) {
      double a = double(i)/(nfreq_-1);
      double factor = 4*(s-1)*(a*a-a)+s;
      factor /= double(nfreq_)/(3.*(nfreq_-1))*((nfreq_)*(2+s)-4+s);
      double delta_t = factor;
      t += delta_t;
      temp[i] = t;
    }
    t_array_[0] = 0.;
    for (int i=1; i<nfreq_+1; ++i) 
      t_array_[i]  = temp[i-1]/temp[temp.size()-1];
  }
  else if (p["FREQUENCY_GRID"]=="linear") {
    for (int i=0; i<nfreq_+1; ++i) 
      t_array_[i] = double(i)/(nfreq_);
  }
  else 
    boost::throw_exception(std::invalid_argument("No valid frequency grid specified"));
  for (int i=0; i<ndat(); ++i) 
    y_(i) = static_cast<double>(p["X_"+boost::lexical_cast<std::string>(i)]);
}




void ContiParameters::setup_kernel(const alps::Parameters& p, const int ntab, const vector_type& freq)
{
  using namespace boost::numeric;
  K_.resize(ndat_, ntab);
  for (int i=0; i<ndat(); ++i) 
    y_(i) = static_cast<double>(p["X_"+boost::lexical_cast<std::string>(i)]);
  if(p["DATASPACE"]=="time") { 
    if (alps::is_master())
      std::cerr << "assume time space data" << std::endl;
    if (p["KERNEL"] == "fermionic") {
      if (alps::is_master())
        std::cerr << "Using fermionic kernel" << std::endl;
      for (int i=0; i<ndat(); ++i) {
    double tau;
    if (p.defined("TAU_"+boost::lexical_cast<std::string>(i)))
      tau = p["TAU_"+boost::lexical_cast<std::string>(i)]; 
    else 
      tau = i / (ndat() * T_);
        for (int j=0; j<ntab; ++j) {
          double omega = freq[j]; //Default().omega_of_t(double(j)/(ntab-1));
          K_(i,j) =  -1. / (std::exp(omega*tau) + std::exp(-omega*(1./T_-tau)));
        }
      }
    }
    else if (p["KERNEL"] == "bosonic") {
      if (alps::is_master())
        std::cerr << "Using bosonic kernel" << std::endl;
      for (int i=0; i<ndat(); ++i) {
    double tau;
    if (p.defined("TAU_"+boost::lexical_cast<std::string>(i)))
      tau = p["TAU_"+boost::lexical_cast<std::string>(i)]; 
    else 
      tau = i / (ndat() * T_);
        K_(i,0) = T_;
        for (int j=1; j<ntab; ++j) {
          double omega = freq[j];
          K_(i,j) = 0.5*omega * (std::exp(-omega*tau) + std::exp(-omega*(1./T_-tau))) / (1 - std::exp(-omega/T_));
        }
      }    
    }
    else if (p["KERNEL"] == "Boris") {
      if (alps::is_master())
        std::cerr << "Using Boris' kernel" << std::endl;
      for (int i=0; i<ndat(); ++i) {
    double tau = p["TAU_"+boost::lexical_cast<std::string>(i)]; 
        for (int j=0; j<ntab; ++j) {
          double omega = freq[j];
          K_(i,j) = std::exp(-omega*tau);
        }
      }    
    }
    else 
      boost::throw_exception(std::invalid_argument("unknown integration kernel"));
  } 
  else if (p["DATASPACE"]=="frequency" && p["KERNEL"] == "fermionic" &&
           p.value_or_default("PARTICLE_HOLE_SYMMETRY", false)) {
    std::cerr << "using particle hole symmetric fermionic frequency space kernel" << std::endl;
    for (int i=0; i<ndat(); ++i) {
      double omegan = (2*i+1)*M_PI*T_;
      for (int j=0; j<ntab; ++j) {
        double omega = freq[j]; 
        K_(i,j) =  -omegan / (omegan*omegan + omega*omega);
      }
    }
  } 
  else if (p["DATASPACE"]=="frequency") { 
    if (alps::is_master())
      std::cerr << "assume frequency space data" << std::endl;
    ublas::matrix<std::complex<double>, ublas::column_major> Kc(ndat_/2, ntab);
    if (p["KERNEL"] == "fermionic") {
      if (alps::is_master())
        std::cerr << "Using fermionic kernel" << std::endl;
      for (int i=0; i<ndat()/2; ++i) {
        std::complex<double> iomegan(0, (2*i+1)*M_PI*T_);
        for (int j=0; j<ntab; ++j) {
          double omega = freq[j]; 
          Kc(i,j) =  1. / (iomegan - omega);
        }
      }
    }
    else if (p["KERNEL"] == "bosonic") {
      if (alps::is_master())
        std::cerr << "Using bosonic kernel" << std::endl;
      for (int i=0; i<ndat()/2; ++i) {
        std::complex<double> iomegan(0, 2*i*M_PI*T_);
        for (int j=1; j<ntab; ++j) {
          double omega = freq[j]; 
          Kc(i,j) =  -1. / (iomegan - omega);
        }
      }    
    }
    else 
      boost::throw_exception(std::invalid_argument("unknown integration kernel"));    
    for (int i=0; i<ndat(); i+=2) {
      for (int j=1; j<ntab; ++j) {
        K_(i,j) = Kc(i/2,j).real();
        K_(i+1,j) = Kc(i/2,j).imag();
      }
    }
  }
  else
    boost::throw_exception(std::invalid_argument("unknown value for parameter DATASPACE"));
  vector_type sigma(ndat());
  if (p.defined("COVARIANCE_MATRIX")) {
    matrix_type cov(ndat(), ndat());
    if (alps::is_master())
      std::cerr << "Reading covariance matrix\n";
    std::ifstream covstream(static_cast<std::string>(p["COVARIANCE_MATRIX"]).c_str());
    if (!covstream)
      boost::throw_exception(std::invalid_argument("could not open covariance matrix file: "+p["COVARIANCE_MATRIX"]));
    int i, j;
    double covariance;
    while (covstream) {
      covstream >> i >> j >> covariance;
      if (i<ndat() && j<ndat())
        cov(i,j) = covariance;
    }
    vector_type var(ndat());
    bindings::lapack::syev('V', 'U', cov , var, bindings::lapack::optimal_workspace());
    matrix_type cov_trans = ublas::trans(cov);
    K_ = ublas::prec_prod(cov_trans, K_);
    y_ = ublas::prec_prod(cov_trans, y_);
    if (alps::is_master()) 
      std::cout << "# Eigenvalues of the covariance matrix:\n";
    for (int i=0; i<ndat(); ++i) { 
      sigma[i] = std::sqrt(var(i));
      if (alps::is_master())
        std::cout << "# " << var(i) << "\n";
    }
  } 
  else {
    for (int i=0; i<ndat(); ++i) 
      sigma[i] = static_cast<double>(p["SIGMA_"+boost::lexical_cast<std::string>(i)]);
  }
  for (int i=0; i<ndat(); ++i) {
    y_[i] /= sigma[i];
    for (int j=0; j<ntab; ++j) 
      K_(i,j) /= sigma[i];
  }
}




MaxEntParameters::MaxEntParameters(const alps::Parameters& p) :
  ContiParameters(p), 
  U_(ndat(), ndat()), Vt_(ndat(), nfreq()), Sigma_(ndat(), ndat()), 
  omega_coord_(nfreq()), delta_omega_(nfreq()), ns_(0)
{
  using namespace boost::numeric;
  if (ndat() > nfreq()) 
    boost::throw_exception(std::invalid_argument("NDAT should be smaller than NFREQ"));
  for (int i=0; i<nfreq(); ++i) {
    omega_coord_[i] = (Default().omega_of_t(t_array_[i]) + Default().omega_of_t(t_array_[i+1]))/2.;
    delta_omega_[i] = Default().omega_of_t(t_array_[i+1]) - Default().omega_of_t(t_array_[i]);
  }
  setup_kernel(p, nfreq(), omega_coord_);
  vector_type S(ndat());
  matrix_type Kt = K_; // gesvd destroys K!
  bindings::lapack::gesvd(Kt, S, U_, Vt_);
  std::cout << "# Singular values of the Kernel:\n";
  const double prec = sqrt(std::numeric_limits<double>::epsilon())*nfreq()*S[0];
  for (int s=0; s<S.size(); ++s) 
    ns_ = (S[s] >= prec) ? s : ns_;
  if (ns() == 0) 
    boost::throw_exception(std::logic_error("all singular values smaller than the precision"));
  U_.resize(ndat(), ns_, true);
  Vt_.resize(ns_, nfreq(), true);
  Sigma_.resize(ns_, ns_);
  Sigma_.clear();
  for (int s=0; s<ns_; ++s) {
    std::cout << "# " << s << "\t" << S[s] <<"\n";
    Sigma_(s,s) = S[s];
  }
  matrix_type Ut = ublas::trans(U_);
  vector_type t = ublas::prec_prod(Ut, y_);
  vector_type y2 = ublas::prec_prod(U_, t);
  double chi = ublas::norm_2(y_-y2);
  std::cout << "minimal chi2: " << chi*chi/y_.size() << std::endl;
}




