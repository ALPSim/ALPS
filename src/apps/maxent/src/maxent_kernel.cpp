/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2026 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Legacy ALPS time and frequency kernels from src/apps/maxent/src/maxent_parms.cpp
// at d4533fd30c65656cf5ab7fa448dd6b921b39aa6f.

// Adapted from CQMP/Maxent src/maxent_kernel.cpp, revision
// e727650a1f0a427e95efebc409d907b4c84011ec, on 2026-10-14.
// Uses ALPS's HDF5 archive and retains ALPS input/output conventions.
#include "maxent_kernel.hpp"
#include <cmath>
#include "maxent_string.hpp"


kernel::kernel(alps::params &p, const vector_type& freq, vector_type &inputGrid):
ndat_(p["NDAT"]),
nfreq_(p["NFREQ"]),
T_(1./static_cast<double>(p["BETA"])),
K_(ndat_,nfreq_)
{
  //K_.clear();
  K_=matrix_type::Zero(ndat_,nfreq_);
  std::string dataspace_name = p["DATASPACE"];
  std::string kernel_name = p["KERNEL"];
  to_lower(dataspace_name);
  to_lower(kernel_name);
  bool ph_symmetry=p["PARTICLE_HOLE_SYMMETRY"];
  std::cout<<"using kernel "<<kernel_name<<" in domain "<<dataspace_name;
  if(ph_symmetry) std::cout<<" with ph symmetry"; else std::cout<<" without ph symmetry"; std::cout<<std::endl;

  set_kernel_type(dataspace_name,kernel_name, ph_symmetry);

  if (dataspace_name == "time") {
    tau_points_.resize(ndat_);
    for (int i=0; i<ndat_; ++i) {
      const std::string key = "TAU_" + std::to_string(i);
      // ALPS: TAU_i override a uniform grid, point by point, for every input.
      if (!p.defined(key) && p.supplied(key)) p.define<double>(key, "");
      tau_points_[i] = p.exists(key) ? p[key].as<double>() : i / ((ndat_-1)*T_);
    }
    inputGrid = tau_points_;
  }
    
  if(ktype_==legendre_fermionic_kernel || ktype_==legendre_bosonic_kernel){
      setup_legendre_kernel(p,freq,ndat_);
  }else if(ktype_==time_fermionic_kernel){
      for (int i=0; i<ndat_; ++i) {
        double tau=tau_points_[i];
        for (int j=0; j<nfreq_; ++j) {
          double omega = freq[j]; //Default().omega_of_t(double(j)/(nfreq_-1));
          K_(i,j) =  -1. / (std::exp(omega*tau) + std::exp(-omega*(1./T_-tau)));
        }
      }
    }
    else if (ktype_==time_bosonic_kernel) {
      for (int i=0; i<ndat_; ++i) {
        double tau=tau_points_[i];
        for (int j=0; j<nfreq_; ++j) {
          // ALPS adaptation: the kernel is even in omega, so evaluate it at
          // |omega|, where no exponential overflows. Upstream's form overflows
          // for beta*omega < -709 and maps the NaN to 0, which is wrong near
          // tau = 0 and beta: there the kernel approaches |omega|/2.
          double omega = std::abs(freq[j]);
          K_(i,j) = omega == 0.0 ? T_ :
            0.5*omega * (std::exp(-omega*tau) + std::exp(-omega*(1./T_-tau))) / -std::expm1(-omega/T_);
        }
      }
    }
    //for zero temperature, only positive frequency matters
    else if (ktype_ == time_boris_kernel) {
      for (int i=0; i<ndat_; ++i) {
        double tau=tau_points_[i];
        for (int j=0; j<nfreq_; ++j) {
          double omega = freq[j];
          K_(i,j) = -std::exp(-omega*tau);
        }
      }
    }
    else if(ktype_==frequency_fermionic_ph_kernel) {
    for (int i=0; i<ndat_; ++i) {
      double omegan = (2*i+1)*M_PI*T_;
      inputGrid(i) = omegan;
      for (int j=0; j<nfreq_; ++j) {
        double omega = freq[j];
        K_(i,j) =  -omegan / (omegan*omegan + omega*omega);
      }
    }
  }
  else if (ktype_==frequency_bosonic_ph_kernel) {
    for (int i=0; i<ndat_; ++i) {
      double Omegan = (2*i)*M_PI*T_;
      inputGrid(i) = Omegan;
      for (int j=0; j<nfreq_; ++j) {
        double Omega = freq[j];
        if(Omega ==0) throw std::runtime_error("Bosonic kernel is singular at frequency zero. Please use grid w/o evaluation at zero.");
        // Preserve the sign convention of ALPS bosonic frequency data.
        K_(i,j) =  -Omega*Omega / (Omegan*Omegan + Omega*Omega);
      }
    }
  }else if (ktype_==frequency_anomalous_ph_kernel) {
    for(int i=0;i<ndat_;++i){
      double omegan = (2*i+1)*M_PI*T_;
      inputGrid(i) = omegan;
      for (int j=0; j<nfreq_; ++j) {
        double omega = freq[j];
        K_(i,j) =  omega*omega / (omegan*omegan + omega*omega);
      }
    }
  }
  else{
    complex_matrix_type Kc(ndat_/2, nfreq_);
    if (ktype_==frequency_fermionic_kernel) {
      ///ndat/2 is defined as such below because ndat=number of points inputed
      ///if ph symmetry, then ndat=number of imag points
      ///otherwise ndat=total number of real+imag points = ndat/2 data points
      for (int i=0; i<ndat_/2; ++i) {
        std::complex<double> iomegan(0, (2*i+1)*M_PI*T_);
        inputGrid(i) = iomegan.imag();
        for (int j=0; j<nfreq_; ++j) {
          double omega = freq[j];
          Kc(i,j) =  1. / (iomegan - omega);
        }
      }
    }
    else if (ktype_==frequency_bosonic_kernel){
      for (int i=0; i<ndat_/2; ++i) {
        std::complex<double> iomegan(0, 2*i*M_PI*T_);
        inputGrid(i) = iomegan.imag();
        inputGrid(i+1) = iomegan.imag();
        for (int j=0; j<nfreq_; ++j) {
          double omega = freq[j];
          // Retain the ALPS complex bosonic convention as well as its PH limit.
          Kc(i,j) =  omega / (iomegan - omega);
        }
      }
    }
    else if (ktype_==frequency_anomalous_kernel){
      for (int i=0; i<ndat_/2; ++i) {
        std::complex<double> iomegan(0, (2*i+1)*M_PI*T_);
        inputGrid(i) = iomegan.imag();
        for (int j=0; j<nfreq_; ++j) {
          double omega = freq[j];
          Kc(i,j) =  -omega / (iomegan - omega);
        }
      }
    }else
      throw std::logic_error("unknown kernel type");

    for (int i=0; i<ndat_; i+=2) {
      //TODO: understand the j=1 here
      for (int j=0; j<nfreq_; ++j) {
        K_(i,j) = Kc(i/2,j).real();
        K_(i+1,j) = Kc(i/2,j).imag();
      }
    }
  }
}

void kernel::set_kernel_type(const std::string &dataspace_name, const std::string &kernel_name,
                             bool ph_symmetry){
  if(dataspace_name=="time"){
    dtype_=time_dataspace;
     
  }else if(dataspace_name=="frequency"){
    dtype_=frequency_dataspace;
  }else if(dataspace_name=="legendre"){
    dtype_=legendre_dataspace;
  }
  else
    throw std::invalid_argument("unknown dataspace name. it should be time, frequency, or legendre");

  if(dtype_==time_dataspace){
    if(kernel_name=="fermionic")
      ktype_=time_fermionic_kernel;
    else if(kernel_name=="bosonic")
      ktype_=time_bosonic_kernel;
    else if(kernel_name=="tzero" || kernel_name=="boris")
      ktype_=time_boris_kernel;
    else throw std::invalid_argument("unknown kernel name. In the time domain it should be fermionic, bosonic, or tzero.");
  }else if(dtype_ == legendre_dataspace){
      if(kernel_name=="fermionic")
          ktype_=legendre_fermionic_kernel;
      else if(kernel_name=="bosonic")
          // The imported closed form implements only the fermionic kernel.
          throw std::invalid_argument("bosonic Legendre input is not supported");
      else throw std::invalid_argument("unknown kernel name. In the legendre domain it should be fermionic or bosonic");
  }
  else{
    if(ph_symmetry){
      if(kernel_name== "fermionic")
        ktype_=frequency_fermionic_ph_kernel;
      else if (kernel_name=="bosonic")
        ktype_=frequency_bosonic_ph_kernel;
      else if (kernel_name=="anomalous")
        ktype_=frequency_anomalous_ph_kernel;
      else throw std::invalid_argument("unknown kernel name. In the particle hole symmetric frequency domain it should be fermionic, bosonic, or anomalous.");
    }else{
      if(kernel_name== "fermionic")
        ktype_=frequency_fermionic_kernel;
      else if (kernel_name=="bosonic")
        ktype_=frequency_bosonic_kernel;
      else if (kernel_name=="anomalous")
        ktype_=frequency_anomalous_kernel;
      else throw std::invalid_argument("unknown kernel name. In the non-particle hole symmetric frequency domain it should be fermionic, bosonic, or anomalous.");
    }
  }

}
namespace {

/// Scaled modified spherical Bessel functions of the first kind,
/// s[l] = i_l(x) exp(-x) for l = 0..lmax-1 and x >= 0.
/// The ratios f_l = i_l/i_{l-1} satisfy f_l = x/(2l+1 + x f_{l+1}); this
/// continued fraction is evaluated downward from an order N that is doubled
/// until f_1..f_{lmax-1} converge, and i_0(x) exp(-x) = (1 - exp(-2x))/(2x).
/// Only products and sums of positive numbers occur, so there is no
/// cancellation for any x.
std::vector<double> scaled_spherical_bessel_i(int lmax, double x) {
  std::vector<double> s(lmax, 0.);
  if (lmax == 0) return s;
  s[0] = x == 0. ? 1. : -std::expm1(-2. * x) / (2. * x);
  if (x == 0. || lmax == 1) return s;
  std::vector<double> f(lmax, 0.), f_prev;
  for (int n = lmax + 20 + 2 * static_cast<int>(std::ceil(x));; n *= 2) {
    double fl = 0.;
    for (int l = n; l >= 1; --l) {
      fl = x / (2 * l + 1 + x * fl);
      if (l < lmax) f[l] = fl;
    }
    bool converged = !f_prev.empty();
    for (int l = 1; converged && l < lmax; ++l)
      converged = std::abs(f[l] - f_prev[l]) <= 1e-16 * std::abs(f[l]);
    if (converged) break;
    f_prev = f;
  }
  for (int l = 1; l < lmax; ++l) s[l] = s[l - 1] * f[l];
  return s;
}

}  // namespace

/// Legendre kernel
///   K(l, omega) = -sqrt(2l+1) int_0^beta P_l(2 tau/beta - 1) exp(-tau omega) / (1 + exp(-beta omega)) dtau
/// in closed form: with a = beta omega / 2,
///   int_0^beta P_l(2 tau/beta - 1) exp(-tau omega) dtau = beta exp(-a) (-1)^l i_l(a),
/// so K(l, omega) = -sqrt(2l+1) beta (-1)^l i_l(a) / (2 cosh a). With i_l(-x) = (-1)^l i_l(x)
/// this is evaluated as -sqrt(2l+1) beta sigma_l [i_l(|a|) exp(-|a|)] / (1 + exp(-2|a|)),
/// sigma_l = (-1)^l for a > 0 and 1 otherwise, which stays finite for any beta omega.
/// ALPS rejects bosonic Legendre input: upstream applies this fermionic kernel
/// to bosonic data (its open item B17).
void kernel::setup_legendre_kernel(const alps::params &/*p*/, const vector_type& freq, const int lmax){
    const double beta = 1. / T_;
    for (int j = 0; j < nfreq_; ++j) {
        const double a = 0.5 * beta * freq[j];
        const double x = std::abs(a);
        const std::vector<double> s = scaled_spherical_bessel_i(lmax, x);
        const double denominator = 1. + std::exp(-2. * x);
        for (int l = 0; l < lmax; ++l) {
            const double sign = (a > 0. && l % 2 == 1) ? -1. : 1.;
            K_(l, j) = -std::sqrt(2. * l + 1.) * beta * sign * s[l] / denominator;
        }
    }
}
