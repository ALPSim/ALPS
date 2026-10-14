/*
 * Copyright (C) 1998-2018 ALPS Collaboration.
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

// Legacy ALPS time and frequency kernels from src/apps/maxent/src/maxent_parms.cpp
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

// Adapted from CQMP/Maxent src/maxent_kernel.cpp, revision
// a3992534228b0de732e033ae23154969724814be, on 2026-10-14.
// Uses ALPS parameter APIs and retains ALPS input/output conventions;
// see ../README.md for provenance and licensing.
#include "maxent_kernel.hpp"
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/legendre.hpp>    //for Legendre kernel
#include <boost/math/special_functions/factorials.hpp>  //for Legendre kernel
#include <boost/math/quadrature/gauss_kronrod.hpp>

namespace bmth = boost::math;

kernel::kernel(alps::params &p, const vector_type& freq, vector_type &inputGrid):
ndat_(p["NDAT"]),
nfreq_(p["NFREQ"]),
T_(p["T"] | (1./static_cast<double>(p["BETA"]))),
K_(ndat_,nfreq_)
{
  //K_.clear();
  K_=matrix_type::Zero(ndat_,nfreq_);
  std::string dataspace_name = p["DATASPACE"];
  std::string kernel_name = p["KERNEL"];
  boost::to_lower(dataspace_name);
  boost::to_lower(kernel_name);
  bool ph_symmetry=p["PARTICLE_HOLE_SYMMETRY"];
  std::cout<<"using kernel "<<kernel_name<<" in domain "<<dataspace_name;
  if(ph_symmetry) std::cout<<" with ph symmetry"; else std::cout<<" without ph symmetry"; std::cout<<std::endl;

  set_kernel_type(dataspace_name,kernel_name, ph_symmetry);

  if (dataspace_name == "time") {
    tau_points_.resize(ndat_);
    for (int i=0; i<ndat_; ++i) {
      const std::string key = "TAU_" + boost::lexical_cast<std::string>(i);
      // Retain ALPS's per-point overrides and uniform default for all inputs.
      tau_points_[i] = p[key] | (i / ((ndat_-1)*T_));
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
        K_(i,0) = T_;
        for (int j=0; j<nfreq_; ++j) {
          double omega = freq[j];
          K_(i,j) = 0.5*omega * (std::exp(-omega*tau) + std::exp(-omega*(1./T_-tau))) / (1 - std::exp(-omega/T_));
        }
      }
    }
    else if(ktype_== time_fermionic_legendre_kernel || ktype_==time_bosonic_legendre_kernel)
        setup_legendre_kernel(p,freq,ndat_);
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
    else if(ktype_==time_fermionic_kernel){
        for (int i=0; i<ndat_; ++i) {
						double tau=tau_points_[i];
            for (int j=0; j<nfreq_; ++j) {
                double omega = freq[j];
                K_(i,j) =  -1.;
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
            if(dtype_==legendre_dataspace)
                ktype_=time_fermionic_legendre_kernel;
            else
                ktype_=time_fermionic_kernel;
    else if(kernel_name=="bosonic")
            if(dtype_==legendre_dataspace)
                ktype_=time_bosonic_legendre_kernel;
            else
                ktype_=time_bosonic_kernel;
    else if(kernel_name=="tzero" || kernel_name=="boris")
      ktype_=time_boris_kernel;
    else throw std::invalid_argument("unknown kernel name. In the time domain it should be fermionic, bosonic, or tzero.");
  }else if(dtype_ == legendre_dataspace){
      if(kernel_name=="fermionic")
          ktype_=legendre_fermionic_kernel;
      else if(kernel_name=="bosonic")
          // The imported quadrature implements only the fermionic kernel.
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
void kernel::setup_legendre_kernel(const alps::params &, const vector_type& freq, const int lmax) {
  if (lmax > boost::math::max_factorial<double>::value)
    throw std::runtime_error("lmax is greater than boost factorial precision");

  const double beta = 1. / T_;
  for (int l = 0; l < lmax; ++l) {
    for (int j = 0; j < nfreq_; ++j) {
      const double z = beta * freq[j];
      // Adapt CQMP's GSL Gauss-61 integration to Boost.Math's adaptive
      // Gauss-Kronrod rule. Dimensionless time and nonpositive exponents
      // preserve the fermionic kernel without overflow for negative beta*omega.
      if (z == 0.) {
        K_(l,j) = l == 0 ? -beta / 2. : 0.;
        continue;
      }
      const auto integrand = [l, z](double t) {
        const double weight = z >= 0.
            ? std::exp(-t*z) / (1. + std::exp(-z))
            : std::exp((1.-t)*z) / (1. + std::exp(z));
        return bmth::legendre_p(l, 2.*t-1.) * weight;
      };
      const double integral = bmth::quadrature::gauss_kronrod<double, 61>::integrate(
          integrand, 0., 1., 10, 1.49e-8);
      K_(l,j) = -std::sqrt(2.*l+1.) * beta * integral;
    }
  }
}
