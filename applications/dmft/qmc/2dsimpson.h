 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *               2012        by Jakub Imriska <jimriska@phys.ethz.ch>
 *
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
#ifndef INTEGRANDH
#define INTEGRANDH
#include<complex>
#include<vector>
#include<stdexcept>
#include<iostream>

#include <alps/parameter.h>


class BetheBandstructure {
  /* does load the hopping amplitudes 't' for different flavors;
     4t=2W=D, where D is the bandwidth;
     assumes, that the flavors (2m) and (2m+1) have the same bandwidth;
     thus 't0' sets the bandwidth for flavors 0 and 1
          't1' sets the bandwidth for flavors 2 and 3, etc.;
     second moment of the semicicular DOS is 't^2'
  */
  public:
    BetheBandstructure(const alps::Parameters& parms, bool verbose=false) {
      unsigned int n_flavor = parms.value_or_default("FLAVORS", 2);
      tsq_.resize(n_flavor/2);
      for(unsigned int i=0; i<n_flavor/2; ++i){
        std::stringstream t_i; t_i<<"t"<<i;  // flavors (2m) and (2m+1) assumed to have the same bandwidth
        double t = (parms.defined(t_i.str()) ? static_cast<double>(parms[t_i.str()]) : static_cast<double>(parms["t"]));
        tsq_[i]=t*t;   // second moment of the Bethe DOS is t^2 [4t=2W=D]
        if (verbose)
          std::cout<<"Semicircular DOS: for flavors "<<2*i<<" and "<<2*i+1<<" using hopping t = "<<t<<std::endl;
      }
    }
    double tsq(unsigned int f=0) const { return tsq_[f/2]; }
  private:
    std::vector<double> tsq_;
};


class TwoDBandstructure {
  /* This class contains all the non-interacting 2-dimensional lattice information for DMFT with option TWODBS, which 
     specifies that the Hilbert transformation is performed by integration over k-space.
     
     Currently, we support:
       square lattice with nearest-neighbor and next nearest-neighbor hoppings
       hexagonal lattice with nearest-neighbor hoppings
     
     Square lattice:
       BZ is chosen to be:  <-pi,pi)x<-pi,pi)
       
     Hexagonal lattice:
       Here we work in reduced coordinates:  kx_tilde = 3/2*kx;  ky_tilde = sqrt(3)/2*ky
       BZ is chosen rectangular: kx_tilde in <-pi,pi);  ky in <-pi/2,pi/2)   */
       
  public:
  enum lattice_t {square, hexagonal};
  
  // constructor
  TwoDBandstructure(const alps::Parameters & parms)
    : t_(parms.value_or_default("t",1.)),
      tprime_(parms.value_or_default("tprime",0))
  {
    if (parms.defined("TWODBS") && parms["TWODBS"]=="hexagonal") {
      lattice_ = hexagonal;
      if (tprime_!=0) 
        throw std::runtime_error("Unsupported option! For TWODBS==hexagonal only nearest-neighbor hoppings are implemented.");
    }
    else
      lattice_ = square;
    std::cout << "TwoDBandstructure: " << (lattice_==hexagonal ? "hexagonal lattice" : "square lattice") << "; n.n. hopping : " << t_ << "; n.n.n. hopping : " << tprime_ << (parms.defined("tprime") ? "" : " (unset)") << std::endl;
  }
  
  double dispersion(double kx, double ky) const {
    if (lattice_==hexagonal) {
      double c = cos(ky);
      return t_*sqrt(1. + 4.*c*(c+cos(kx)));
    }
    else { //square
      return -2.*t_*(cos(kx)+cos(ky)) -4.*tprime_*cos(kx)*cos(ky);
    }
  }
  bool single_band() const { return (lattice_==hexagonal ? false : true); }
  double k_area() const { return (kx_max()-kx_min())*(ky_max()-ky_min()); }
  double kx_min() const { return -M_PI; }
  double ky_min() const { return (lattice_==hexagonal ? -M_PI/2. : -M_PI); }
  double kx_max() const { return M_PI; }
  double ky_max() const { return (lattice_==hexagonal ? M_PI/2. : M_PI); }
  double second_moment() const { return (lattice_==hexagonal ? 3.*t_*t_ : 4.*(t_*t_ + tprime_*tprime_) ); }
  
  private:
  lattice_t lattice_;
  const double t_, tprime_;
};

double dispersion(double kx, double ky, double t, double tprime);
double dispersion_hexagonal(double kx, double ky, double t);

class DOS_integrand {
  public:
  // AFM version:
  DOS_integrand(const std::vector<double> &dos, const std::vector<double> &eps, 
                const std::complex<double> &zeta_A, const std::complex<double> &zeta_B)
    : dos_(dos), eps_(eps), zeta_A_(zeta_A), zeta_B_(zeta_B), mode_(-2)
    {}
  // PM version:
  DOS_integrand(const std::vector<double> &dos, const std::vector<double> &eps, 
                const std::complex<double> &zeta_A)
    : dos_(dos), eps_(eps), zeta_A_(zeta_A), mode_(-1)
    {}
  // version to calculate moments of DOS:
  DOS_integrand(const unsigned int & power, const std::vector<double> &dos, const std::vector<double> &eps) 
    : dos_(dos), eps_(eps), mode_(power)
    { if (mode_>2) throw std::logic_error("DOS_power_of_eps: unsupported power"); }
    
  std::complex<double> operator()(unsigned int n) const {
    if (mode_<0)
      return (mode_==-1 ? dos_[n]/(zeta_A_-eps_[n]) : dos_[n]/(zeta_A_*zeta_B_ - eps_[n]*eps_[n]) );
    else
      return dos_[n] * (mode_==0 ? 1. : (mode_==1 ? eps_[n] : eps_[n]*eps_[n]));
  }
    unsigned int size() const { return dos_.size(); }

  private:
    const std::vector<double> & dos_, eps_;
    const std::complex<double> zeta_A_, zeta_B_;
    int mode_;
};

class integrand{
  public:
    virtual std::complex<double>operator()(double kx, double ky) const=0;
    virtual ~integrand(){};
};

class pm_band_integrand:public integrand{
  public:
    pm_band_integrand(const TwoDBandstructure & bandstruct, std::complex<double> zeta_A, double band=1.)
      : bandstruct_(bandstruct), zeta_A_(zeta_A), band_(band)
      {}
    std::complex<double>operator()(double kx, double ky) const
      { return 1./(zeta_A_-band_*bandstruct_.dispersion(kx,ky)); }
  protected:
    const TwoDBandstructure & bandstruct_;
    const std::complex<double> zeta_A_;
    const double band_;
};

class afm_band_integrand:public pm_band_integrand{
  public:
    afm_band_integrand(const TwoDBandstructure & bandstruct, std::complex<double> zeta_A, std::complex<double> zeta_B)
      : pm_band_integrand(bandstruct, zeta_A), zeta_B_(zeta_B)
      {}
    std::complex<double>operator()(double kx, double ky) const
      { return 1./(zeta_A_*zeta_B_-bandstruct_.dispersion(kx,ky)*bandstruct_.dispersion(kx,ky)); }
  private:
    const std::complex<double> zeta_B_;
};

class pm_integrand:public integrand{
  public:
    pm_integrand(double t, double tprime, std::complex<double> zeta_A)
      : t_(t), tprime_(tprime), zeta_A_(zeta_A)
      {}
    std::complex<double> operator()(double kx, double ky) const
      { return 1./(zeta_A_-dispersion(kx,ky,t_,tprime_ )); }
  protected:
    const double t_, tprime_;
    const std::complex<double> zeta_A_;
};

class afm_integrand:public pm_integrand{ //class for the integrand with an epsilon(k) that uses the t-band structure.
  public:
    afm_integrand(double t, double tprime, std::complex<double> zeta_A, std::complex<double> zeta_B)
      : pm_integrand(t,tprime,zeta_A), zeta_B_(zeta_B)
      {}
    std::complex<double> operator()(double kx, double ky) const 
      { return 1./(zeta_A_*zeta_B_-dispersion(kx,ky,t_,tprime_)*dispersion(kx,ky,t_,tprime_)); }
  private:
    const std::complex<double> zeta_B_;
};

std::complex<double> simpson_integrate(const DOS_integrand &);
std::complex<double> twodsimpson(const integrand &, double ax, double ay, double bx, double by, int N);
#endif
