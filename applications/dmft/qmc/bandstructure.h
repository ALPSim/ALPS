 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2013        by Jakub Imriska <jimriska@phys.ethz.ch>
 *
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/
 
#ifndef BANDSTRUCTUREH__
#define BANDSTRUCTUREH__

#include<complex>
#include<vector>
#include<stdexcept>
#include<iostream>

#include <alps/parameter.h>

class Bandstructure {
  /*
   * class providing interface for a universal bandstructure
   * for 1 site only
   */
public:
  
  void set_parms(alps::Parameters& parms) const;
  
  Bandstructure(const alps::Parameters& parms)
    : epssq_(static_cast<unsigned>(parms.value_or_default("FLAVORS", 2)),-1.),
      eps_(static_cast<unsigned>(parms.value_or_default("FLAVORS", 2)),0.)
    {}
  
  virtual std::complex<double> HilbertIntegral_PM(const std::complex<double>& zeta_A, const unsigned flavor) const;
  virtual std::complex<double> HilbertIntegral_AFM(const std::complex<double>& zeta_A_times_zeta_B, const unsigned flavor) const;
  
  double first_moment(int flavor) const { return eps_[flavor]; }
  double second_moment(int flavor) const { return epssq_[flavor]; }
  
  virtual ~Bandstructure() {}
  
protected:
  std::vector<double> epssq_;
  std::vector<double> eps_;
};

// ---------------------------------------------------------------------------------------------

class DOSBandstructure : public Bandstructure {
  /*
   * class for bandstructure set by density of states (DOS)
   * DOS is loaded from a textfile
   * flavors 2n and 2n+1 share the same DOS
   */
public:

  DOSBandstructure(const alps::Parameters& parms, bool verbose=false);
  
  std::complex<double> HilbertIntegral_PM(const std::complex<double>& zeta_A, const unsigned flavor) const;
  std::complex<double> HilbertIntegral_AFM(const std::complex<double>& zeta_A_times_zeta_B, const unsigned flavor) const;
  
private:
  std::vector<std::vector<double> > dos_;
  std::vector<double> dos_min_, dos_step_;
};

// ---------------------------------------------------------------------------------------------

class SemicircleBandstructure : public Bandstructure {
  /*
   * for semicircular DOS
   * does load the hopping amplitudes 't' for different flavors;
     4t=2W=D, where D is the bandwidth;
     assumes, that the flavors (2m) and (2m+1) have the same bandwidth;
     thus 't0' sets the bandwidth for flavors 0 and 1
          't1' sets the bandwidth for flavors 2 and 3, etc.;
     second moment of the semicicular DOS is 't^2'
     first moment is 0
   */
public:
  
  SemicircleBandstructure(const alps::Parameters& parms, bool verbose=false)
    : Bandstructure(parms)
  {
    unsigned n_flavor=parms.value_or_default("FLAVORS",2);
    if (n_flavor%2!=0)
      throw std::logic_error("SemicircleBandstructure: current implementation does not allow for odd n_flavor.");
    for(unsigned int f=0; f<n_flavor/2; ++f){
      std::stringstream t_f; t_f<<"t"<<f;  // flavors (2m) and (2m+1) assumed to have the same bandwidth
      double t = (parms.defined(t_f.str()) ? static_cast<double>(parms[t_f.str()]) : static_cast<double>(parms["t"]));
      epssq_[2*f]=t*t;   // second moment of the Bethe DOS is t^2 [4t=2W=D]
      epssq_[2*f+1]=t*t;
      if (verbose)
        std::cout<<"BANDSTRUCTURE:"<<std::endl<<"Semicircular DOS: for flavors "<<2*f<<" and "<<2*f+1<<" using hopping t ="<<t<<std::endl;
    }
  }

  std::complex<double> HilbertIntegral_PM(const std::complex<double>& zeta_A, const unsigned flavor) const;
  std::complex<double> HilbertIntegral_AFM(const std::complex<double>& zeta_A_times_zeta_B, const unsigned flavor) const;
  
};

// ---------------------------------------------------------------------------------------------

class SquareLatticeBandstructure : public Bandstructure {
  /* This class contains all the non-interacting 2-dimensional lattice information for DMFT
   * for square lattice with nearest-neighbor and next nearest-neighbor hoppings
   * it is used with option TWODBS, which specifies that the Hilbert transformation is performed by integration over k-space.
   * 
   * BZ is chosen to be:  <-pi,pi)x<-pi,pi)
   */
  
public:
  
  // constructor
  SquareLatticeBandstructure(const alps::Parameters & parms, bool verbose=false)
    : Bandstructure(parms),
      t_(static_cast<double>(parms.value_or_default("t",1.))),
      tprime_(static_cast<double>(parms.value_or_default("tprime",0.))),
      L_(static_cast<int>(parms.value_or_default("L",128)))
  {
    if (static_cast<int>(parms.value_or_default("FLAVORS",2))!=2)
      throw std::logic_error("TwoDBandstructure: current implementation does not allow for n_flavor!=2.");
    epssq_[0]=4.*(t_*t_ + tprime_*tprime_);
    epssq_[1]=epssq_[0];
    if (verbose)
      std::cout<<"BANDSTRUCTURE:"<<std::endl<<"Square lattice: n.n. hopping = " << t_ << "; n.n.n. hopping = " << std::endl;
  }
  
  double operator()(double kx, double ky) const {
    return -2.*t_*(cos(kx)+cos(ky)) -4.*tprime_*cos(kx)*cos(ky);
  }
  
  std::complex<double> HilbertIntegral_PM(const std::complex<double>& zeta_A, const unsigned flavor) const;
  std::complex<double> HilbertIntegral_AFM(const std::complex<double>& zeta_A_times_zeta_B, const unsigned flavor) const;
  
private:
  static const double kx_min;
  static const double ky_min;
  static const double kx_max;
  static const double ky_max;
  static double k_area() { return (kx_max-kx_min)*(ky_max-ky_min); }
  
  const double t_, tprime_;
  const int L_;
};

// ---------------------------------------------------------------------------------------------

class HexagonalLatticeBandstructure : public Bandstructure {
  /* This class contains all the non-interacting 2-dimensional lattice information for DMFT
   * for hexagonal lattice with nearest-neighbor hoppings
   * it is used with option TWODBS, which specifies that the Hilbert transformation is performed by integration over k-space
   *
   * Here we work in reduced coordinates:  kx_tilde = 3/2*kx;  ky_tilde = sqrt(3)/2*ky
   * BZ is chosen rectangular: kx_tilde in <-pi,pi);  ky in <-pi/2,pi/2)   
   */
  
public:
  
  // constructor
  HexagonalLatticeBandstructure(const alps::Parameters & parms, bool verbose=false)
    : Bandstructure(parms),
      t_(static_cast<double>(parms.value_or_default("t",1.))),
      L_(static_cast<int>(parms.value_or_default("L",128)))
  {
    if (static_cast<double>(parms.value_or_default("tprime",0.))!=0.)
      std::cout<<"WARNING: for the hexagonal lattice the parameter 'tprime' is ignored."<<std::endl;
    if (static_cast<int>(parms.value_or_default("FLAVORS",2))!=2)
      throw std::logic_error("TwoDBandstructure: current implementation does not allow for n_flavor!=2.");
    epssq_[0]=3.*t_*t_;
    epssq_[1]=epssq_[0];
    if (verbose)
      std::cout<<"BANDSTRUCTURE:"<<std::endl<<"Hexagonal lattice: n.n. hopping = " << t_ << std::endl;
  }
  
  double operator()(double kx, double ky) const {
    /* Note, that this is a 2-band model.
     * However, the dispersion is simply
     *   epsilon_{n=1,kx,ky}=-epsilon_{n=2,kx,ky}
     * thus we do not essentially need the band index as a parameter and we return just the positive value
     */
    double c = cos(ky);
    return t_*sqrt(1. + 4.*c*(c+cos(kx)));
  }
  
  std::complex<double> HilbertIntegral_PM(const std::complex<double>& zeta_A, const unsigned flavor) const;
  std::complex<double> HilbertIntegral_AFM(const std::complex<double>& zeta_A_times_zeta_B, const unsigned flavor) const;
  
private:
  static const double kx_min;
  static const double ky_min;
  static const double kx_max;
  static const double ky_max;
  static double k_area() { return (kx_max-kx_min)*(ky_max-ky_min); }
  
  const double t_;
  const int L_;
};

// ---------------------------------------------------------------------------------------------

boost::shared_ptr<Bandstructure> BandstructureFactory(const alps::Parameters& parms, bool verbose=false);

#endif //BANDSTRUCTUREH__
