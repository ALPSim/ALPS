 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *               2012 - 2013 by Jakub Imriska <jimriska@phys.ethz.ch>
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
 
#include<iostream>
#include "bandstructure.h"


// implementation of 2d simpson
// for integration over rectangular area
// using 2N+1 mesh points in each direction
template<class integrand>
std::complex<double> twodsimpson(const integrand& f, double ax, double ay, double bx, double by, int N){
  double h=(bx-ax)/(2.*N);
  double k=(by-ay)/(2.*N);

  std::complex<double> result = f(ax,ay) + f(ax,by) + f(bx, ay) + f(bx,by);
  
  // values between boundaries
  for ( int i = 1; i <= N; ++i ) {
    result += 4.*f(ax,ay+k*(2*i-1));
  }
  for ( int i = 1; i <= N-1; ++i ) {
    result += 2.*f(ax,ay+k*2*i);
  }
  for ( int i = 1; i <= N; ++i ) {
    result += 4.*f(bx,ay+k*(2*i-1));
  }
  for ( int i = 1; i <= N-1; ++i ) {
    result += 2.*f(bx,ay+k*2*i);
  }
  
  
  for ( int i = 1; i <= N; ++i ) {
    result += 4.*f(ax+h*(2*i-1),ay);
  }
  for ( int i = 1; i <= N-1; ++i ) {
    result += 2.*f(ax+h*(2*i),ay);
  }
  for ( int i = 1; i <= N; ++i ) {
    result += 4.*f(ax+h*(2*i-1),by);
  }
  for ( int i = 1; i <= N-1; ++i ) {
    result += 2.*f(ax+h*(2*i),by);
  }

  //inner part
  //
  for(int i=1;i<=N;++i){
    for(int j=1;j<=N;++j){
      result+=16.*f(ax+h*(2*i-1), ay+k*(2*j-1));
    }
  }
  for(int i=1;i<=N;++i){
    for(int j=1;j<N;++j){
      result+=8.*f(ax+h*(2*i-1), ay+k*(2*j));
    }
  }
  for(int i=1;i<N;++i){
    for(int j=1;j<=N;++j){
      result+=8.*f(ax+h*(2*i), ay+k*(2*j-1));
    }
  }
  for(int i=1;i<N;++i){
    for(int j=1;j<N;++j){
      result+=4.*f(ax+h*(2*i), ay+k*(2*j));
    }
  }
  
  result *= h*k/9.;
  return result;
}


// 1D Simpson integration
template <class Integrand>
typename Integrand::result_type simpson_integrate(const Integrand& integrand) {
  if (integrand.size() % 2 != 1) { throw std::runtime_error("simpson_integrate: ERROR: use odd number of DOS points for higher precision (due to Simpson integration)"); }

  typename Integrand::result_type sum1=0, sum2=0;
  // Simpson integration
  for(unsigned i=1;i<integrand.size()-2;i+=2){
    sum1+=integrand(i);
    sum2+=integrand(i+1);                   // Assuming equidistant epsilon points
  }
  return (4.*sum1+2.*sum2+integrand(0)+integrand(integrand.size()-1)+4.*integrand(integrand.size()-2))/3.;
}


// ---------------------------------------------------------------------------------------------


std::complex<double> Bandstructure::HilbertIntegral_PM(const std::complex<double>& zeta_A, const unsigned flavor) const {
  throw std::logic_error("This should not be executed at any time.");
}

std::complex<double> Bandstructure::HilbertIntegral_AFM(const std::complex<double>& zeta_A_times_zeta_B, const unsigned flavor) const {
  throw std::logic_error("This should not be executed at any time.");
}

void Bandstructure::set_parms(alps::Parameters& parms) const {
  for (unsigned int f=0; f<epssq_.size(); ++f) {
    if (parms.defined("EPS_"+boost::lexical_cast<std::string>(f))) {
      double value_by_user=static_cast<double>(parms["EPS_"+boost::lexical_cast<std::string>(f)]);
      if (std::abs(value_by_user-eps_[f])/(std::abs(value_by_user)+std::abs(eps_[f]))>1e-3 && std::abs(value_by_user-eps_[f])>1e-10) {
        std::cout<<"WARNING: Bandstructure: inconsistency in parameter 'EPS_"<<f<<"': value set by user = "<<value_by_user<<"; value expected/computed by program = "<<eps_[f]<<std::endl;
      }
    } else {
      parms["EPS_"+boost::lexical_cast<std::string>(f)] = eps_[f];
    }
    if (parms.defined("EPSSQ_"+boost::lexical_cast<std::string>(f))) {
      double value_by_user=static_cast<double>(parms["EPSSQ_"+boost::lexical_cast<std::string>(f)]);
      if (std::abs(value_by_user-epssq_[f])/(std::abs(value_by_user)+std::abs(epssq_[f]))>1e-3 && std::abs(value_by_user-epssq_[f])>1e-10) {
        std::cout<<"WARNING: Bandstructure: inconsistency in parameter 'EPSSQ_"<<f<<"': value set by user = "<<value_by_user<<"; value expected/computed by program = "<<epssq_[f]<<std::endl;
      }
    } else {
      parms["EPSSQ_"+boost::lexical_cast<std::string>(f)] = epssq_[f];
    }
  }
}


// ---------------------------------------------------------------------------------------------


std::complex<double> SemicircleBandstructure::HilbertIntegral_PM(const std::complex<double>& zeta_A, const unsigned flavor) const {
  return zeta_A*(1.-std::sqrt(1.-4.*epssq_[flavor]/zeta_A/zeta_A))/2./epssq_[flavor];
}

std::complex<double> SemicircleBandstructure::HilbertIntegral_AFM(const std::complex<double>& zeta_A_times_zeta_B, const unsigned flavor) const {
  return (1.-std::sqrt(1.-4.*epssq_[flavor]/zeta_A_times_zeta_B))/2./epssq_[flavor];
}


// ---------------------------------------------------------------------------------------------


template<class T>
class DOS_integrand {
public:
  typedef T result_type;
  DOS_integrand(const std::vector<double>& dos, double low_e, double bin_size) : dos_(dos), low_e_(low_e), bin_size_(bin_size) {}
  unsigned size() const { return dos_.size(); }
  virtual result_type operator()(unsigned) const { throw std::logic_error("ERROR: the program shall never come here."); }
  virtual ~DOS_integrand() {}
protected:
  const std::vector<double>& dos_;
  const double low_e_, bin_size_;
};

class DOS_integrand_0th_moment : public DOS_integrand<double> {
public:
  DOS_integrand_0th_moment(const std::vector<double>& dos, double low_e, double bin_size) : DOS_integrand<double>(dos,low_e,bin_size) {}
  double operator()(unsigned i) const { return dos_[i]; }
};
class DOS_integrand_1st_moment : public DOS_integrand<double> {
public:
  DOS_integrand_1st_moment(const std::vector<double>& dos, double low_e, double bin_size) : DOS_integrand<double>(dos,low_e,bin_size) {}
  double operator()(unsigned i) const { return dos_[i]*(low_e_+i*bin_size_); }
};
class DOS_integrand_2nd_moment : public DOS_integrand<double> {
public:
  DOS_integrand_2nd_moment(const std::vector<double>& dos, double low_e, double bin_size) : DOS_integrand<double>(dos,low_e,bin_size) {}
  double operator()(unsigned i) const { return dos_[i]*(low_e_+i*bin_size_)*(low_e_+i*bin_size_); }
};
class DOS_HilbertIntegrand_PM : public DOS_integrand<std::complex<double> > {
public:
  DOS_HilbertIntegrand_PM(const std::vector<double>& dos, double low_e, double bin_size, std::complex<double> zeta)
    : DOS_integrand<std::complex<double> >(dos,low_e,bin_size), zeta_(zeta) {}
  std::complex<double> operator()(unsigned i) const { return dos_[i]/(zeta_-low_e_-i*bin_size_); }
private:
  std::complex<double> zeta_;
};
class DOS_HilbertIntegrand_AFM : public DOS_integrand<std::complex<double> > {
public:
  DOS_HilbertIntegrand_AFM(const std::vector<double>& dos, double low_e, double bin_size, std::complex<double> zeta0_times_zeta1)
    : DOS_integrand<std::complex<double> >(dos,low_e,bin_size), zeta0_times_zeta1_(zeta0_times_zeta1) {}
  std::complex<double> operator()(unsigned i) const { return dos_[i]/(zeta0_times_zeta1_-(low_e_+i*bin_size_)*(low_e_+i*bin_size_)); }
private:
  std::complex<double> zeta0_times_zeta1_;
};


DOSBandstructure::DOSBandstructure(const alps::Parameters& parms, bool verbose)
  : Bandstructure(parms),
    dos_(static_cast<unsigned>(parms.value_or_default("FLAVORS", 2))/2),
    dos_min_(static_cast<unsigned>(parms.value_or_default("FLAVORS", 2))/2),
    dos_step_(static_cast<unsigned>(parms.value_or_default("FLAVORS", 2))/2)
{
  std::ifstream dos_file(parms["DOSFILE"].c_str());
  if(!dos_file.good()) {
    std::cerr<<"ERROR: DOSBandstructure: problem to read DOSFILE = " << parms["DOSFILE"] << std::endl;
    throw std::runtime_error("DOSFILE is not good!");
  }
  unsigned n_bands = static_cast<unsigned>(parms.value_or_default("FLAVORS", 2))/2;
  double eps, d;
  if (verbose) std::cout<<"BANDSTRUCTURE:"<<std::endl<<"Using density of states loaded from: "<<parms["DOSFILE"]<<std::endl;
  std::vector<std::vector<double> > e_(n_bands);
  if (verbose) std::cout<<"(Note: Assuming equidistant energy intervals.)"<<std::endl;
  if (verbose && n_bands>1) std::cout<<"(Note: Assuming the same number of bins for all bands.)"<<std::endl;
  while(dos_file>>eps>>d){   
    e_[0].push_back(eps);
    dos_[0].push_back(d);
    for(unsigned band=1; band<n_bands; ++band){
      if(!(dos_file>>eps>>d))throw std::runtime_error("DOSFILE is not corrupt!");
      e_[band].push_back(eps);
      dos_[band].push_back(d);
    }
  }
  
  for(unsigned band=0; band<n_bands; ++band){
    dos_min_[band]=e_[band][0];
    dos_step_[band]=(e_[band][e_[band].size()-1]-dos_min_[band])/static_cast<double>(dos_[band].size()-1);
    double s, S;
    s=simpson_integrate(DOS_integrand_0th_moment(dos_[band], dos_min_[band], dos_step_[band]));
    S=0;
    for(unsigned i=0; i<dos_[band].size(); ++i){
      dos_[band][i]/=s;
      S+=dos_[band][i];
    }
    if (verbose) std::cout<<"check: total sum after normalization is: "<<S<<" (should be close to 1)"<<std::endl;

    //find first and second moment of band structure: simpson integrate: (e.g. for the 2nd) dos[n]*epsilon[n]*epsilon[n] 
    eps_[2*band]=simpson_integrate(DOS_integrand_1st_moment(dos_[band], dos_min_[band], dos_step_[band]));
    eps_[2*band+1]=eps_[2*band];
    epssq_[2*band]=simpson_integrate(DOS_integrand_2nd_moment(dos_[band], dos_min_[band], dos_step_[band]));
    epssq_[2*band+1]=epssq_[2*band];
    if (verbose) {
      std::cout<<"Flavors "<<2*band<<" and "<<2*band+1<<": first moment of bandstructure: "<<eps_[band]<<std::endl;
      std::cout<<"Flavors "<<2*band<<" and "<<2*band+1<<": second moment of bandstructure: "<<epssq_[band]<<std::endl;
    }
  }
}

std::complex<double> DOSBandstructure::HilbertIntegral_PM(const std::complex<double>& zeta_A, const unsigned flavor) const {
  return simpson_integrate(DOS_HilbertIntegrand_PM(dos_[flavor/2],dos_min_[flavor/2],dos_step_[flavor/2],zeta_A));
}

std::complex<double> DOSBandstructure::HilbertIntegral_AFM(const std::complex<double>& zeta_A_times_zeta_B, const unsigned flavor) const {
  return simpson_integrate(DOS_HilbertIntegrand_AFM(dos_[flavor/2],dos_min_[flavor/2],dos_step_[flavor/2],zeta_A_times_zeta_B));
}


// ---------------------------------------------------------------------------------------------


template <class Dispersion>
class TwoD_HilbertIntegrand_PM {
public:
  typedef std::complex<double> result_type;
  TwoD_HilbertIntegrand_PM(const Dispersion& f, std::complex<double> zeta) : f_(f), zeta_(zeta) {}
  result_type operator()(double kx, double ky) const { return 1./(zeta_-f_(kx,ky)); }
private:
  const Dispersion& f_;
  std::complex<double> zeta_;
};
template <class Dispersion>
class TwoD_HilbertIntegrand_AFM {
public:
  typedef std::complex<double> result_type;
  TwoD_HilbertIntegrand_AFM(const Dispersion& f, std::complex<double> zeta0_times_zeta1) : f_(f), zeta0_times_zeta1_(zeta0_times_zeta1) {}
  result_type operator()(double kx, double ky) const { double tmp=f_(kx,ky); return 1./(zeta0_times_zeta1_-tmp*tmp); }
private:
  const Dispersion& f_;
  std::complex<double> zeta0_times_zeta1_;
};


std::complex<double> SquareLatticeBandstructure::HilbertIntegral_PM(const std::complex<double>& zeta_A, const unsigned flavor) const {
  return twodsimpson(TwoD_HilbertIntegrand_PM<SquareLatticeBandstructure>(*this,zeta_A), kx_min, ky_min, kx_max, ky_max, L_) / k_area();
}

std::complex<double> SquareLatticeBandstructure::HilbertIntegral_AFM(const std::complex<double>& zeta_A_times_zeta_B, const unsigned flavor) const {
  return twodsimpson(TwoD_HilbertIntegrand_AFM<SquareLatticeBandstructure>(*this,zeta_A_times_zeta_B), kx_min, ky_min, kx_max, ky_max, L_) / k_area();
}



std::complex<double> HexagonalLatticeBandstructure::HilbertIntegral_PM(const std::complex<double>& zeta_A, const unsigned flavor) const {
  /* Note: as the model has 2 bands with dispersion +e(kx,ky) and -e(kx,ky), we may compute the PM HilbertIntegral as follows
   *      0.5 * sum_{n=1,2} int_{BZ} de kx de ky 1/(zeta_A-e(n,kx,ky))
   *    = 0.5 * int_{BZ} de kx de ky [1/(zeta_A-e(n=1,kx,ky)) + 1/(zeta_A+e(n=1,kx,ky))]
   *    = zeta_A int_{BZ} de kx de ky 1/(zeta_A^2 - e^2)
   * so we may use the AFM integral
   */
  return zeta_A*twodsimpson(TwoD_HilbertIntegrand_AFM<HexagonalLatticeBandstructure>(*this,zeta_A*zeta_A), kx_min, ky_min, kx_max, ky_max, L_) / k_area();
}

std::complex<double> HexagonalLatticeBandstructure::HilbertIntegral_AFM(const std::complex<double>& zeta_A_times_zeta_B, const unsigned flavor) const {
  /* Note: as the model has 2 bands with dispersion +e(kx,ky) and -e(kx,ky), we do not need the handle the 2 bands separately
   */
  return twodsimpson(TwoD_HilbertIntegrand_AFM<HexagonalLatticeBandstructure>(*this,zeta_A_times_zeta_B), kx_min, ky_min, kx_max, ky_max, L_) / k_area();
}


// ---------------------------------------------------------------------------------------------


boost::shared_ptr<Bandstructure> BandstructureFactory(const alps::Parameters& parms, bool verbose) {
  if (parms.defined("DOSFILE")) {
    boost::shared_ptr<Bandstructure> bandstruct(new DOSBandstructure(parms,verbose));
    return bandstruct;
  } else if (parms.defined("TWODBS") && parms["TWODBS"]=="hexagonal") {
    boost::shared_ptr<Bandstructure> bandstruct(new HexagonalLatticeBandstructure(parms,verbose));
    return bandstruct;
  } else if (parms.defined("TWODBS")) {
    boost::shared_ptr<Bandstructure> bandstruct(new SquareLatticeBandstructure(parms,verbose));
    return bandstruct;
  } else {
    boost::shared_ptr<Bandstructure> bandstruct(new SemicircleBandstructure(parms,verbose));
    return bandstruct;
  }
}


// static data member initialization

const double SquareLatticeBandstructure::kx_min=-M_PI;
const double SquareLatticeBandstructure::ky_min=-M_PI;
const double SquareLatticeBandstructure::kx_max=M_PI;
const double SquareLatticeBandstructure::ky_max=M_PI;


const double HexagonalLatticeBandstructure::kx_min=-M_PI;
const double HexagonalLatticeBandstructure::ky_min=-M_PI/2.;
const double HexagonalLatticeBandstructure::kx_max=M_PI;
const double HexagonalLatticeBandstructure::ky_max=M_PI/2.;
