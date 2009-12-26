 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/
#ifndef INTEGRANDH
#define INTEGRANDH
#include<complex>
double dispersion(double kx, double ky, double t, double tprime);
class integrand{
  public:
    virtual std::complex<double>operator()(double kx, double ky) const=0;
    virtual ~integrand(){};
};
class afm_integrand:public integrand{ //class for the integrand with an epsilon(k) that uses the t-band structure.
  public:
  afm_integrand(double t, double tprime, std::complex<double> zeta_A, std::complex<double> zeta_B){
    zeta_A_=zeta_A;
    zeta_B_=zeta_B;
    t_=t;
    tprime_=tprime;
  }
  std::complex<double> operator()(double kx, double ky) const{
    double eps=dispersion(kx,ky,t_,tprime_);
    return 1./(zeta_A_*zeta_B_-eps*eps);
  }
  private:
  double t_;
  double tprime_;
  std::complex<double> zeta_A_;
  std::complex<double> zeta_B_;
};


std::complex<double> twodsimpson(const integrand &, double ax, double ay, double bx, double by, int N);
#endif
