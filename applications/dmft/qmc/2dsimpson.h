 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
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
