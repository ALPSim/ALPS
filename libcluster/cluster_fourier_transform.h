/*****************************************************************************
 *
 * Copyright (C) 2005 - 2015 Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/
#pragma once
#include <complex>

#include "fourier_transform.h"
#include "cluster.h"

class ClusterG0FourierTransformer : public FourierTransformer
{
public: 
  ClusterG0FourierTransformer(const double beta, const double mu, const double h, const int n_flavor, const int n_site, 
                              const std::vector<std::vector<double> >& eps, const std::vector<std::vector<double> >& epssq)
  : FourierTransformer(beta, n_flavor, n_site)
  { 
    if (n_site%2!=0 && h!=0) 
      boost::throw_exception(std::logic_error("FourierTransformer: cannot handle antiferromagnetism with odd cluster size\n"));
    for(int f=0; f<n_flavor; f++) {
      int s = f % 2 ? -1 : 1;
      for(int i=0; i<n_site; i++) {
        c1_[f][i][i] = 1.;
        c2_[f][i][i] = eps[f][i] - mu;
        c3_[f][i][i] = epssq[f][i] - 2*mu*eps[f][i] + mu*mu +h*h/4.;
      }
      if (n_site%2==0) {
        for(int i=0; i<n_site; i+=2) {
          c2_[f][i+1][i] = s*h/2;
          c2_[f][i][i+1] = c2_[f][i+1][i];
          c3_[f][i+1][i] = s*h/2*(eps[f][i] + eps[f][i+1] - 2*mu);
          c3_[f][i][i+1] = c3_[f][i+1][i];
        }
      }
    }
   /* for(int f=0;f<n_flavor;f++){
      for(int i=0;i<n_site;i++){
        for(int j=0;j<n_site;j++){
          std::cout<<f<<i<<j<<" "<<c1_[f][i][j]<<" "<<c2_[f][i][j]<<" "<<c3_[f][i][j]<<std::endl;
        }
      }
    }*/
  }
};

class ClusterAntiFerroGFourierTransformer : public FourierTransformer
{
public: 
  ClusterAntiFerroGFourierTransformer(const double beta, const double mu, const double h, const double U, 
                                      const int n_flavor, const int n_site, 
                                      const std::vector<double>& density, const std::vector<double>& magnetization, 
                                      const std::vector<std::vector<double> >& eps, const std::vector<std::vector<double> >& epssq)
  : FourierTransformer(beta, n_flavor, n_site)
  { 
    for(int f=0;f<n_flavor;++f){
      int fbar = f%2==0 ? f+1 : f-1;
      int s = f % 2 ? -1 : 1;
      for(int i=0;i<n_site;++i){
        c1_[f][i][i] = 1.;
        c2_[f][i][i] = eps[f][i] - mu + U*density[fbar]; 
        c3_[f][i][i] = epssq[f][i] - 2.*mu*eps[f][i] + mu*mu + h*h/4 
        + 2.*U*density[fbar]*(eps[f][i]-mu) + U*U*density[fbar] + U*s*h*magnetization[fbar];
        Sc0_[f][i][i] = U * (density[fbar]-0.5);
        Sc1_[f][i][i] = U*U * (density[fbar]*(1-density[fbar]) - magnetization[fbar]*magnetization[fbar]);
        /*double m = magnetization[fbar];
         double n = density[fbar];
         double e = eps[f][i];
         double eq = i%2 ? eps[f][i-1] : eps[f][i+1];*/
        Sc2_[f][i][i] = 0.;/*U*U*n*(1-n) * (U*(1-n) + (e-mu))
                            + U*U*m * (2*U*(n-1)*m + U*m*n - (eq-mu)*n + eq)
                            + U*s*h*m* ( U*(0.5-n) + (e-eq) );*/
      }
      for(int i=0; i<n_site; i+=2) {
        c2_[f][i+1][i] = s*h/2 + U*magnetization[fbar];
        c2_[f][i][i+1] = c2_[f][i+1][i];
        c3_[f][i+1][i] = (U*magnetization[fbar]+s*h/2)*(eps[f][i] + eps[f][i+1] - 2*mu)
        + U*U*magnetization[fbar] + U*s*h*density[fbar];
        c3_[f][i][i+1] = c3_[f][i+1][i];
        Sc0_[f][i+1][i] = U*magnetization[fbar];
        Sc0_[f][i][i+1] = Sc0_[f][i+1][i];
        Sc1_[f][i+1][i] = U*U*magnetization[fbar]*(1-2*density[fbar]);
        /*double m = magnetization[fbar];
         double n = density[fbar];
         double e = eps[f][i];
         double eq = i%2 ? eps[f][i-1] : eps[f][i+1];*/
        Sc2_[f][i][i+1] = 0;/*U*U*U * (m*m*m + m*(3*n*n-4*n+1)) 
                             - 0.5*U*U*s*h* (n*(1-n) + m*m)
                             +U*U*(((2*mu-eq-e)*n-mu)*m + 0*(e+eq));*/
        Sc2_[f][i][i+1] = Sc2_[f][i+1][i];	
      }
    }
    /*std::ofstream forstr("for.dat");
     for(int f=0;f<n_flavor;++f){
     std::cout << "density: " << density[f] << std::endl; 
     for(int s=0;s<n_site;s+=1){
     for (int p1=0; p1<=0; ++p1) {
     for (int p2=0; p2<=0; ++p2) {
     for(uint k=0; k<100; k++){
     std::complex<double> iw(0,(2*k+1)*M_PI/beta_);
     forstr << iw.imag() << "\t" << f << "\t" << s+p1 << "\t" << s+p2 << "\t"
     << f_omega(iw, c1_[f][s+p1][s+p2],c2_[f][s+p1][s+p2], c3_[f][s+p1][s+p2]).real() << "\t"
     << f_omega(iw, c1_[f][s+p1][s+p2],c2_[f][s+p1][s+p2], c3_[f][s+p1][s+p2]).imag()
     << std::endl;
     }
     }
     }
     }
     }*/
  }
  
};
