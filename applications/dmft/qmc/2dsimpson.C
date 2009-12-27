 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/
 
#include<iostream>
#include "2dsimpson.h"
//implementation of 2d simpson for DCA and Single site AFM self consistency
std::complex<double> twodsimpson(const integrand &f, double ax, double ay, double bx, double by, int N){
  std::complex<double>       result=0.;
  double       h=(bx-ax)/(2.*N);
  double       k=(by-ay)/(2.*N);

  result =  f(ax,ay) + f(ax,by) + f(bx, ay) +f(bx,by) ;
            
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
  for(int i=1;i<N;++i){
    for(int j=1;j<=N;++j){
      result+=8.*f(ax+h*(2*i-1), ay+k*(2*j));
    }
  }
  for(int i=1;i<=N;++i){
    for(int j=1;j<N;++j){
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

double dispersion(double kx, double ky, double t, double tprime){
 return -2.*t*(cos(kx)+cos(ky)) -4.*tprime*cos(kx)*cos(ky);
}

