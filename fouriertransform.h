/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
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

#include "hirschfyesim.h"
#include <complex>
#include "solver.h" // for the multiple_vector_type


void generate_spline_matrix(dense_matrix & spline_matrix);

//example:  V = vector_type
template <class V> void evaluate_second_derivatives(double dt, double bw, double mu, double h, dense_matrix & spline_matrix, V & g, V & second_derivatives) {

  // g, rhs and second_derivatives have dimension N+1
  int Np1 = spline_matrix.size1();
  double dtsq=dt*dt;
  double bwsq=bw*bw;
  // rhs is the vector containing the data of the curve y = g(tau), which allows to 
  // compute the vector of second derivatives y'' at times tau_n by evaluating
  // y'' = spline_matrix * rhs(y)
  //
  //                   (G''(0)+G''(beta))*(dt)^2                                
  //                   y0 - 2*y1 + y2
  //                   y1 - 2*y2 + y3
  // rhs = 6/(dt)^2 * ...
  //
  //                   yNp1-3 - 2*yNp1-2 + yNp1-1
  //                   (G'(0)+G'(beta))*(dt)^2 + y0 - y1 + yNp1-2 - yNp1-1     

  V rhs;
  
  rhs.reserve(Np1);

  rhs[0] = bwsq*dtsq; 
  for (int i=1; i<Np1-1; i++) {
    rhs[i] = g[i-1]-2*g[i]+g[i+1];
  }
  rhs[Np1-1] = dtsq*(mu+h/2.) + g[0] - g[1] + g[Np1-2] - g[Np1-1];

  for (int i=0; i<Np1; i++) {
    for (int j=0; j<Np1; j++) {
      second_derivatives[i] += spline_matrix(i,j)*rhs[j];
	}
	second_derivatives[i] *= 6./dtsq;
  }
}

//example class V type:  vector_type.
//example class W type:  complex_vector_type
//is it really meaningful to have this as a template function?!?
template <class V, class W> void forward_FT(double beta, double bw, double mu, double h, V & v, W & v_omega) {


  int N=v.size();
  double dt = beta/N;
  dense_matrix spline_matrix(N+1, N+1);
  generate_spline_matrix(spline_matrix);
  
  // v2 is the matrix containing the second derivatives y'' of interpolated y=v[tau] at points tau_n	 
  V v2;
  v2.reserve(N+1);	//allocate memory for both vector. We will need a vector class that implements reserve()
 
  v.push_back(1-v[0]); 
  evaluate_second_derivatives(dt, bw, mu, h, spline_matrix, v, v2);

  for (int k=0; k<N; k++) {
	std::complex<double> iw(0, (2*k+1)*M_PI/beta);
	v_omega[k] = 0;
	
	for (int n=1; n<N; n++) {
	  v_omega[k] += exp(iw*(n*dt))*(v2[n+1]-2*v2[n]+v2[n-1]);
	}
	v_omega[k] += v2[1] - v2[0] + v2[N] - v2[N-1];
	v_omega[k] *= 1./(dt*pow(iw,4));
  }
  v.pop_back();
} 

inline std::complex<double> f_omega(std::complex<double> iw, double mu) {
  return 1./iw + mu/(iw*iw);
}

inline double f_tau(double tau_beta, double mu) {
  return 1 + 0.25*(1-2*tau_beta);
}

//example class V type:  vector_type.
//example class W type:  complex_vector_type
template <class V, class W> void backward_FT(double beta, double bw, double mu, double h, V & v, W v_omega) {
  
  int N = v.size();
  double dt = beta/N;
  
  for (int k=0; k<N; k++) {
    std::complex<double> iw(0,(2*k+1)*M_PI/beta);
    v_omega[k] -= f_omega(iw, mu+h/2.);
  }	
  
  for (int i=0; i<N; i++) {
    v[i] = f_tau(i*dt, mu+h/2.);
	for (int k=0; k<N; k++) {
	  double wt((2*k+1)*i*M_PI/N);
	  v[i] += 2/beta*(cos(wt)*v_omega[k].real()+sin(wt)*v_omega[k].imag());
	}
  }  
}
