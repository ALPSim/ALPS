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


void generate_spline_matrix(dense_matrix & spline_matrix);


//void evaluate_second_derivatives(double dt, double bw, double mu, double h, dense_matrix & spline_matrix, multiple_vector_type & g, multiple_vector_type & second_derivatives) {
template <class V> void evaluate_second_derivatives(double dt, double bw, double mu, double h, dense_matrix & spline_matrix, V & g, V & second_derivatives) {

  // g, rhs and second_derivatives have dimension N+1
  int Np1 = spline_matrix.size1();
  
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

  V rhs(Np1);
  //multiple_vector_type rhs(Np1);
  
  rhs.first[0] = std::pow(bw*dt,2); 
  rhs.second[0] = rhs.first()[0];
  for (int i=1; i<Np1-1; i++) {
    rhs.first[i] = g.first[i-1]-2*g.first[i]+g.first[i+1];
	rhs.second[i]= g.second[i-1]-2*g.second[i]+g.second[i+1];
  }
  rhs.first[Np1-1] = std::pow(dt,2)*(mu+h/2.) + g.first[0] - g.first[1] + g.first[Np1-2] - g.first[Np1-1];
  rhs.second[Np1-1] = std::pow(dt,2)*(mu-h/2.) + g.second[0] - g.second[1] + g.second[Np1-2] - g.second[Np1-1];

  for (int i=0; i<Np1; i++) {
    for (int j=0; j<Np1; j++) {
      second_derivatives.first[i] += spline_matrix(i,j)*rhs.first()[j];
	  second_derivatives.second[i] += spline_matrix(i,j)*rhs.second()[j];
	}
	second_derivatives.first[i] *= 6./std::pow(dt,2);
	second_derivatives.second[i] *= 6./std::pow(dt,2);
  }
}



template <class V, class W> void forward_FT(double BETA, double bw, double mu, double h, V & v, W & v_omega) {

  int N=v.size();
  double dt = BETA/N;
  dense_matrix & spline_matrix(N+1, N+1);
  generate_spline_matrix(spline_matrix);
  
  // v2 is the matrix containing the second derivatives y'' of interpolated y=v[tau] at points tau_n	 
  V v2(N+1);
  //multiple_vector_type v2(N+1);
  std::pair<double, double> last(1-v.first[0], 1-v.second[0]);
  v.push_back(last); 
  evaluate_second_derivatives(dt, bw, mu, h, spline_matrix, v, v2);

  for (int k=0; k<N; k++) {
	std::complex<double> iw(0, (2*k+1)*M_PI/BETA);
	v_omega.first[k] = 0;
	v_omega.second[k] = 0;
	
	for (int n=1; n<N; n++) {
	  v_omega.first[k] += exp(iw*(n*dt))*(v2.first[n+1]-2*v2.first[n]+v2.first[n-1]);
	  v_omega.second[k] += exp(iw*(n*dt))*(v2.second[n+1]-2*v2.second[n]+v2.second[n-1]);
	}
	v_omega.first[k] += v2.first[1] - v2.first[0] + v2.first[N] - v2.first[N-1];
	v_omega.second[k] += v2.second[1] - v2.second[0] + v2.second[N] - v2.second[N-1];
	v_omega.first[k] *= 1./(dt*pow(iw,4));
	v_omega.second[k] *= 1./(dt*pow(iw,4));
  }
  v.pop_back();
} 



template <class V, class W> void backward_FT(double BETA, double bw, double mu, double h, V & v, W & v_omega) {
  
  int Np1 = v.size();
  int N = Np1-1;
  double dt = BETA/N;
  
  for (int i=0; i<Np1; i++) {
    v.first[i] = 0;
	v.second[i] = 0;
	for (int k=0; k<N; k++) {
	  double wt((2*k+1)*i*M_PI/N);
	  v.first[i] += 2/BETA*(cos(wt)*v_omega.first[k].real()+sin(wt)*v_omega.first[k].imag());
	  v.second[i] += 2/BETA*(cos(wt)*v_omega.second[k].real()+sin(wt)*v_omega.second[k].imag());
	}
  }  
}
