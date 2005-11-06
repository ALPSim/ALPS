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

/* $Id$ */

#ifndef ALPS_DMFT_HILBERTTRANSFORMER_H
#define ALPS_DMFT_HILBERTTRANSFORMER_H

/// @file hilberttransformer.h
/// @brief Hilbert transformations
///
/// declares the abstract base class and concrete realizations for the Hilbert transformations
/// @sa HilbertTransformer, SemicircleHilbertTransformer, SquarelatticeHilbertTransformer, GeneralHilbertTransformer
///
/// @todo the transformations need to be implemented. responsible: Philipp Werner and Hartmut Monien

#include "solver.h"
#include <functional>

#include "fouriertransform.h"

/// @brief performs a Hilbert transformation
///
/// The HilbertTransformer performs a Hilbert transformation for the self energy and density of states
//class HilbertTransformer : public std::binary_function<const vector_type&,double,vector_type>
class HilbertTransformer : public std::binary_function<const pair_vector_type&,double,pair_vector_type>
{
public:
/// the function call operator performs a Hilbert transformation of the self energy 
/// and chemical potential given as parameters. The density of states is constant for
/// each object and usually specified in the constructor
///
/// @param Sigma the frequency dependent self energy as a function of Matsubara frequency
/// @param mu the chemical potential, h the magnetic field
/// @return the result of the Hilbert transform of the self energy (equation to be added)
  //virtual result_type operator()(first_argument_type Sigma, second_argument_type mu)=0;
  virtual result_type operator()(first_argument_type G, second_argument_type mu, third_argument_type h)=0; 
  virtual ~HilbertTransformer() {}
};

/// A Hilbert transformation for a semcircle density of states
/// It receives G(\tau) as input and returns G0(\tau)
class SemicircleHilbertTransformer : public HilbertTransformer {
public:
/// the constructor accepts the bandwidth
  SemicircleHilbertTransformer(double bandwidth) : bw_(bandwidth) {}
  //result_type operator()(first_argument_type Sigma, second_argument_type mu)
  result_type operator()(first_argument_type G, second_argument_type mu, third_argument_type h)
  {

  double BETA=1; // shall we use this convention?
  std::pair<std::vector<std::complex<double> >, std::vector<std::complex<double> > > G_omega(G.size());
  std::pair<std::vector<std::complex<double> >, std::vector<std::complex<double> > > G0_omega(G.size());
  
  forward_FT(BETA, bw_, mu, h, G, G_omega);
  for(int i=0; i<G_omega.size(); i++) {
    complex<double> iw(0.,(2*i+1)*M_PI/BETA);
	G0_omega.first[i] = iw + mu + h/2. - std::pow(bw_,2)*G_omega.second[i];
	G0_omega.second[i] = iw + mu - h/2. - std::pow(bw_,2)*G_omega.first[i];
  }
  
  first_argument_type G0(G.size());
  backward_FT(BETA, bw_, mu, h, G0, G0_omega);
   
  return G0;
   
  // to implement
  //return vector_type();
  }
private:
  double bw_; // or band width
};

/// A Hilbert transformation for a square lattice density of states
/// @todo the transformation need to be implemented. responsible: Philipp Werner and Hartmut Monien
class SquarelatticeHilbertTransformer : public HilbertTransformer {
public:
/// the constructor accepts the bandwidth
  SquarelatticeHilbertTransformer(double bandwidth) : bw_(bandwidth) {}
  result_type operator()(first_argument_type Sigma, second_argument_type mu)
  {
    //to implement
  
  return vector_type();
  }
private:
  double bw_; // or band width
};

// ...

/// A Hilbert transformation for a general density of states, specified in the constructor
/// @todo the transformation need to be implemented. responsible: Philipp Werner
class GeneralHilbertTransformer : public HilbertTransformer {
public:
/// the constructor accepts the desity of states
/// @param D the density of states. 
///          Note that issues such as units are still undecided and might require further 
///          additions to the interface
  GeneralHilbertTransformer(const vector_type& D) : D_(D) {}
  result_type operator()(first_argument_type Sigma, second_argument_type mu)
  {
    //to implement
  
  return vector_type();
  }
private:
  vector_type D_; // or band width
};


#endif
