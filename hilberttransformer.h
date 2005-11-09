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
class HilbertTransformer
{
public:
/// the function call operator performs a Hilbert transformation of the self energy 
/// and chemical potential given as parameters. The density of states is constant for
/// each object and usually specified in the constructor
///
/// @param G the Greens function as a function of tau 
/// @param mu the chemical potential, h the magnetic field, beta the inverse temperature
/// @return the result of the Hilbert transform (G0)
  virtual multiple_vector_type operator()(multiple_vector_type& G, double mu, double h, double beta)=0; 
  virtual ~HilbertTransformer() {}
};

/// A Hilbert transformation for a semcircle density of states
/// It receives G(\tau) as input and returns G0(\tau)
class SemicircleHilbertTransformer : public HilbertTransformer {
public:
/// the constructor accepts the bandwidth
  SemicircleHilbertTransformer(double bandwidth) : bw_(bandwidth) {}
  multiple_vector_type operator()(multiple_vector_type& G, double mu, double h, double beta)
  {

  multiple_vector_type G0(G);
/*
  /// this needs more work: do you really want complex vectors here?
  multiple_vector_type G_omega(G); // remember these are pairs now!,
*/  
  multiple_complex_vector_type G_omega;
  G_omega.first.resize(G.first.size()); // remember these are pairs now!
  G_omega.second.resize(G.second.size()); // remember these are pairs now!
  
  multiple_complex_vector_type G0_omega(G_omega); 
  
  forward_FT(beta, bw_, mu, h, G, G_omega);
  for(int i=0; i<G_omega.first.size(); i++) {
    std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
	G0_omega.first[i] = iw + mu + h/2. - std::pow(bw_,2)*G_omega.second[i];
	G0_omega.second[i] = iw + mu - h/2. - std::pow(bw_,2)*G_omega.first[i];
  }
  
  backward_FT(beta, bw_, mu, h, G0, G0_omega);
   
  return G0;
   
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
  multiple_vector_type operator()(multiple_vector_type Sigma, double mu)
  {
    //to implement
  
  return multiple_vector_type();
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
  multiple_vector_type operator()(multiple_vector_type Sigma, double mu)
  {
    //to implement
  
  return multiple_vector_type();
  }
private:
  vector_type D_; // or band width
};


#endif
