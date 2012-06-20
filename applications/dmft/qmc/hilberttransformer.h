 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *               2012        by Jakub Imriska <jimriska@phys.ethz.ch>
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

/* $Id: hilberttransformer.h 367 2009-08-05 10:04:31Z fuchs $ */

#ifndef ALPS_DMFT_HILBERTTRANSFORMER_H
#define ALPS_DMFT_HILBERTTRANSFORMER_H


/// @file hilberttransformer.h
/// @brief Hilbert transformations
///
/// declares the abstract base class and concrete realizations for the Hilbert transformations
/// @sa HilbertTransformer, SemicircleHilbertTransformer, SquarelatticeHilbertTransformer, 
///     GeneralHilbertTransformer, DCATransformer
///
/// @todo the transformations need to be implemented. responsible: Philipp Werner and Hartmut Monien

#include "types.h"
#include "2dsimpson.h"
#include "fouriertransform.h"

/// @brief performs a Hilbert transformation
///
/// The HilbertTransformer performs a Hilbert transformation for the self energy and density of states
class HilbertTransformer
{
public:
  /// the function call operator performs a Hilbert transformation of the self energy 
  /// and chemical potential given as parameters. The density of states is constant for
  /// each object and usually specified in the constructor of a derived class
  ///
  /// @param G_tau the Greens function as a function of imaginary time tau 
  /// @param mu the chemical potential, h the magnetic field, beta the inverse temperature
  /// @return the result of the Hilbert transform (G0_tau)
  virtual itime_green_function_t operator()(const itime_green_function_t& G_tau, 
                                            double mu, double h, double beta)=0; 
  itime_green_function_t symmetrize(const itime_green_function_t& G_tau, const bool symmetrization) const;
  virtual itime_green_function_t initial_G0(const alps::Parameters& parms);
  virtual ~HilbertTransformer() {}
};



/// A Hilbert transformation for a semcircle density of states
/// It receives G(\tau) as input and returns G0(\tau)
class SemicircleHilbertTransformer : public HilbertTransformer 
{
public:
  /// the constructor accepts the bandwidth
  SemicircleHilbertTransformer(double t) : t_(t) {}
  
  ///operator() implements abstract virtual operator() of base class HilbertTransformer 
  ///and performs the actual Hilbert transformation.
  itime_green_function_t operator()(const itime_green_function_t& G_tau, 
                                    double mu, double h, double beta);
  itime_green_function_t initial_G0(const alps::Parameters& parms);
  
private:
  double t_; 
};



/// A Hilbert transformation for a square lattice density of states
/// @todo the transformation need to be implemented. responsible: Philipp Werner and Hartmut Monien
class SquarelatticeHilbertTransformer : public HilbertTransformer {
public:
  /// the constructor accepts the bandwidth
  SquarelatticeHilbertTransformer(double bandwidth) : bw_(bandwidth) {}
  virtual itime_green_function_t operator()(const itime_green_function_t& G_tau, 
                                            double mu, double h, double beta)
  {
    //to implement
    std::cout<<"Square Lattice Hilbert transformer"<<G_tau<<std::endl;
    std::cout<<mu;
    std::cout<<beta;
    std::cout<<h;
    throw std::logic_error("Square Lattice Hilbert transformer: this function is NOT implemented. aborting.");
    return itime_green_function_t(0,0,0);
  }
private:
  double bw_; // or band width
};



/// A Hilbert transformation for a general density of states, specified in the constructor
/// @todo the transformation need to be implemented. responsible: Philipp Werner
class GeneralHilbertTransformer : public HilbertTransformer {
public:
  /// the constructor accepts the desity of states
  /// @param D the density of states. 
  ///          Note that issues such as units are still undecided and might require further 
  ///          additions to the interface
  GeneralHilbertTransformer(const vector_type& D) : D_(D) {}
  virtual itime_green_function_t operator()(const itime_green_function_t& G_tau, 
                                            double mu, double h, double beta)
  {
    //to implement
    std::cout<<G_tau<<std::endl;
    std::cout<<mu;
    std::cout<<beta;
    std::cout<<h;
    throw std::logic_error("General Hilbert transformer: this function is NOT implemented. aborting.");
    return itime_green_function_t(0,0,0);
  }
private:
  vector_type D_; // or band width
};



// definition of Hilbert transformers living in Fourier (Matsubara) space



/// @brief performs a Hilbert transformation
///
/// The FrequencySpaceHilbertTransformer performs a Hilbert transformation for the self energy and density of states.
/// Arguments are expected to be in Matsubara Frequencies.
/// See FrequencySpaceHilbertTransformer for a class that takes its arguments in imaginary time space.
class FrequencySpaceHilbertTransformer{
public:
  /// the function call operator performs a Hilbert transformation of the self energy 
  /// and chemical potential given as parameters. The density of states is constant for
  /// each object and usually specified in the constructor of a derived class
  ///
  /// @param G_omega the Greens function as a function of Matsubara Frequency omega 
  /// @param mu the chemical potential, h the magnetic field, beta the inverse temperature
  /// @return the result of the Hilbert transform: the bare Green's function G0 in Matsubara frequencies
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t & G_omega, 
                                                matsubara_green_function_t &G0_omega, 
                                                const double mu, const double h, 
                                                const double beta)=0;
  virtual matsubara_green_function_t initial_G0(const alps::Parameters& parms);
  virtual ~FrequencySpaceHilbertTransformer() {}
  template <class T>
  green_function<T> symmetrize(const green_function<T>& G, const bool symmetrization) const
  {
    green_function<T> G_new(G);
    if (symmetrization) {
      assert(G_new.nflavor()%2==0);
      for(spin_t flavor=0;flavor<G_new.nflavor(); flavor+=2){
        for(itime_index_t tau=0;tau<G_new.ntime();++tau){
          G_new(tau, flavor  )=0.5*(G_new(tau, flavor)+G_new(tau, flavor+1));
          G_new(tau, flavor+1)=G_new(tau, flavor);
        }
      }
    }
    return G_new;
  }
};



/// A Hilbert transformation for a semcircle density of states
/// It receives the dressed Green's function G(\omega) as input and returns tha bare GF G0(\omega)
class FSSemicircleHilbertTransformer : public FrequencySpaceHilbertTransformer {
public:
  /// the constructor accepts the bandwidth
  FSSemicircleHilbertTransformer(double t) : t_(t) {}
  ///operator() implements abstract virtual operator() of base class HilbertTransformer 
  ///and performs the actual Hilbert transformation.
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t& G_omega, 
                                                matsubara_green_function_t &G0_omega_ignored, 
                                                const double mu, const double h, const double beta);
  
private:
  double t_; 
};



class FSHamiltonianHilbertTransformer: public FrequencySpaceHilbertTransformer{
public:
  /// the constructor reads in the LDA Hamiltonian from a file
  FSHamiltonianHilbertTransformer(const alps::Parameters& parms);
  /// take G_iomegan as the input and construct G0_iomegan out of it. Do this
  /// via summation in k-space. The k-points are defined in the Hamiltonian.
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega, 
                                                matsubara_green_function_t &G0_omega, 
                                                const double mu, const double h, const double beta);
};



class FSDOSHilbertTransformer: public FrequencySpaceHilbertTransformer{
public:
  /// the constructor reads in the DOS from a file
  FSDOSHilbertTransformer(alps::Parameters& parms);
  /// take G_iomegan as the input and construct G0_iomegan out of it. Do this
  /// using integration over D(epsilon).
  /// We perform the HT G <- \int \frac{D(e)}{A - e} de
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega, 
                                                matsubara_green_function_t &G0_omega, 
                                                const double mu, const double h, const double beta);
  
protected:

  std::vector<double> epsilon;
  std::vector<double> dos;
};


//this can do Antiferromagnetism according to Georges et. al.
class AFM_FSDOSHilbertTransformer: public FSDOSHilbertTransformer{
public:
  AFM_FSDOSHilbertTransformer(alps::Parameters &parms)
  :FSDOSHilbertTransformer(parms){}
  virtual ~AFM_FSDOSHilbertTransformer(){}
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                                matsubara_green_function_t &G0_omega, 
                                                const double mu, const double h, const double beta);
};


// for 2-dimensional Hubbard (more than nearest-neighbor hoppings possible), currently for square and hexagonal lattice
/// We perform the HT G <- \int \frac{1}{A - e(kx,ky)} dkx dky
class TwoDHilbertTransformer: public FrequencySpaceHilbertTransformer{
public:
  TwoDHilbertTransformer(alps::Parameters& parms);
  
  matsubara_green_function_t operator() (const matsubara_green_function_t & G_omega, 
                                        matsubara_green_function_t &G0_omega, 
                                        const double mu, const double h, const double beta);
  
protected:
  TwoDBandstructure bandstruct_;   // there is the bandstructure information
  const int L_;
};


class TwoDAFMHilbertTransformer: public TwoDHilbertTransformer{
public:
  TwoDAFMHilbertTransformer(alps::Parameters& parms)
    : TwoDHilbertTransformer(parms)
    {}
  
  matsubara_green_function_t operator() (const matsubara_green_function_t & G_omega, 
                                        matsubara_green_function_t &G0_omega, 
                                        const double mu, const double h, const double beta);
};


#endif /*ALPS_DMFT_HILBERTTRANSFORMER_H*/
