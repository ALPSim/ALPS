/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *               2012 - 2013 by Jakub Imriska <jimriska@phys.ethz.ch>
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
/// @sa HilbertTransformer, SemicircleHilbertTransformer, FrequencySpaceHilbertTransformer
///

#include "types.h"
#include "bandstructure.h"
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
                                            double mu, double h, double beta) const=0;
  itime_green_function_t symmetrize(const itime_green_function_t& G_tau, const bool symmetrization) const;
  virtual itime_green_function_t initial_G0(const alps::Parameters& parms) const;
  virtual ~HilbertTransformer() {}
};



/// A Hilbert transformation for a semicircle density of states
/// It receives G(\tau) as input and returns G0(\tau)
class SemicircleHilbertTransformer : public HilbertTransformer 
{
public:
  /// the constructor accepts the bandwidth
  SemicircleHilbertTransformer(alps::Parameters& parms) 
    : bethe_parms(parms,true)
  {
    bethe_parms.set_parms(parms);
  }
  
  ///operator() implements abstract virtual operator() of base class HilbertTransformer 
  ///and performs the actual Hilbert transformation.
  itime_green_function_t operator()(const itime_green_function_t& G_tau, 
                                    double mu, double h, double beta) const;
  itime_green_function_t initial_G0(const alps::Parameters& parms) const;
  
private:
  SemicircleBandstructure bethe_parms;
};






/// @brief performs a Hilbert transformation
///
/// The FrequencySpaceHilbertTransformer performs a Hilbert transformation for the self energy and density of states.
/// Arguments are expected to be in Matsubara Frequencies.
class FrequencySpaceHilbertTransformer{
public:
  
  virtual ~FrequencySpaceHilbertTransformer() {}
  
  /// the function call operator performs a Hilbert transformation of the self energy 
  /// and chemical potential given as parameters.
  /// The density of states is specific for each of the derived classes
  ///
  /// @param G_omega the Greens function as a function of Matsubara Frequency omega 
  /// @param mu the chemical potential, h the magnetic field, beta the inverse temperature
  /// @return the result of the Hilbert transform: the bare Green's function G0 in Matsubara frequencies
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t & G_omega, 
                                                matsubara_green_function_t &G0_omega, 
                                                const double mu, const double h, const double beta) const=0;
  virtual matsubara_green_function_t initial_G0(const alps::Parameters& parms) const=0;
  
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


/// @brief performs a Hilbert transformation
///
/// The density of states is handled via class Bandstructure: 
///    semicircle DOS
///    user-defined via DOS histogram
///    tight-binding for square or hexagonal lattice
/// Arguments are expected to be in Matsubara Frequencies.
class GeneralFSHilbertTransformer : public FrequencySpaceHilbertTransformer {
public:
  
  GeneralFSHilbertTransformer(const alps::Parameters& parms, bool /*ignored*/);
  GeneralFSHilbertTransformer(alps::Parameters& parms);
  virtual ~GeneralFSHilbertTransformer() {}
  
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t & G_omega, 
                                                matsubara_green_function_t &G0_omega, 
                                                const double mu, const double h, const double beta) const;
  virtual matsubara_green_function_t initial_G0(const alps::Parameters& parms) const;

private:
  bool AFM;
  boost::shared_ptr<Bandstructure> bandstruct;
};


/// A Hilbert transformation for a semicircle density of states
/// Currently UNUSED: the FrequencySpaceHilbertTransformer is able to handle the semicircle DOS
/// NOTE: that the SemicircleFSHilbertTransformer::operator() uses the equation t^2*Delta=G for converged solution,
///       thus the result of the Hilbert transformation does not necesarilly equal to that of the 
///       FrequencySpaceHilbertTransformer::operator() in the not-yet-converged case
///       The effect on the convergency rate not fully explored.
class SemicircleFSHilbertTransformer : public FrequencySpaceHilbertTransformer {
public:
  SemicircleFSHilbertTransformer(alps::Parameters& parms)
    : bandstruct(parms) 
  {
    bandstruct.set_parms(parms);
  }

  /// It receives the dressed Green's function G(\omega) as input and returns tha bare GF G0(\omega)
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t& G_omega, 
                                                matsubara_green_function_t &G0_omega_ignored, 
                                                const double mu, const double h, const double beta) const;
  
  virtual matsubara_green_function_t initial_G0(const alps::Parameters& parms) const;
  
private:
  SemicircleBandstructure bandstruct;
};



#endif /*ALPS_DMFT_HILBERTTRANSFORMER_H*/
