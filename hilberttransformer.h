/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
*                       Matthias Troyer <troyer@comp-phys.org>
*                       Emanuel Gull <gullc@comp-phys.org>
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

#include "types.h"
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
	/// each object and usually specified in the constructor of a derived class
	///
	/// @param G_tau the Greens function as a function of imaginary time tau 
	/// @param mu the chemical potential, h the magnetic field, beta the inverse temperature
	/// @return the result of the Hilbert transform (G0_tau)
	virtual itime_green_function_t operator()(itime_green_function_t& G_tau, FourierTransformer &fourier_ignored, double mu, double h, double beta, bool paramagnet)=0; 
	virtual ~HilbertTransformer() {}
};

/// A Hilbert transformation for a semcircle density of states
/// It receives G(\tau) as input and returns G0(\tau)
class SemicircleHilbertTransformer : public HilbertTransformer {
public:
	/// the constructor accepts the bandwidth
	SemicircleHilbertTransformer(double t) : t_(t) {}
	///operator() implements abstract virtual operator() of base class HilbertTransformer and performs the actual Hilbert transformation.
	itime_green_function_t operator()(itime_green_function_t& G_tau, FourierTransformer &fourier, double mu, double h, double beta, bool paramagnet)
	{
		
		double tsq=t_*t_;
                matsubara_green_function_t G_omega(G_tau.ntime()-1, G_tau.nsite(), G_tau.nflavor());
                matsubara_green_function_t G0_omega(G_tau.ntime()-1, G_tau.nsite(), G_tau.nflavor());
                itime_green_function_t G0_tau(G_tau.ntime(), G_tau.nsite(), G_tau.nflavor());
		if (paramagnet) {
                  assert(G_tau.nflavor()%2==0);
                  for(spin_t flavor=0;flavor<G_tau.nflavor(); flavor+=2){
                    for(itime_index_t tau=0;tau<G_tau.ntime();++tau){
                      G_tau(tau, flavor  )=0.5*(G_tau(tau, flavor)+G_tau(tau, flavor+1));
                      G_tau(tau, flavor+1)=G_tau(tau, flavor);
                    }
                  }
		}
                fourier.forward_ft(G_tau, G_omega);
		for(uint i=0; i<G_omega.nfreq(); i++) {
		  std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
                  for(spin_t flavor=0;flavor<G_omega.nflavor(); flavor++){
                    G0_omega(i, flavor) = 1./(iw + mu - tsq*G_omega(i, flavor%2==0?flavor+1:flavor-1));
                  }
		}
                fourier.backward_ft(G0_tau, G0_omega);
		return G0_tau;
	}
private:
	double t_; 
};

/// A Hilbert transformation for a square lattice density of states
/// @todo the transformation need to be implemented. responsible: Philipp Werner and Hartmut Monien
class SquarelatticeHilbertTransformer : public HilbertTransformer {
public:
	/// the constructor accepts the bandwidth
	SquarelatticeHilbertTransformer(double bandwidth) : bw_(bandwidth) {}
	itime_green_function_t operator()(itime_green_function_t Sigma, double mu)
	{
		//to implement
          std::cout<<Sigma;
          std::cout<<mu;
		return itime_green_function_t(0,0,0);
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
	itime_green_function_t operator()(itime_green_function_t Sigma, double mu)
	{
		//to implement
          std::cout<<Sigma;
          std::cout<<mu;
		return itime_green_function_t(0,0,0);
	}
private:
	vector_type D_; // or band width
};



// definition of Hilbert transformers living in Fourier (Matsubara) space

/// @brief performs a Hilbert transformation
///
/// The FrequencySpaceHilbertTransformer performs a Hilbert transformation for the self energy and density of states. Arguments are expected to be in Matsubara Frequencies.
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
	virtual matsubara_green_function_t operator()(matsubara_green_function_t & G_omega, matsubara_green_function_t &G0_omega, const double mu, const double h, const double beta,const bool paramagnet)=0; 
	virtual ~FrequencySpaceHilbertTransformer() {}
	
};

/// A Hilbert transformation for a semcircle density of states
/// It receives the dressed Green's function G(\omega) as input and returns tha bare GF G0(\omega)
class FSSemicircleHilbertTransformer : public FrequencySpaceHilbertTransformer {
public:
	/// the constructor accepts the bandwidth
	FSSemicircleHilbertTransformer(double t) : t_(t) {}
	///operator() implements abstract virtual operator() of base class HilbertTransformer and performs the actual Hilbert transformation.
	matsubara_green_function_t operator()(matsubara_green_function_t& G_omega, matsubara_green_function_t &G0_omega_ignored, const double mu, const double h, const double beta, const bool paramagnet)
	{
                std::cout<<"SCHT: using: mu:"<<mu<<" h: "<<h<<"beta: "<<beta<<" paramagnet: "<<paramagnet<<std::endl;
		double tsq=t_*t_;
                matsubara_green_function_t G0_omega(G_omega);
		//formula according to review, p. 61, formula 221. Enforce symmetry up = down; 
                //std::cout<<G_omega<<std::endl;
		if (paramagnet) {
                  if(G_omega.nflavor()==1){ for(uint i=0;i<G_omega.nfreq();++i) G_omega(i,0)=std::complex<double>(0, G_omega(i,0).imag());} //fix for only 1 zone and half filling.
                  if(G_omega.nflavor()>1 && G_omega.nflavor()%2==0){
                    for(uint flavor=0;flavor<G_omega.nflavor();flavor+=2){
                      for(uint i=0;i<G_omega.nfreq();++i){
                        G_omega(i,flavor  )= 0.5*(G_omega(i,flavor)+G_omega(i,flavor+1));
                        G_omega(i,flavor+1)=G_omega(i,flavor);
                      }
                    }
                  }
                }
                if(G_omega.nflavor()==1){ //special case. Paramagnet built in.
                  for(uint i=0; i<G_omega.ntime(); i++) {
                    std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
                    G0_omega(i,0) =1./(iw + mu - tsq*G_omega(i,0));
                  }
                }
                else{
                  assert(G_omega.nflavor()%2==0);
                  for(uint flavor=0;flavor<G_omega.nflavor();flavor+=2){
                    for(uint i=0; i<G_omega.nfreq(); i++) {
                      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
                      G0_omega(i,flavor  ) =1./(iw + mu - tsq*G_omega(i,flavor+1)); 
                      G0_omega(i,flavor+1) =1./(iw + mu - tsq*G_omega(i,flavor  )); 
                    }
                  }
                }
		return G0_omega;
	}
private:
	double t_; 
};
class FSHamiltonianHilbertTransformer: public FrequencySpaceHilbertTransformer{
  public:
    /// the constructor reads in the LDA Hamiltonian from a file
    FSHamiltonianHilbertTransformer(const alps::Parameters& params);
    /// take G_iomegan as the input and construct G0_iomegan out of it. Do this
    /// via summation in k-space. The k-points are defined in the Hamiltonian.
    matsubara_green_function_t operator()(matsubara_green_function_t &G_omega, matsubara_green_function_t &G0_omega, const double mu, const double h, const double beta, const bool paramagnet);
};
class FSDOSHilbertTransformer: public FrequencySpaceHilbertTransformer{
  public:
    /// the constructor reads in the DOS from a file
    FSDOSHilbertTransformer(const alps::Parameters& params){
      assert(params.defined("DOSFILE"));
      std::ifstream dos_file(params["DOSFILE"].c_str());
      assert(dos_file.is_open());
      double eps, d;
      std::cout<<"using density of states: "<<std::endl;
      while(dos_file>>eps>>d){
        epsilon.push_back(eps);
        dos.push_back(d);
        std::cout<<eps<<" "<<d<<std::endl;
      }
      //normalize DOS to one:
      double s=0;
      for(uint i=0;i<epsilon.size()-1;++i){
        s+=dos[i]*(epsilon[i+1]-epsilon[i]); //later: use simpson if this turns out to be bad.
      }
      //std::cout<<"normalization constant: "<<s<<std::endl;
      for(uint i=0;i<epsilon.size()-1;++i){
        dos[i]/=s;
      }
    }
    /// take G_iomegan as the input and construct G0_iomegan out of it. Do this
    /// using integration over D(epsilon).
    /// We perform the HT G <- \int \frac{D(e)}{A - e} de
    matsubara_green_function_t operator()(matsubara_green_function_t &G_omega, matsubara_green_function_t &G0_omega, const double mu, const double h, const double beta, const bool paramagnet){
      assert(G_omega.nsite()==1);
      assert(h==0);
      if(paramagnet && G_omega.nflavor()==1){ for(uint i=0;i<G_omega.nfreq();++i) G_omega(i,0)=std::complex<double>(0, G_omega(i,0).imag());}
      if(paramagnet && G_omega.nflavor()==2){ for(uint i=0;i<G_omega.nfreq();++i){ G_omega(i,0)=0.5*(G_omega(i,0)+G_omega(i,1)); G_omega(i,1)=G_omega(i,0);}}
      //matsubara_green_function_t G0_omega(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
      matsubara_green_function_t G_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
      matsubara_green_function_t G0_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
      matsubara_green_function_t sigma(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
      for(spin_t f=0;f<G_omega.nflavor();++f){
        for(frequency_t w=0;w<G_omega.nfreq();++w){
          std::complex<double> g=0;
          sigma(w,f)= 1./G0_omega(w,f)-1./G_omega(w,f); //self energy
          std::complex<double> A=std::complex<double>(0, (2*w+1)*M_PI/beta) + mu -sigma(w, f);
          for(uint n=0;n<dos.size()-1;++n){
            double delta_epsilon = epsilon[n]-epsilon[n+1];
            double delta_dos     = dos[n]    -dos    [n+1];
            //perform analytic integration of D(e)/(A-e) for D known at points
            //e0, e1 at the end of the interval.
            g+=1./(delta_epsilon)*(delta_dos*delta_epsilon+(A*delta_dos +dos[n+1]*epsilon[n]-dos[n]*epsilon[n+1])*(log(A-epsilon[n]) - log(A-epsilon[n+1])));
          }
          G_omega_new(w, f)=g;
        }
      }
      for(spin_t f=0;f<G_omega.nflavor();++f){
        for(frequency_t w=0;w<G_omega.nfreq();++w){
          G0_omega_new(w,f)=1./(sigma(w,f)+1./G_omega_new(w,f));
        }
      }
      return G0_omega_new;
    }
  private:
    std::vector<double> epsilon;
    std::vector<double> dos;
};
#endif /*ALPS_DMFT_HILBERTTRANSFORMER_H*/
