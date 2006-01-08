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
	virtual itime_multiple_vector_t operator()(itime_multiple_vector_t& G_tau, double mu, double h, double beta, bool paramagnet)=0; 
	virtual ~HilbertTransformer() {}
};

/// A Hilbert transformation for a semcircle density of states
/// It receives G(\tau) as input and returns G0(\tau)
class SemicircleHilbertTransformer : public HilbertTransformer {
public:
	/// the constructor accepts the bandwidth
	SemicircleHilbertTransformer(double t) : t_(t) {}
	///operator() implements abstract virtual operator() of base class HilbertTransformer and performs the actual Hilbert transformation.
	multiple_vector_type operator()(multiple_vector_type& G, double mu, double h, double beta, bool paramagnet)
	{
		
		double tsq=t_*t_;
		
		if (paramagnet) {
			
			vector_type Gpara(G.first.size());
			complex_vector_type G_omega(G.first.size());
			complex_vector_type G0_omega(G.first.size());
			
			for (uint i=0; i<Gpara.size(); i++)
				Gpara[i] = 0.5*(G.first[i]+G.second[i]);
			
			forward_FT(beta, -mu, mu*mu+tsq, Gpara, G_omega);
			
			for(uint i=0; i<G_omega.size(); i++) {
				std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
				G0_omega[i] = 1./(iw + mu - tsq*G_omega[i]);
			}
			
			backward_FT(beta, -mu, mu*mu+tsq, Gpara, G0_omega);
			
			multiple_vector_type G0(G);
			G0.first = Gpara;
			G0.second = Gpara;
			
			return G0;
		}
		else {
			
			double mu_up = mu+h/2.;
			double mu_down = mu-h/2.;
			
			multiple_vector_type G0(G);
			multiple_complex_vector_type G_omega;
			G_omega.first.resize(G.first.size());
			G_omega.second.resize(G.second.size());
			
			multiple_complex_vector_type G0_omega(G_omega); 
			
			forward_FT(beta, -mu_up, mu_up*mu_up+tsq, G.first, G_omega.first);
			forward_FT(beta, -mu_down, mu_down*mu_down+tsq, G.second, G_omega.second);
			
			for(uint i=0; i<G_omega.first.size(); i++) {
				std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
				G0_omega.first[i] = 1./(iw + mu_up - tsq*G_omega.second[i]);
				G0_omega.second[i] = 1./(iw + mu_down - tsq*G_omega.first[i]);
			}
			
			backward_FT(beta, -mu_up, mu_up*mu_up+tsq, G0.first, G0_omega.first);
			backward_FT(beta, -mu_down, mu_down*mu_down+tsq, G0.second, G0_omega.second);
			
			return G0;
		}
		
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

//-------------------------------------------------------------------------------------------------

// Hilbert transformers for the F-functions 

/// A Hilbert transformation for a semcircle density of states
/// It receives G(\tau) as input and returns F(\tau)
class F_SemicircleHilbertTransformer : public HilbertTransformer {
public:
	/// the constructor accepts the bandwidth
	F_SemicircleHilbertTransformer(double t) : t_(t) {}
	///operator() implements abstract virtual operator() of base class HilbertTransformer and performs the actual Hilbert transformation.
	multiple_vector_type operator()(multiple_vector_type& G, double dummy1, double dummy2, double dummy3, bool paramagnet)
	{
	
	    int Np1=G.first.size();
	
		if (paramagnet) {
			
			vector_type Gpara(G.first.size());
			
			for (uint i=0; i<Np1; i++)
				Gpara[i] = 0.5*(G.first[i]+G.second[i]);
						
			multiple_vector_type F(G);
			for (int i=0; i<Np1; i++) {
				F.first[i] = -t_*t_*Gpara[Np1-1-i];		// Green-function must include \tau=\beta !!!
				F.second[i] = -t_*t_*Gpara[Np1-1-i];
			}
			return F;
		}
		else {

			multiple_vector_type F(G);
			for (int i=0; i<Np1; i++) {
				F.first[i] = -t_*t_*G.second[Np1-1-i];
				F.second[i] = -t_*t_*G.first[Np1-1-i];
			}			
			return F;
		}
		
	}
private:
	double t_; 
};

//-------------------------------------------------------------------------------------------------

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
	virtual freq_multiple_vector_t operator()(const freq_multiple_vector_t& G_omega, const double mu, const double h, const double beta,const bool paramagnet)=0; 
	virtual ~FrequencySpaceHilbertTransformer() {}
	
};

/// A Hilbert transformation for a semcircle density of states
/// It receives the dressed Green's function G(\omega) as input and returns tha bare GF G0(\omega)
class FSSemicircleHilbertTransformer : public FrequencySpaceHilbertTransformer {
public:
	/// the constructor accepts the bandwidth
	FSSemicircleHilbertTransformer(double t) : t_(t) {}
	///operator() implements abstract virtual operator() of base class HilbertTransformer and performs the actual Hilbert transformation.
	freq_multiple_vector_t operator()(const freq_multiple_vector_t& G_omega,const double mu, const double h, const double beta, const bool paramagnet)
	{
		
		double tsq=t_*t_;
		//formula according to review, p. 61, formula 221. Enforce symmetry up = down; G(tau)=G(beta-tau).
		if (paramagnet) {
		        freq_multiple_vector_t G0_omega(G_omega);                
                        for(uint i=0; i<G_omega.first.size(); i++) {
                                std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
				//minus sign necessary in order to get reasonable results. check, we might have the wrong convention!
				G0_omega.first[i]=std::complex<double>(0, - 1./(iw.imag() + mu - tsq*0.5*(G_omega.first[i].imag()+G_omega.second[i].imag())));
				G0_omega.second[i]=std::complex<double>(0, - 1./(iw.imag() + mu - tsq*0.5*(G_omega.first[i].imag()+G_omega.second[i].imag())));
                        }
                        return G0_omega;
		}
		else {
		//formula according to review, p. 61, formula 223 'antiferromagnetic Neel phase'	
		//derived properly in formula 95 on page 32. Minus sign still not clear.
			double mu_up = mu+h/2.;
			double mu_down = mu-h/2.;
			
			freq_multiple_vector_t G0_omega(G_omega);
			
			for(uint i=0; i<G_omega.first.size(); i++) {
				std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
				//check: is the zero justified in the AFM case?
				G0_omega.first[i]= std::complex<double>(0,- 1./(iw.imag() + mu_up   - tsq*G_omega.second[i].imag()));
				G0_omega.second[i]= std::complex<double>(0,- 1./(iw.imag() + mu_down - tsq*G_omega.first[i].imag() ));
			}
			return G0_omega;
		}
		
	}
private:
	double t_; 
};

#endif /*ALPS_DMFT_HILBERTTRANSFORMER_H*/
