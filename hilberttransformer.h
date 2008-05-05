/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005-2007 by Philipp Werner <werner@comp-phys.org>,
*                         Matthias Troyer <troyer@comp-phys.org>,
*                         Emanuel Gull <gullc@comp-phys.org>,
*                         Sebastian Fuchs <fuchs@comp-phys.org>
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
/// @sa HilbertTransformer, SemicircleHilbertTransformer, SquarelatticeHilbertTransformer, 
///     GeneralHilbertTransformer, DCATransformer
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
  virtual itime_green_function_t operator()(const itime_green_function_t& G_tau, FourierTransformer &fourier_ignored, 
                                            double mu, double h, double beta)=0; 
  

  itime_green_function_t symmetrize(const itime_green_function_t& G_tau, const bool paramagnet) const {
    itime_green_function_t G(G_tau);
    if (paramagnet) {
      assert(G.nflavor()%2==0);
      for(spin_t flavor=0;flavor<G.nflavor(); flavor+=2){
        for(itime_index_t tau=0;tau<G.ntime();++tau){
          G(tau, flavor  )=0.5*(G(tau, flavor)+G(tau, flavor+1));
          G(tau, flavor+1)=G(tau, flavor);
        }
      }
    }
    return G;
  }

  virtual itime_green_function_t initial_G0(const alps::Parameters& parms) {
    std::cout<<"reading initial G0_tau"<<std::endl;
    int n_tau=boost::lexical_cast<int>(parms["N"]);
    int n_orbital=parms.value_or_default("FLAVORS", 2);
    itime_green_function_t G0_tau(n_tau+1, 1, n_orbital);
    G0_tau.read(parms["G0TAU_INPUT"].c_str());
    return G0_tau;
  }
  
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
  itime_green_function_t operator()(const itime_green_function_t& G_tau, FourierTransformer &fourier, 
                                    double mu, double h, double beta)
  {
    if(h!=0){ std::cerr<<"ignoring magnetic field!"<<std::endl; abort(); }
    double tsq=t_*t_;
    matsubara_green_function_t G_omega(G_tau.ntime()-1, G_tau.nsite(), G_tau.nflavor());
    matsubara_green_function_t G0_omega(G_tau.ntime()-1, G_tau.nsite(), G_tau.nflavor());
    itime_green_function_t G0_tau(G_tau.ntime(), G_tau.nsite(), G_tau.nflavor());
    fourier.forward_ft(G_tau, G_omega);
    //std::cout<<"G omega real: "<<std::endl;
    print_real_green_matsubara(std::cout, G_omega, beta);
    for(uint i=0; i<G_omega.nfreq(); i++) {
      std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
      for(spin_t flavor=0;flavor<G_omega.nflavor(); flavor++){
        G0_omega(i, flavor) = 1./(iw + mu - tsq*G_omega(i, flavor%2==0?flavor+1:flavor-1));
      }
    }
    //std::cout<<"symmetrized G0 omega real: "<<std::endl;
    print_real_green_matsubara(std::cout, G0_omega, beta);
    fourier.backward_ft(G0_tau, G0_omega);
    return G0_tau;
  }
  
  itime_green_function_t initial_G0(const alps::Parameters& parms) 
  {
    if (parms.defined("G0TAU_INPUT"))
      return HilbertTransformer::initial_G0(parms);
    else {
      std::cout<<"reading initial G0_tau"<<std::endl;
      boost::shared_ptr<FourierTransformer> fourier_ptr;
      FourierTransformer::generate_transformer(parms, fourier_ptr);
      int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
      int n_time=boost::lexical_cast<int>(parms["N"]);
      int n_flavor=parms.value_or_default("FLAVORS", 2);
      double beta = static_cast<double>(parms["BETA"]);
      double mu = static_cast<double>(parms["MU"]);
      double h = static_cast<double>(parms["H"]);
      if(h!=0) {std::cerr<<"we don't yet handle magnetic fields!"<<std::endl; abort(); }
      double h_init = parms.value_or_default("H_INIT", 0.);
      double tsq=t_*t_;
      matsubara_green_function_t G0_omega(n_matsubara,n_flavor);  
      itime_green_function_t G0_tau(n_time, n_flavor);
      for(uint i=0; i<G0_omega.nfreq(); i++) {
        std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
        for(spin_t flavor=0;flavor<G0_omega.nflavor(); flavor++) {
          std::complex<double> zeta = iw+mu+(flavor%2 ? -h_init : h_init);
          G0_omega(i, flavor) = (zeta - sqrt(zeta*zeta-4*tsq))/(2*tsq);
        }
      }
      G0_tau.write("G0_tau_input");
      G0_omega.write("G0_omega_input");
      return G0_tau;
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
  virtual itime_green_function_t operator()(const itime_green_function_t& G_tau, FourierTransformer &fourier_ignored,
                                                    double mu, double h, double beta)
  {
    //to implement
    std::cerr<<"Square Lattice Hilbert transformer"<<G_tau<<std::endl;
    std::cout<<&fourier_ignored<<std::endl;
    std::cout<<mu;
    std::cout<<beta;
    std::cout<<h;
    abort();
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
  virtual itime_green_function_t operator()(const itime_green_function_t& G_tau, FourierTransformer &fourier_ignored,
                                                  double mu, double h, double beta)
  {
    //to implement
    //std::cout<<Sigma;
    std::cerr<<G_tau<<std::endl;
    std::cerr<<&fourier_ignored<<std::endl;
    std::cerr<<mu;
    std::cerr<<beta;
    std::cerr<<h;
    abort();
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
  virtual matsubara_green_function_t initial_G0(const alps::Parameters& parms) {
    std::cout<<"reading initial G0_omega"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    int n_orbital=parms.value_or_default("FLAVORS", 2);
    matsubara_green_function_t G0_omega(n_matsubara, n_orbital);
    G0_omega.read(parms["G0OMEGA_INPUT"].c_str());
    return G0_omega;
  }
  
  virtual ~FrequencySpaceHilbertTransformer() {}
  template <class T>
  green_function<T> symmetrize(const green_function<T>& G, const bool paramagnet) const
  {
    green_function<T> G_new(G);
    if (paramagnet) {
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
                                                const double mu, const double h, const double beta)
  {
    std::cout<<"SCHT: using: mu:"<<mu<<" h: "<<h<<"beta: "<<beta<<" "<<G0_omega_ignored.nfreq()<<std::endl;
    double tsq=t_*t_;
    matsubara_green_function_t G0_omega(G_omega);
    //formula according to review, p. 61, formula 221. 
    if(G_omega.nflavor()==1){ //special case. 
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
  
  matsubara_green_function_t initial_G0(const alps::Parameters& parms)
  {
    if (parms.defined("G0OMEGA_INPUT"))
      return FrequencySpaceHilbertTransformer::initial_G0(parms);
    else {
      std::cout<<"calculating initial G0_omega"<<std::endl;
      int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
      int n_time=boost::lexical_cast<int>(parms["N"]);
      int n_flavor=parms.value_or_default("FLAVORS", 2);
      double beta = static_cast<double>(parms["BETA"]);
      double mu = static_cast<double>(parms["MU"]);
      double h = static_cast<double>(parms["H"]);
      double h_init = parms.value_or_default("H_INIT", 0.);
      double tsq=t_*t_;
      matsubara_green_function_t G0_omega(n_matsubara,n_flavor);
      for(uint i=0; i<G0_omega.nfreq(); i++) {
        std::complex<double> iw(0.,(2*i+1)*M_PI/beta);
        for(spin_t flavor=0;flavor<G0_omega.nflavor(); flavor++) {
          std::complex<double> zeta = iw+mu+(flavor%2 ? -h_init-h : h_init-h);
          std::complex<double> tmp = sqrt(zeta*zeta-4*tsq);
          tmp *= tmp.imag()<0 ? -1 : 1;
          G0_omega(i, flavor) = (zeta - tmp)/(2*tsq);
          //std::cout << zeta << " " << sqrt(zeta*zeta) << std::endl;
        }
      }
      G0_omega.write("G0_omega_input");
      boost::shared_ptr<FourierTransformer> fourier_ptr;
      FourierTransformer::generate_transformer(parms, fourier_ptr);
      itime_green_function_t G0_tau(n_time,n_flavor);
      fourier_ptr->backward_ft(G0_tau, G0_omega);
      G0_tau.write("G0_tau_input");
      return G0_omega;
    }
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
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega, 
                                                matsubara_green_function_t &G0_omega, 
                                                const double mu, const double h, const double beta);
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
    if(dos.size()%2!=0){ std::cerr<<"please use even number of DOS points"<<std::endl; abort(); }
    //normalize DOS to one:
    double s=0;
    for(uint i=1;i<epsilon.size()-2;i+=2){
      s+=4.*dos[i]+2*dos[i+1];
    }
    s+=dos[0]+dos[epsilon.size()-1]+4.*dos[epsilon.size()-2];
    s/=(3.);
    std::cout<<"normalization constant: "<<s<<std::endl;
    std::cout<<"h is: "<<(epsilon[1]-epsilon[0])<<std::endl;
    for(uint i=0;i<epsilon.size()-1;++i){
      dos[i]/=s;
    }
    double S=0;
    for(uint i=0;i<epsilon.size()-1;++i){
      S+=dos[i];
    }
    std::cout<<"check: total sum is: "<<S<<std::endl;
  }
  /// take G_iomegan as the input and construct G0_iomegan out of it. Do this
  /// using integration over D(epsilon).
  /// We perform the HT G <- \int \frac{D(e)}{A - e} de
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega, 
                                                matsubara_green_function_t &G0_omega, 
                                                const double mu, const double h, const double beta)
  {
    assert(G_omega.nsite()==1);
    if(h!=0){std::cerr<<" we don't yet handle magnetic fields!"<<std::endl; abort(); }
    //matsubara_green_function_t G0_omega(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
    matsubara_green_function_t G_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
    matsubara_green_function_t G0_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
    matsubara_green_function_t sigma(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
    for(spin_t f=0;f<G_omega.nflavor();++f){
      //std::cout<<"check of G0omega: "<<f<<" "<<G0_omega(0,f)<<" "<<G_omega(0,f);
      for(frequency_t w=0;w<G_omega.nfreq();++w){
        std::complex<double> g=0;
        sigma(w,f)= 1./G0_omega(w,f)-1./G_omega(w,f); 
        std::complex<double> A=std::complex<double>(0, (2*w+1)*M_PI/beta) + mu -sigma(w, f);
        for(uint n=0;n<dos.size()-1;++n){
          //double delta_epsilon = epsilon[n]-epsilon[n+1];
          //double delta_dos     = dos[n]    -dos    [n+1];
          //perform analytic integration of D(e)/(A-e) for D known at points
          //e0, e1 at the end of the interval.
          //g+=1./(delta_epsilon)*(delta_dos*delta_epsilon+(A*delta_dos +dos[n+1]*epsilon[n]-dos[n]*epsilon[n+1])
          //                       *(log(A-epsilon[n]) - log(A-epsilon[n+1])));
          g+=0.5*(dos[n+1]-dos[n])*(log(A-epsilon[n]) + log(A-epsilon[n+1]));
          //g-=delta_epsilon*dos[n]/(A-epsilon[n]);
        }
        G_omega_new(w, f)=g;
      }
      //std::cout<<"G omega new is: "<<G_omega_new(0,f)<<std::endl;
    }
    if(G_omega.nflavor() !=2){std::cerr<<"don't know how to handle != 2 flavors."<<std::endl; abort();}
    for(spin_t f=0;f<G_omega.nflavor();++f){
      for(frequency_t w=0;w<G_omega.nfreq();++w){
        G0_omega_new(w,f)=1./(sigma(w,f)+1./G_omega_new(w,f));
      }
    }
    //std::cout<<"G0_omega_new"<<" 0 "<<G0_omega_new(0,0)<<std::endl;
    //std::cout<<"G0_omega_new"<<" 1 "<<G0_omega_new(0,1)<<std::endl;
    return G0_omega_new;
  }
  
protected:
  std::vector<double> epsilon;
  std::vector<double> dos;
};
inline std::complex<double> integrand(int n, const std::vector<double> &dos, const std::vector<double> &epsilon, const std::complex<double> &zeta_A,
    const std::complex<double> &zeta_B){ return dos[n]/(zeta_A*zeta_B - epsilon[n]*epsilon[n]); }
//this can do Antiferromagnetism according to Georges et. al.
class AFM_FSDOSHilbertTransformer: public FSDOSHilbertTransformer{
  public:
  AFM_FSDOSHilbertTransformer(const alps::Parameters &parms)
  :FSDOSHilbertTransformer(parms){}
  virtual ~AFM_FSDOSHilbertTransformer(){}
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,matsubara_green_function_t &G0_omega, 
      const double mu, const double h, const double beta)
  {
    if(G_omega.nsite()!=1){ std::cerr<<"don't know how to handle n_site >1. Use Sebastian's cluster loop for that. aborting."<<std::endl; abort();}
    if(G_omega.nflavor()!=2){ std::cerr<<"don't know how to handle n_flavor!=2. aborting."<<std::endl; abort();}
    matsubara_green_function_t Sigma(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
    matsubara_green_function_t G_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
    matsubara_green_function_t G0_omega_new(G_omega.nfreq(), G_omega.nsite(), G_omega.nflavor());
    //compute self energy
    for(frequency_t w=0;w<G_omega.nfreq();++w){
      double wn=(2.*w+1)*M_PI/beta;
      for(spin_t f=0;f<G_omega.nflavor();++f){
        Sigma(w,f)=1./G0_omega(w,f)-1./G_omega(w,f);
      }
      //std::cout<<"sigma: "<<Sigma(w,0)<<" "<<Sigma(w,1)<<std::endl;
      std::complex<double>zeta_A=std::complex<double>(mu-h,wn)-Sigma(w,0);
      std::complex<double>zeta_B=std::complex<double>(mu+h,wn)-Sigma(w,1);

      //Simpson: integrate dos(e)/(zeta_A zeta_B-e^2)
      if(dos.size()%2 !=0){std::cerr<<"for Simpson precision: use a DOS with even number of integration points. "<<std::endl; abort(); }
      std::complex<double> I=0;
      for(int i=1;i<dos.size()-2;i+=2){
        I+=4.*integrand(i,dos,epsilon,zeta_A,zeta_B)+2.*integrand(i+1,dos,epsilon,zeta_A,zeta_B);
      }
      I+=integrand(0,dos,epsilon,zeta_A,zeta_B)+integrand(dos.size()-1,dos,epsilon,zeta_A,zeta_B)+4.*integrand(dos.size()-2,dos,epsilon,zeta_A,zeta_B);
      I/=3.;
      /*for(int i=0;i<dos.size();++i){
        I+=dos[i]/(zeta_A*zeta_B - epsilon[i]*epsilon[i]);
      }*/
      //compute the new G's:
      G_omega_new(w,0)=zeta_B*I; //formula 97 in review
      G_omega_new(w,1)=zeta_A*I;
      G0_omega_new(w,0)=1./(1./G_omega_new(w,0)+Sigma(w,0));
      G0_omega_new(w,1)=1./(1./G_omega_new(w,1)+Sigma(w,1));
      //std::cout<<"G0: "<<G0_omega_new(w,0)<<" "<<G0_omega_new(w,1)<<std::endl;
    }
    //std::cout<<std::endl<<std::endl;
    std::cout<<"sigma: "<<std::endl;
    for(int i=0;i<Sigma.nfreq();++i){
      std::cout<<(2.*i+1)*M_PI/beta<<" "<<Sigma(i,0).real()<<" "<<Sigma(i,0).imag()<<" "<<Sigma(i,1).real()<<" "<<Sigma(i,1).imag()<<std::endl;
    }
    std::cout<<"g0 omega new: "<<std::endl;
    for(int i=0;i<G0_omega_new.nfreq();++i){
      std::cout<<(2.*i+1)*M_PI/beta<<" "<<G0_omega_new(i,0).real()<<" "<<G0_omega_new(i,0).imag()<<" "<<G0_omega_new(i,1).real()<<" "<<G0_omega_new(i,1).imag()<<std::endl;
    }
    return G0_omega_new;
  }
};
#endif /*ALPS_DMFT_HILBERTTRANSFORMER_H*/
