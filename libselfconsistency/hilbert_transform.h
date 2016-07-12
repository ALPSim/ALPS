/*****************************************************************************
 *
 * Copyright (C) 2005-2015 Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#pragma once
#include "alps/params.hpp"
#include "alps/gf/tail.hpp"
typedef alps::gf::omega_sigma_gf_with_tail single_site_matsubara_gf;
typedef alps::gf::itime_sigma_gf_with_tail single_site_itime_gf;


///Matsubara space Hilbert transforms
class FrequencySpaceHilbertTransformer{
public:
  ///define the parameters that this class uses
  static void define_parameters(alps::params &p);
  ///initialize a single site Hilbert transformer from a parameter file
  FrequencySpaceHilbertTransformer(const alps::params &p);
  ///run the Hilbert transform
  virtual single_site_matsubara_gf operator()(const single_site_matsubara_gf & G_omega, single_site_matsubara_gf &G0_omega)=0;
  ///Create an initial Green's function for starting a simulation
  virtual single_site_matsubara_gf initial_G0();
  ///virtual destructor
  virtual ~FrequencySpaceHilbertTransformer() {}
  ///function to extract the average of the dispersion, needed for the high frequency terms
  virtual double epsilon_av(int orbital) const=0;
  ///function to extract the average of the square of the dispersion, needed for the high frequency terms
  virtual double epsilonsq_av(int orbital) const=0;

  ///symmetrization function for performing a spin symmetrization
  template<class T, class MESH1> void spin_symmetrize(alps::gf::two_index_gf<T, MESH1, alps::gf::index_mesh> &GF) const
  {
    if(GF.mesh2().extent()%2!=0) throw std::invalid_argument("for spin symmetrization you have to have an even number of orbitals.");
    for(typename MESH1::index_type i(0); i<GF.mesh1().extent();++i){
      for(alps::gf::index_mesh::index_type flavor(0);flavor<GF.mesh2().extent(); flavor+=2){
        T symmetrized_value=0.5*(GF(i,flavor)+GF(i,flavor+1));
        GF(i,flavor  )=symmetrized_value;
        GF(i,flavor+1)=symmetrized_value;
      }
    }
  }
protected:
  ///read an input G0 file
  bool read_initial_G0(single_site_matsubara_gf &G0_omega);
  ///write a newly generated input G0 file so that it can be used as an input next time.
  void write_initial_G0(const single_site_matsubara_gf &G0_omega);
  ///set the three highfrequency tails
  void set_tails(single_site_matsubara_gf &G0_omega) const{ set_c1_unity(G0_omega); set_c2(G0_omega); set_c3(G0_omega);}
  ///For Green's functions, the c1 high frequency term is 1
  void set_c1_unity(single_site_matsubara_gf &G0_omega) const;
  ///For Green's functions, the c2 high frequency term is mu - epsilon +/- h
  void set_c2(single_site_matsubara_gf &G0_omega) const;
  ///For Green's functions, the c2 high frequency term is (mu +/- h)^2 -2 (mu +/- h) *epsilon + epsilon^2
  void set_c3(single_site_matsubara_gf &G0_omega) const;

  //number of Matsubara frequencies
  const int n_matsubara_;
  const int n_orbital_;
  const double beta_;
  const double mu_;
  const std::string input_filename_;
  const bool atomic_initial_G0_;
};



/// A Hilbert transform for a semcircle density of states
class FSSemicircleHilbertTransformer : public FrequencySpaceHilbertTransformer {
public:
  /// the constructor accepts the bandwidth
  FSSemicircleHilbertTransformer(const alps::params &p) : FrequencySpaceHilbertTransformer(p),t_(p["t"]) {}
  ///function to extract the average of the dispersion, needed for the high frequency terms
  double epsilon_av(int orbital) const{(void) orbital; return 0.;}
  ///function to extract the average of the square of the dispersion, needed for the high frequency terms
  double epsilonsq_av(int orbital) const{(void) orbital; return t_*t_;}
  ///operator() implements abstract virtual operator() of base class HilbertTransformer and performs the actual Hilbert transform.
  virtual single_site_matsubara_gf operator()(const single_site_matsubara_gf & G_omega, single_site_matsubara_gf &G0_omega);
  ///Create an initial Green's function for starting a simulation
  single_site_matsubara_gf initial_G0();

private:
  double t_; 
};

/// A Hilbert transform for a general density of states
class FSDOSHilbertTransformer: public FrequencySpaceHilbertTransformer{
public:
  /// the constructor reads in the DOS from a file
  FSDOSHilbertTransformer(const alps::params& params);
  ///function to extract the average of the dispersion, needed for the high frequency terms
  double epsilon_av(int i) const{(void) i; return eps_av_;}
  ///function to extract the average of the square of the dispersion, needed for the high frequency terms
  double epsilonsq_av(int i) const{(void) i; return epssq_av_;}
  ///operator() implements abstract virtual operator() of base class HilbertTransformer and performs the actual Hilbert transform.
  virtual single_site_matsubara_gf operator()(const single_site_matsubara_gf & G_omega, single_site_matsubara_gf &G0_omega);

protected:
  ///make sure that the DOS is normalized to 1.
  void normalize_dos();
  ///find <epsilon> and <epsilon^2> via integration over the BZ.
  void compute_moments();

  std::vector<double> epsilon_;
  std::vector<double> dos_;
  double eps_av_;
  double epssq_av_;
};

/// A Hilbert transform for a two-dimensional band structure. Can do both PM and AFM solutions
class TwoDAFMHilbertTransformer: public FrequencySpaceHilbertTransformer{
public:
  ///constructor: set up the 2d AFM integration
  TwoDAFMHilbertTransformer(const alps::params& params);
  ///function to extract the average of the dispersion, needed for the high frequency terms
  double epsilon_av(int i) const{(void) i; return 0.;}
  ///function to extract the average of the square of the dispersion, needed for the high frequency terms
  double epsilonsq_av(int i) const{(void) i; return 4*(t_*t_+tprime_*tprime_);}
  ///operator() implements abstract virtual operator() of base class HilbertTransformer and performs the actual Hilbert transform.
  virtual single_site_matsubara_gf operator()(const single_site_matsubara_gf & G_omega, single_site_matsubara_gf &G0_omega);

  ///class for doing integration in 2d
  class integrand{
  public:
    virtual std::complex<double>operator()(double kx, double ky) const=0;
    virtual ~integrand(){};
  };
  ///class implementing the integral that needs to be done: see equation 97 in the Georges review. We have the 2d dispersion hard coded here.
  class afm_integrand:public integrand{ //class for the integrand with an epsilon(k) that uses the t-band structure.
  public:
    afm_integrand(double t, double tprime, std::complex<double> zeta_A, std::complex<double> zeta_B){
      zeta_A_=zeta_A;
      zeta_B_=zeta_B;
      t_=t;
      tprime_=tprime;
    }
    std::complex<double> operator()(double kx, double ky) const{
      double eps=twod_square_lattice_dispersion(kx,ky);
      //std::cout<<"disp: "<<eps<<" t: "<< t_<<" tprime: "<<tprime_<<" zeta A: "<<zeta_A_<<" zeta B: "<<zeta_B_<<" kx: "<<kx<<" ky "<<ky<<std::endl;
      return 1./(zeta_A_*zeta_B_-eps*eps);
    }
  private:
    double twod_square_lattice_dispersion(double kx, double ky) const{
      return -2.*t_*(cos(kx)+cos(ky)) -4.*tprime_*cos(kx)*cos(ky);
    }
    double t_;
    double tprime_;
    std::complex<double> zeta_A_;
    std::complex<double> zeta_B_;
  };

private:
  std::complex<double> twodsimpson(const integrand &, double ax, double ay, double bx, double by, int N);
  const double t_;
  const double tprime_;
  const int L_;
  bool paramagnet_;

};
