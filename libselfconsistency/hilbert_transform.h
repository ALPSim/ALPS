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


// definition of Matsubara Hilbert transforms

class FrequencySpaceHilbertTransformer{
public:
  virtual single_site_matsubara_gf operator()(const single_site_matsubara_gf & G_omega,
      single_site_matsubara_gf &G0_omega,
      const double mu, const double h,
      const double beta)=0;
  static void define_parameters(alps::params &p);
  virtual single_site_matsubara_gf initial_G0();
  FrequencySpaceHilbertTransformer(const alps::params &p);
  virtual ~FrequencySpaceHilbertTransformer() {}
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
private:
  bool read_initial_G0(single_site_matsubara_gf &G0_omega);
  void write_initial_G0(const single_site_matsubara_gf &G0_omega);

  //number of Matsubara frequencies
  const int n_matsubara_;
  const int n_orbital_;
  const double beta_;
  const double mu_;
  const double h_init_;
  const std::string input_filename_;
  const bool atomic_initial_G0_;
};



/// A Hilbert transformation for a semcircle density of states
/// It receives the dressed Green's function G(\omega) as input and returns tha bare GF G0(\omega)
class FSSemicircleHilbertTransformer : public FrequencySpaceHilbertTransformer {
public:
  /// the constructor accepts the bandwidth
  FSSemicircleHilbertTransformer(const alps::params &p) : FrequencySpaceHilbertTransformer(p),t_(p["t"]) {}
  ///operator() implements abstract virtual operator() of base class HilbertTransformer 
  ///and performs the actual Hilbert transformation.
  virtual single_site_matsubara_gf operator()(const single_site_matsubara_gf& G_omega,
      single_site_matsubara_gf &G0_omega_ignored,
      const double mu, const double h, const double beta);
private:
  double t_; 
};



/*class FSHamiltonianHilbertTransformer: public FrequencySpaceHilbertTransformer{
public:
  /// the constructor reads in the LDA Hamiltonian from a file
  FSHamiltonianHilbertTransformer(const alps::params& params);
  /// take G_iomegan as the input and construct G0_iomegan out of it. Do this
  /// via summation in k-space. The k-points are defined in the Hamiltonian.
  virtual single_site_matsubara_gf operator()(const single_site_matsubara_gf &G_omega,
      single_site_matsubara_gf &G0_omega,
      const double mu, const double h, const double beta);
};*/



class FSDOSHilbertTransformer: public FrequencySpaceHilbertTransformer{
public:
  /// the constructor reads in the DOS from a file
  FSDOSHilbertTransformer(const alps::params& params);
  /// take G_iomegan as the input and construct G0_iomegan out of it. Do this
  /// using integration over D(epsilon).
  /// We perform the HT G <- \int \frac{D(e)}{A - e} de
  virtual single_site_matsubara_gf operator()(const single_site_matsubara_gf &G_omega,
      single_site_matsubara_gf &G0_omega,
      const double mu, const double h, const double beta);

protected:
  std::vector<double> epsilon;
  std::vector<double> dos;
};



inline std::complex<double> integrand(int n, const std::vector<double> &dos, const std::vector<double> &epsilon, 
    const std::complex<double> &zeta_A,
    const std::complex<double> &zeta_B)
    {
  return dos[n]/(zeta_A*zeta_B - epsilon[n]*epsilon[n]); 
    }



//this can do Antiferromagnetism according to Georges et. al.
class AFM_FSDOSHilbertTransformer: public FSDOSHilbertTransformer{
public:
  AFM_FSDOSHilbertTransformer(const alps::params &parms)
:FSDOSHilbertTransformer(parms){}
  virtual ~AFM_FSDOSHilbertTransformer(){}
  virtual single_site_matsubara_gf operator()(const single_site_matsubara_gf &G_omega,
      single_site_matsubara_gf &G0_omega,
      const double mu, const double h, const double beta);
};


class TwoDAFMHilbertTransformer: public FrequencySpaceHilbertTransformer{
public:
  TwoDAFMHilbertTransformer(const alps::params& params):FrequencySpaceHilbertTransformer(params),
    t_(params["t"]),
    tprime_(params["tprime"]),
    L_(params["L"]){
    if(L_%2 !=0) throw std::runtime_error("please choose L to be even for Simpson integration for fast convergence!!");
  }


  single_site_matsubara_gf operator()(const single_site_matsubara_gf & G_omega,
      single_site_matsubara_gf &G0_omega,
      const double mu, const double h, const double beta);
  class integrand{
  public:
    virtual std::complex<double>operator()(double kx, double ky) const=0;
    virtual ~integrand(){};
  };
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
};
