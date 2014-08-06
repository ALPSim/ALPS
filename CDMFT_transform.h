#ifndef CDMFT_TRANSFORM_H
#define CDMFT_TRANSFORM_H

#include"cluster.h"
#include "boost/numeric/ublas/io.hpp"
typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major> dense_matrix;
typedef boost::numeric::ublas::matrix<std::complex<double>,boost::numeric::ublas::column_major> complex_dense_matrix;

//The CDMFT transormer does the Fourier transforms, dispersions, etc that are needed for the CDMFT.
class CDMFTTransformer : public ClusterTransformer
{
public:
  CDMFTTransformer(const alps::params& params);
  virtual ~CDMFTTransformer(){} //pro forma virtual destructor
  
  //symmetrization in real momentum space - for the normal GF
  template<typename T> green_function<T> symmetrize_real_space(const green_function<T>& G) const
  { return symmetrize_helper(G, pair_classes_real_space()); }
  template<typename T> green_function<T> symmetrize_k_space(const green_function<T>& G) const
  { return symmetrize_helper(G, pair_classes_k_space()); }
  
  //symmetrization in real and momentum space - for the NAMBU GF
  template<typename T> nambu_green_function<T> symmetrize_real_space(const nambu_green_function<T>& G) const
  { return symmetrize_helper(G, pair_classes_real_space()); }
  template<typename T> nambu_green_function<T> symmetrize_k_space(const nambu_green_function<T>& G) const
  { return symmetrize_helper(G, pair_classes_k_space()); }
  
  
  //Real space and k-space transformers for the Green's functions
  virtual matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const=0;
  virtual matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const=0;
  virtual itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const=0;
  virtual itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const=0;
  //Real space and k-space transformers for the Nambu Green's function
  virtual nambu_matsubara_green_function_t transform_into_real_space(const nambu_matsubara_green_function_t& G) const=0;
  virtual nambu_matsubara_green_function_t transform_into_k_space(const nambu_matsubara_green_function_t& G) const=0;
  virtual nambu_itime_green_function_t transform_into_real_space(const nambu_itime_green_function_t& G) const=0;
  virtual nambu_itime_green_function_t transform_into_k_space(const nambu_itime_green_function_t& G) const=0;
  
  //Hilbert transform operator
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                                const matsubara_green_function_t &G0_omega,
                                                double &mu, const double h, const double beta)=0;
  //Hilbert transform operator for Nambu GF
  virtual nambu_matsubara_green_function_t operator()(const nambu_matsubara_green_function_t &G_omega,
                                                      const nambu_matsubara_green_function_t &G0_omega,
                                                      nambu_matsubara_green_function_t &selfenergy,
                                                      double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy)=0;
  
  int* momentum_translation() const;
  double dispersion_moments(const unsigned int i) const { return dispersion_moments_[i]; }
  double bandwidth() const { return bandwidth_; };
  
  double epsilon(const vector_type k) const
  {
    throw std::logic_error("CDMFT dispersion not implemented.");
    return 0.;
  }
  virtual std::vector<double> epsav() const{throw std::logic_error("CDMFT dispersion not implemented.");}
  virtual std::vector<double> epssqav() const {throw std::logic_error("CDMFT dispersion not implemented.");};

protected:
  complex_dense_matrix U_;
  
private:
  template<typename T>
  green_function<T> symmetrize_helper(const green_function<T>& G,
                                      const std::pair<std::set<pair_class_type>::const_iterator,
                                      std::set<pair_class_type>::const_iterator> pair_classes) const;
  template<typename T>
  nambu_green_function<T> symmetrize_helper(const nambu_green_function<T>& G,
                                            const std::pair<std::set<pair_class_type>::const_iterator,
                                            std::set<pair_class_type>::const_iterator> pair_classes) const;
  
  
  
  const double t_, t2_, t3_, t4_;
  double bandwidth_;
  const int L_;
  double bandwidth_scale_factor_;
  double half_filling_;
  std::vector<double> dispersion_moments_;
};




//does the CDMFT cluster transform for the paramagnetic phase (no AFM order)
class LadderCDMFTTransformer : public CDMFTTransformer
{
public:
  LadderCDMFTTransformer(const alps::params& params) :
  CDMFTTransformer(params),
  compute_filling_(params.defined("FILLING")),
  n_target_(compute_filling_?(double)(params["FILLING"]):-1),
  epsav_(n_cluster_sites(),0.),
  epssqav_(n_cluster_sites(),0.),
  L_(params["L"]),
  t_(params["t"]),
  tprime_(params["tprime"]|0),
  tstar_(params["tstar"]),
  params_(params)
  {
    set_U_(n_cluster_sites(),U_);

    /*U_.resize(4,4);
    U_(0,0)=1.; U_(0,1)=1.; U_(0,2)=1.; U_(0,3)=1.;
    U_(1,0)=1.; U_(1,1)=-1.; U_(1,2)=1.; U_(1,3)=-1.;
    U_(2,0)=1.; U_(2,1)=1.; U_(2,2)=-1.; U_(2,3)=-1.;
    U_(3,0)=1.; U_(3,1)=-1.; U_(3,2)=-1.; U_(3,3)=1.;
    U_/=2.;*/
     
     /*epsav_[0]=-8*(M_PI*t+2*tprime)/(M_PI*M_PI);
     epsav_[1]= 8*(M_PI*t-2*tprime)/(M_PI*M_PI);
     epsav_[2]=16*tprime/(M_PI*M_PI);
     epsav_[3]=16*tprime/(M_PI*M_PI);
     
     epssqav_[0]=4*((8+M_PI*M_PI)*t*t+8*M_PI*t*tprime+M_PI*M_PI*tprime*tprime)/(M_PI*M_PI);
     epssqav_[1]=4*((8+M_PI*M_PI)*t*t-8*M_PI*t*tprime+M_PI*M_PI*tprime*tprime)/(M_PI*M_PI);
     epssqav_[2]=4*((1-8/(M_PI*M_PI))*t*t+tprime*tprime);
     epssqav_[3]=4*((1-8/(M_PI*M_PI))*t*t+tprime*tprime);*/
  }
  //Hilbert transform operator
  matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega,
                                        const matsubara_green_function_t &G0_omega,
                                        double &mu, const double h, const double beta);
  //Hilbert transform operator for Nambu GF
  nambu_matsubara_green_function_t operator()(const nambu_matsubara_green_function_t &G_omega,
                                              const nambu_matsubara_green_function_t &G0_omega,
                                              nambu_matsubara_green_function_t &selfenergy,
                                              double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy);
  
  //spin symmetrization functions
  itime_green_function_t symmetrize(const itime_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t symmetrize(const matsubara_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }
  //Fourier transform functions
  matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const;
  matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const ;
  itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const ;
  itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const ;
  
  //spin symmetrization functions in Nambu space
  nambu_itime_green_function_t symmetrize(const nambu_itime_green_function_t& G, const bool paramagnet) const
  {
    return spin_symmetrizer(G, paramagnet);
  }
  nambu_matsubara_green_function_t symmetrize(const nambu_matsubara_green_function_t& G, const bool paramagnet) const
  {
    return spin_symmetrizer(G, paramagnet);
  }
  nambu_matsubara_green_function_t transform_into_real_space(const nambu_matsubara_green_function_t& G) const ;
  nambu_matsubara_green_function_t transform_into_k_space(const nambu_matsubara_green_function_t& G) const;
  nambu_itime_green_function_t transform_into_real_space(const nambu_itime_green_function_t& G) const ;
  nambu_itime_green_function_t transform_into_k_space(const nambu_itime_green_function_t& G) const ;
  
  //initial G0
  virtual matsubara_green_function_t initial_G0      (const alps::params& parms, shape_t shape=nondiagonal);
  virtual nambu_matsubara_green_function_t initial_G0_nambu(const alps::params& parms);
  std::vector<double> epsav() const;
  std::vector<double> epssqav() const;
private:
  template<typename T, class D>
  green_function<T> transform_coordinates(const       green_function<T>& G) const;
  template<typename T, class D>
  nambu_green_function<T> transform_coordinates(const nambu_green_function<T>& G) const;
  template <class T>
  green_function<T> spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const;
  
  //in nambu case: <double> very different from <complex<double> > -> no template
  nambu_green_function<std::complex<double> > spin_symmetrizer(const nambu_green_function<std::complex<double> >& G, const bool paramagnet) const;
  nambu_green_function<double> spin_symmetrizer(const nambu_green_function<double>& G, const bool paramagnet)const {throw std::logic_error("not implemented for obsevaluators. ");}
  
  matsubara_green_function_t fourpoint_green_function_inverse(const matsubara_green_function_t& G, const unsigned int nt) const;
  template<typename T>
  green_function<T> fourpoint_green_function_add(const green_function<T>& G1,
                                                 const green_function<T>& G2) const;
  template<typename T>
  green_function<T> fourpoint_green_function_subtract(const green_function<T>& G1,
                                                      const green_function<T>& G2) const;
  template<typename T>
  green_function<T> fourpoint_green_function_divide(const green_function<T>& G1,
                                                    const T d) const;
  void set_Tij_matrix(double kx, double ky, complex_dense_matrix &H0);
  void set_nambu_Tij_dwave_matrix(double kx, double ky, double eta1,double eta2,double mu, double h,complex_dense_matrix &H0);
  void set_nambu_Tij_dwave_matrix(double kx, double ky, double eta1,double eta2,double mu, complex_dense_matrix &H0);
  void set_U_(unsigned int nsite, complex_dense_matrix &U_);
  struct Forward;
  struct Backward;
  
  bool compute_filling_;
  double n_target_;
  //complex_dense_matrix U_;
  std::vector<double> epsav_;
  std::vector<double> epssqav_;
  const int L_;
  const double t_;
  const double tprime_;
  const double tstar_;

  const alps::params params_;
};

#endif /*ALPS_CLUSTER_H*/
