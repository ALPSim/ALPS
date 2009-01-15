/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2007 by Sebastian Fuchs <fuchs@comp-phys.org>
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

#ifndef ALPS_CLUSTER_H
#define ALPS_CLUSTER_H

#include "types.h"
#include "green_function.h"
#include "hilberttransformer.h"
#include <alps/lattice/graph_helper.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/lapack/gesv.hpp>

extern "C" void dgtsv_(const int*, const int*, double*, double*, double*, double*, const int*, int*);
extern "C" void zgtsv_(const int*, const int*, void*, void*, void*, void*, const int*, int*);


class ClusterHelper: private alps::graph_helper<alps::graph_helper<>::graph_type>
{
public:
  
  template <typename T> struct Matrix
  { typedef boost::numeric::ublas::matrix<T, boost::numeric::ublas::column_major> Type; };
  typedef alps::graph_helper<alps::graph_helper<>::graph_type>::vector_type vector_type;
  typedef std::set<site_t> site_class_type;
  typedef std::pair<site_t, site_t> pair_type;
  typedef std::set<pair_type> pair_class_type;
  typedef std::set<site_class_type> site_class_set_type;
  typedef std::set<pair_class_type> pair_class_set_type;
  typedef std::pair<site_class_set_type::const_iterator, site_class_set_type::const_iterator> site_class_set_itpair_type;
  typedef std::pair<pair_class_set_type::const_iterator, pair_class_set_type::const_iterator> pair_class_set_itpair_type;

  ClusterHelper(const alps::Parameters& params);

  double cluster_coordinate(const unsigned int i, const unsigned int j) const { return cluster_coordinates_(i,j); }
  double cluster_momentum(const unsigned int i, const unsigned int j) const { return cluster_momenta_(i,j); }
  double lattice_momentum(const unsigned int i, const unsigned int j) const { return lattice_momenta_(i,j); }
  bool bipartite() const { return bipartite_; }
  unsigned int n_cluster_sites() const { return n_cluster_sites_; }
  unsigned int n_cluster_cells() const { return alps::graph_helper<alps::graph_helper<>::graph_type>::volume(); }
  unsigned int dimension() const { return alps::graph_helper<alps::graph_helper<>::graph_type>::dimension(); }
  site_class_set_itpair_type site_classes_real_space() const 
  { return site_class_set_itpair_type (site_class_set_real_space_.begin(), site_class_set_real_space_.end()); }
  pair_class_set_itpair_type pair_classes_real_space() const 
  { return pair_class_set_itpair_type(pair_class_set_real_space_.begin(), pair_class_set_real_space_.end()); }
  site_class_set_itpair_type site_classes_k_space() const 
  { return site_class_set_itpair_type(site_class_set_k_space_.begin(), site_class_set_k_space_.end()); }
  pair_class_set_itpair_type pair_classes_k_space() const 
  { return pair_class_set_itpair_type(pair_class_set_k_space_.begin(), pair_class_set_k_space_.end()); }


private:
  void setup_tables();
  void symmetry_analysis(const Matrix<double>::Type& basis, 
			 const Matrix<double>::Type& basis_inv, 
			 const Matrix<double>::Type& reciprocal_basis, 
			 const Matrix<double>::Type& reciprocal_basis_inv );
  std::pair<site_class_set_type, pair_class_set_type> setup_classes(const Matrix<site_t>::Type& symmetry_table) const;
  bool inside_cluster(vector_type v, const Matrix<double>::Type& b) const;
  bool inside_reduced_brillouin_zone(vector_type v, const Matrix<double>::Type& B) const;
  bool append_column_no_duplicates(Matrix<double>::Type& m, const vector_type& v) const;
  bool append_column_no_duplicates(Matrix<site_t>::Type& m, std::map<site_t, site_t>& smap) const;
  vector_type map_back_into_brillouin_zone(vector_type v) const;
  vector_type map_back_into_cluster(vector_type v, 
				    const Matrix<double>::Type& B_inv, const Matrix<double>::Type& B) const;
  unsigned int cluster_momentum_approximation(vector_type v, 
				      const Matrix<double>::Type& mom,
				      const Matrix<double>::Type& B) const;
  unsigned int identify_vertex(const vector_type& v, const Matrix<double>::Type& coordinate_matrix) const;
  unsigned int identify_sublattice(const Matrix<double>::Type& coordinate_matrix, const unsigned int i) const;
  double round(const double d) const;
  vector_type round(vector_type v) const;
  bool almost_equal(const double d1, const double d2) const { return fabs(d1-d2)<epsilon_; }
  vector_type matrix_vector_mult(const Matrix<double>::Type& M, const vector_type& v) const;
  template<typename T> typename Matrix<T>::Type matrix_inverse(typename Matrix<T>::Type M) const;
    
  Matrix<double>::Type cluster_coordinates_;
  Matrix<double>::Type cluster_momenta_;
  Matrix<double>::Type lattice_momenta_;
  site_class_set_type site_class_set_real_space_;
  pair_class_set_type pair_class_set_real_space_;
  site_class_set_type site_class_set_k_space_;
  pair_class_set_type pair_class_set_k_space_;
  unsigned int n_cluster_sites_;
  bool bipartite_;
  const double epsilon_;
};



class ClusterTransformer : protected ClusterHelper
{
public:
  ClusterTransformer(const alps::Parameters& params) :
    ClusterHelper(params) {}
    virtual ~ClusterTransformer() {} //pro forma virtual destructor.
  
  virtual matsubara_green_function_t operator()(const matsubara_green_function_t & G_omega, 
						const matsubara_green_function_t &G0_omega, 
						const double mu, const double h, 
						const double beta)=0; 
  virtual itime_green_function_t symmetrize(const itime_green_function_t& G, const bool paramagnet) const=0;
  virtual matsubara_green_function_t symmetrize(const matsubara_green_function_t& G, const bool paramagnet) const=0;
  virtual matsubara_green_function_t initial_G0(const alps::Parameters& parms); 
};



class DCATransformer : public ClusterTransformer
{
public:
  DCATransformer(const alps::Parameters& params) :
    ClusterTransformer(params) {}
  virtual ~DCATransformer(){} //pro forma virtual destructor
  
  template<typename T> green_function<T> symmetrize_real_space(const green_function<T>& G) const 
  { return symmetrize_helper(G, pair_classes_real_space()); }
  template<typename T> green_function<T> symmetrize_k_space(const green_function<T>& G) const 
  { return symmetrize_helper(G, pair_classes_k_space()); }
  
  virtual matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const=0; 
  virtual matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const=0;
  virtual itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const=0; 
  virtual itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const=0;
  virtual matsubara_green_function_t transform_fourpoint(const matsubara_green_function_t& G_omega_fourpoint,
							 const matsubara_green_function_t& G_omega, 
							 const matsubara_green_function_t& G0_omega, 
							 const double mu, const double h, 
							 const double beta) const =0; 
  std::complex<double>* exp_iksite(const int sign) const;

private:
  template<typename T>
  green_function<T> symmetrize_helper(const green_function<T>& G, 
			       const std::pair<std::set<pair_class_type>::const_iterator, 
			       std::set<pair_class_type>::const_iterator> pair_classes) const;
};



class SimpleDCATransformer : public DCATransformer
{
public:
  SimpleDCATransformer(const alps::Parameters& params) :
    DCATransformer(params), t_(boost::lexical_cast<double>(params["t"])), tprime_(params.defined("t'")?boost::lexical_cast<double>(params["t'"]):0) {}
  
  matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega, 
					const matsubara_green_function_t &G0_omega, 
					const double mu, const double h, const double beta);
  itime_green_function_t symmetrize(const itime_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t symmetrize(const matsubara_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const 
  { return transform_coordinates<std::complex<double>, Backward>(G); }
  matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const 
  { return transform_coordinates<std::complex<double>, Forward>(G); }
  itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const 
  { return transform_coordinates<double, Backward>(G); }
  itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const 
  { return transform_coordinates<double, Forward>(G); }
  
  matsubara_green_function_t initial_G0(const alps::Parameters& parms);  
  matsubara_green_function_t transform_fourpoint(const matsubara_green_function_t& G_omega_foupoint,
						 const matsubara_green_function_t &G_omega, 
						 const matsubara_green_function_t &G0_omega, 
						 const double mu, const double h, 
						 const double beta) const; 
  hifreq_moments k_space_moments(const double &mu) const;
  void compute_epsav(std::vector<double> &epsav, std::vector<double> &epsavsq) const;
 private:
  template<typename T, class D>
  green_function<T> transform_coordinates(const green_function<T>& G) const;
  template <class T>
  green_function<T> spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const;

  double epsilon(const vector_type k) const
  {  //cosine dispersion (nearest nbr hopping) for all dimensions. t' only valid/reasonable in two dimensions for nn hopping.
    double eps = 0;
    for (unsigned int i=0; i<dimension(); ++i) 
      eps -= 2*t_*cos(k[i]);
    /*if (tprime_!=0.) {
      switch(dimension()) {
      case 1 :
	eps -= 2*tprime_*cos(k[0]*2.);
	break;
      case 2 :
	eps -= 2*tprime_*cos(k[0]+k[1]);
	eps -= 2*tprime_*cos(k[0]-k[1]);
	break;
      case 3 :
	eps -= 2*tprime_*cos(k[0]+k[1]);
	eps -= 2*tprime_*cos(k[0]-k[1]);
	eps -= 2*tprime_*cos(k[1]+k[2]);
	eps -= 2*tprime_*cos(k[1]-k[2]);
	eps -= 2*tprime_*cos(k[0]+k[2]);
	eps -= 2*tprime_*cos(k[0]-k[2]);
	break;
      }
      }*/
    if(tprime_!=0){
      assert(dimension()==2);
      eps -= 4*tprime_*cos(k[0])*cos(k[1]);
    }
    return eps; 
  }
  
  const double t_, tprime_;
  struct Forward;
  struct Backward;
  
};



class AntiFerroDCATransformer : public DCATransformer
{
public:
  AntiFerroDCATransformer(const alps::Parameters& params) :
    DCATransformer(params), t_(boost::lexical_cast<double>(params["t"])) {}
  
  matsubara_green_function_t operator()(const matsubara_green_function_t &G_omega, 
					const matsubara_green_function_t &G0_omega, 
					const double mu, const double h, const double beta);
  itime_green_function_t symmetrize(const itime_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }
  matsubara_green_function_t symmetrize(const matsubara_green_function_t& G, const bool paramagnet) const
  { return spin_symmetrizer(G, paramagnet); }    
  matsubara_green_function_t transform_into_real_space(const matsubara_green_function_t& G) const 
  { return transform_coordinates<std::complex<double>, Backward>(G); }
  matsubara_green_function_t transform_into_k_space(const matsubara_green_function_t& G) const 
  { return transform_coordinates<std::complex<double>, Forward>(G); }
  itime_green_function_t transform_into_real_space(const itime_green_function_t& G) const 
  { return transform_coordinates<double, Backward>(G); }
  itime_green_function_t transform_into_k_space(const itime_green_function_t& G) const 
  { return transform_coordinates<double, Forward>(G); }
  
  matsubara_green_function_t initial_G0(const alps::Parameters& parms);  
  matsubara_green_function_t transform_fourpoint(const matsubara_green_function_t& G_omega_fourpoint,
						 const matsubara_green_function_t &G_omega, 
						 const matsubara_green_function_t &G0_omega, 
						 const double mu, const double h, 
						 const double beta) const
  {
    boost::throw_exception(std::logic_error("transform fourpoint function for antiferromagnetic phases not yet implemented"));
    std::cout<<G_omega_fourpoint<<std::endl;
    std::cout<<G_omega<<std::endl;
    std::cout<<G0_omega<<std::endl;
    std::cout<<mu<<std::endl;
    std::cout<<h<<std::endl;
    std::cout<<beta<<std::endl;
    return G_omega;
  } 

private:
  template <class T>
  green_function<T> spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const;
  template<typename T, class D>
  green_function<T> transform_coordinates(const green_function<T>& G) const;
  template<typename T>
  green_function<T> green_function_inverse(const green_function<T>& G) const;
  template<typename T>
  green_function<T> green_function_add(const green_function<T>& G1, const green_function<T>& G2) const;
  template<typename T>
  green_function<T> green_function_subtract(const green_function<T>& G1, const green_function<T>& G2) const;
  template<typename T>
  green_function<T> green_function_divide(const green_function<T>& G1, const T d) const;
  

  double epsilon(const vector_type k) const
  {
    double eps = 0;
    for (unsigned int i=0; i<dimension(); ++i) 
      eps -= 2*t_*cos (k[i]);
    return eps; 
  }
  
  const double t_;
  
  struct Forward;
  struct Backward;
  
};






//implementation (see also clusterhelper.C, clustertransformer.C)




inline double ClusterHelper::round(const double d) const
{
  int rounded = int(d > 0.0 ? d + 0.5 : d - 0.5);
  if (almost_equal(rounded, d))
    return rounded;
  else
    return d;
}



inline ClusterHelper::vector_type ClusterHelper::round(vector_type v) const
{
  for (unsigned int d=0; d<v.size(); ++d) 
    v[d] = round(v[d]);
  return v;
}



inline ClusterHelper::vector_type ClusterHelper::matrix_vector_mult(const Matrix<double>::Type& M, const vector_type& v) const
{
  vector_type vn(v.size());
  for (unsigned int i=0; i<v.size(); ++i) {
    vn[i] = 0.;
    for (unsigned int j=0; j<v.size(); ++j) {
      vn[i] += M(j,i)*v[j];
    }
  }
  return vn;
}


template<typename T>
inline typename ClusterHelper::Matrix<T>::Type ClusterHelper::matrix_inverse(typename Matrix<T>::Type M) const 
{
  Matrix<double>::Type M_inv(M.size1(), M.size2());
  for (unsigned int i=0; i<M.size1(); ++i) 
    for (unsigned int j=0; j<M.size2(); ++j) 
      M_inv(i,j) = i==j;
  boost::numeric::bindings::lapack::gesv(M, M_inv);
  return M_inv;
}




template<typename T>
green_function<T> DCATransformer::symmetrize_helper(const green_function<T>& G, 
					     const std::pair<std::set<pair_class_type>::const_iterator, 
					     std::set<pair_class_type>::const_iterator> pair_classes) const
{
  unsigned int nflavor=G.nflavor();
  unsigned int nsite=G.nsite();
  unsigned int nfreq=G.nfreq();
  green_function<T> G_new(nfreq, nsite, nflavor);
  std::set<pair_class_type>::const_iterator itc, itc_end;
  for(boost::tie(itc, itc_end)=pair_classes; itc!=itc_end; itc++) {
    for (spin_t flavor=0; flavor<nflavor; ++flavor) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
	T av=0;
	for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++)
	  av += G(freq, it->first, it->second, flavor);
	av /= itc->size();
	for(pair_class_type::const_iterator it=itc->begin(); it!=itc->end(); it++)
	  G_new(freq, it->first, it->second, flavor) = av;
      }
    }
  }
  return G_new;
}



struct SimpleDCATransformer::Forward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G, 
		  const itime_index_t t, const site_t i, const site_t j, const site_t k, const spin_t flavor,
		  double factor) {
    G_new(t, k, k, flavor) += factor * G(t, i, j, flavor);
  }
};


  
struct SimpleDCATransformer::Backward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G, 
		  const itime_index_t t, const site_t i, const site_t j, const site_t k, const spin_t flavor,
		  double factor) {
    G_new(t, i, j, flavor) += factor * G(t, k, k, flavor);
  }
};
  


template<typename T, class D>
green_function<T> SimpleDCATransformer::transform_coordinates(const green_function<T>& G) const 
{
  unsigned int nflavor=G.nflavor();
  unsigned int nsite=G.nsite();
  unsigned int nt=G.ntime();
  assert(nsite==n_cluster_sites());
  green_function<T> G_new(nt, nsite, nflavor);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
	for (itime_index_t t=0; t<nt; ++t) {
	  G_new(t, i, j, flavor) = 0.;
	}
      }
    }
  }
  for (site_t k=0; k<nsite; ++k) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
	double phase = 0.;
	for (unsigned int d=0; d<dimension(); ++d) 
	  phase += cluster_momentum(k,d)*(cluster_coordinate(i,d)-cluster_coordinate(j,d));
	for (itime_index_t t=0; t<nt; ++t) {
	  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
	    D::add(G_new, G, t, i, j, k, flavor, cos(phase)/double(nsite));
	  }
	}
      }
    }
  }
  return G_new;
}


template<class T>
green_function<T> SimpleDCATransformer::spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const
{
  green_function<T> G_new(G);
  if (paramagnet) {
    if(G_new.nflavor()>1 && G_new.nflavor()%2==0) {
      for(spin_t flavor=0; flavor<G_new.nflavor(); flavor+=2) {
	for (site_t k=0; k<n_cluster_sites(); ++k) {
	  for(frequency_t freq=0; freq<G_new.nfreq(); ++freq) {
	    G_new(freq, k, k, flavor)= 0.5*(G_new(freq, k, k, flavor)+G_new(freq, k, k, flavor+1));
	    G_new(freq, k, k, flavor+1)=G_new(freq, k, k, flavor);
	  }
	}
      }
    }
  }
  return G_new;
}


template<typename T, class D>
green_function<T> AntiFerroDCATransformer::transform_coordinates(const green_function<T>& G) const 
{
  unsigned int nflavor=G.nflavor();
  unsigned int nsite=G.nsite();
  unsigned int nt=G.ntime();
  assert(nsite==n_cluster_sites());
  green_function<T> G_new(nt, nsite, nflavor);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t i=0; i<nsite; ++i) {
      for (site_t j=0; j<nsite; ++j) {
	for (itime_index_t t=0; t<nt; ++t) {
	  G_new(t, i, j, flavor) = 0.;
	}
      }
    }
  }
  for (site_t k=0; k<nsite; k+=2) {
    for (site_t q1=0; q1<=1; ++q1) {
      for (site_t q2=0; q2<=1; ++q2) {
	for (site_t i=0; i<nsite; ++i) {
	  for (site_t j=0; j<nsite; ++j) {
	    double phase = 0.;
	    for (unsigned int d=0; d<dimension(); ++d) 
	      phase += cluster_momentum(k+q1,d)*cluster_coordinate(i,d)-cluster_momentum(k+q2,d)*cluster_coordinate(j,d);
	    for (itime_index_t t=0; t<nt; ++t) {
	      for (spin_t flavor=0; flavor<nflavor; ++flavor) {
		D::add(G_new, G, t, i, j, k+q1, k+q2, flavor, cos(phase)/double(nsite));
	      }
	    }
	  }
	}
      }
    }
  }
  return G_new;
}


struct AntiFerroDCATransformer::Forward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G, 
		  const itime_index_t t, const site_t i, const site_t j, const site_t k1, const site_t k2, const spin_t flavor,
		  double factor) {
    G_new(t, k1, k2, flavor) += factor * G(t, i, j, flavor);
  }
};


  
struct AntiFerroDCATransformer::Backward {
  template<typename T>
  static void add(green_function<T>& G_new, const green_function<T>& G, 
		  const itime_index_t t, const site_t i, const site_t j, const site_t k1, const site_t k2, const spin_t flavor,
		  double factor) {
    G_new(t, i, j, flavor) += factor * G(t, k1, k2, flavor);
  }
};
  


template<class T>
green_function<T> AntiFerroDCATransformer::spin_symmetrizer(const green_function<T>& G, const bool paramagnet) const
{
  if(paramagnet==1){ std::cerr<<"Sebastian: paramagnet is true, should we really be using the AFM transformer?"<<std::endl; }
  green_function<T> G_new(G);
  if(G_new.nflavor()>1 && G_new.nflavor()%2==0) {
    for(spin_t flavor=0; flavor<G_new.nflavor(); flavor+=2) {
      for (site_t k=0; k<n_cluster_sites(); k+=2) {
	for (site_t q1=0; q1<=1; ++q1) {
	  for (site_t q2=0; q2<=1; ++q2) {
	    double sign = (q1==q2) ? 1. : -1.;
	    for(frequency_t freq=0; freq<G_new.nfreq(); ++freq) {
	      G_new(freq, k+q1, k+q2, flavor) = 0.5 * (G_new(freq, k+q1, k+q2, flavor) + sign * G_new(freq, k+q1, k+q2, flavor+1));
	      G_new(freq, k+q1, k+q2, flavor+1) = sign * G_new(freq, k+q1, k+q2, flavor);
	    }
	  }
	}
      }
    }
  }
  return G_new;
}


template<typename T>
green_function<T> AntiFerroDCATransformer::green_function_inverse(const green_function<T>& G) const
{
  green_function<T> G_inv(G);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (int k=0; k<G.nsite(); k+=2) {
	T detG = G(t, k, k, f)*G(t, k+1, k+1, f) - G(t, k+1, k, f)*G(t, k, k+1, f);
	G_inv(t, k, k, f) = G(t, k+1, k+1, f)/detG;
	G_inv(t, k+1, k+1, f) = G(t, k, k, f)/detG;
	G_inv(t, k+1, k, f) = -G(t, k+1, k, f)/detG;
	G_inv(t, k, k+1, f) = -G(t, k, k+1, f)/detG;
      }
    }
  }
  return G_inv;
}


template <typename T>
green_function<T> AntiFerroDCATransformer::green_function_add(const green_function<T>& G1,
							      const green_function<T>& G2) const
{
  green_function<T> G(G1);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (int k=0; k<G.nsite(); k+=2) {
	G(t, k, k, f) += G2(t, k, k, f);
	G(t, k+1, k+1, f) += G2(t, k+1, k+1, f);
	G(t, k+1, k, f) += G2(t, k+1, k, f);
	G(t, k, k+1, f) += G2(t, k, k+1, f);
      }
    }
  }
  return G;
}



template <typename T>
green_function<T> AntiFerroDCATransformer::green_function_subtract(const green_function<T>& G1,
								   const green_function<T>& G2) const
{
  green_function<T> G(G1);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (int k=0; k<G.nsite(); k+=2) {
	G(t, k, k, f) -= G2(t, k, k, f);
	G(t, k+1, k+1, f) -= G2(t, k+1, k+1, f);
	G(t, k+1, k, f) -= G2(t, k+1, k, f);
	G(t, k, k+1, f) -= G2(t, k, k+1, f);
      }
    }
  }
  return G;
}



template <typename T>
green_function<T> AntiFerroDCATransformer::green_function_divide(const green_function<T>& G1,
								 const T d) const
{
  green_function<T> G(G1);
  for (unsigned int t=0; t<G.ntime(); ++t) {
    for (spin_t f=0; f<G.nflavor(); ++f) {
      for (int k=0; k<G.nsite(); k+=2) {
	G(t, k, k, f) /= d;
	G(t, k+1, k+1, f) /= d;
	G(t, k+1, k, f) /= d;
	G(t, k, k+1, f) /= d;
      }
    }
  }
  return G;
}



#endif /*ALPS_CLUSTER_H*/
