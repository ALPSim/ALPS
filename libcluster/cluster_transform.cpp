/*****************************************************************************
 *
 * Copyright (C) 2007-2015 by Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#include <iostream>
#include <string>

#include "alps/numeric/vector_functions.hpp"
#include "alps/hdf5/vector.hpp"
#include "alps/params.hpp"

#include "cluster_transform.h"

cluster_matsubara_kspace_gf NormalStateClusterTransformer::initial_G0(const alps::params& parms) {
  std::cout<<"reading initial G0_omega from: "<<parms["G0_OMEGA_INPUT"].as<std::string>()<<std::endl;
  int nfreq=parms["NMATSUBARA"];
  int n_orbital=parms["FLAVORS"];
  int n_site=parms["SITES"];
  double beta=parms["BETA"];

  cluster_matsubara_kspace_gf G0_omega(alps::gf::omega_k_sigma_gf(
      alps::gf::matsubara_positive_mesh(beta,nfreq),
      momentum_mesh(),
      alps::gf::index_mesh(n_orbital))
  );

  alps::hdf5::archive ar(parms["G0_OMEGA_INPUT"].as<std::string>());
  G0_omega.load(ar, "/G0");
  return G0_omega;
}
alps::gf::momentum_index_mesh NormalStateClusterTransformer::momentum_mesh() const{
  alps::gf::momentum_index_mesh::container_type points(boost::extents[n_cluster_sites()][dimension()]);
  for(int i=0;i<n_cluster_sites();++i){
    for(int j=0;j<dimension();++j){
      points[i][j]=cluster_momentum(i,j);
    }
  }
  return alps::gf::momentum_index_mesh(points);
}
void NormalStateClusterTransformer::spin_symmetrize(cluster_itime_kspace_gf& gf) const{
  if(gf.mesh3().extent()!=2) throw std::runtime_error("mesh should only have two flavors.");
  for (alps::gf::itime_index om=alps::gf::itime_index(0); om<gf.mesh1().extent(); ++om) {
    for (alps::gf::momentum_index ii=alps::gf::momentum_index(0); ii<gf.mesh2().extent(); ++ii) {
      double mean=(gf(om,ii,alps::gf::index(0))+gf(om,ii,alps::gf::index(1)))/2.;
      gf(om,ii,alps::gf::index(0))=mean;
      gf(om,ii,alps::gf::index(1))=mean;
    }
  }
}
///perform spin symmetrization of matsubara k-space Green's function
void NormalStateClusterTransformer::spin_symmetrize(cluster_matsubara_kspace_gf& gf) const{
  if(gf.mesh3().extent()!=2) throw std::runtime_error("mesh should only have two flavors.");
  for (alps::gf::matsubara_index om=alps::gf::matsubara_index(0); om<gf.mesh1().extent(); ++om) {
    for (alps::gf::momentum_index ii=alps::gf::momentum_index(0); ii<gf.mesh2().extent(); ++ii) {
      std::complex<double> mean=(gf(om,ii,alps::gf::index(0))+gf(om,ii,alps::gf::index(1)))/2.;
      gf(om,ii,alps::gf::index(0))=mean;
      gf(om,ii,alps::gf::index(1))=mean;
    }
  }
}
///perform spin symmetrization of imaginary time k-space cluster Green's function
void NormalStateClusterTransformer::spin_symmetrize(cluster_itime_realspace_gf& gf) const{
  if(gf.mesh3().extent()!=2) throw std::runtime_error("mesh should only have two flavors.");
  for (alps::gf::itime_index om=alps::gf::itime_index(0); om<gf.mesh1().extent(); ++om) {
    for (alps::gf::real_space_index ii=alps::gf::real_space_index(0); ii<gf.mesh2().extent(); ++ii) {
      for (alps::gf::real_space_index jj=alps::gf::real_space_index(0); jj<gf.mesh3().extent(); ++jj) {
        double mean=(gf(om,ii,jj,alps::gf::index(0))+gf(om,ii,jj,alps::gf::index(1)))/2.;
        gf(om,ii,jj,alps::gf::index(0))=mean;
        gf(om,ii,jj,alps::gf::index(1))=mean;
      }
    }
  }
}
///perform spin symmetrization of matsubara k-space Green's function
void NormalStateClusterTransformer::spin_symmetrize(cluster_matsubara_realspace_gf& gf) const{
  if(gf.mesh3().extent()!=2) throw std::runtime_error("mesh should only have two flavors.");
  for (alps::gf::matsubara_index om=alps::gf::matsubara_index(0); om<gf.mesh1().extent(); ++om) {
    for (alps::gf::real_space_index ii=alps::gf::real_space_index(0); ii<gf.mesh2().extent(); ++ii) {
      for (alps::gf::real_space_index jj=alps::gf::real_space_index(0); jj<gf.mesh3().extent(); ++jj) {
        std::complex<double> mean=(gf(om,ii,jj,alps::gf::index(0))+gf(om,ii,jj,alps::gf::index(1)))/2.;
        gf(om,ii,jj,alps::gf::index(0))=mean;
        gf(om,ii,jj,alps::gf::index(1))=mean;
      }
    }
  }
}
struct Forward {
  template<class GF_RS, class GF_KS> static void add(GF_RS & G_new, const GF_KS& G,
      const int t, const int i, const int j, const int k, const int flavor,
      double factor) {
    G_new(GF_KS::mesh1_type::index_type(t), alps::gf::momentum_index(k), alps::gf::index(flavor)) +=
        factor * G(GF_RS::mesh1_type::index_type(t), alps::gf::real_space_index(i), alps::gf::real_space_index(j), alps::gf::index(flavor));
  }
  template<class GF_RS, class GF_KS>
  static void add(GF_RS& G_new, const GF_KS& G,
      const int t, const int i, const int j, const int k,
      double factor) {
    for(alps::gf::index s1(0); s1<2;++s1){
      for(alps::gf::index s2(0); s2<2;++s2){
        G_new(GF_KS::mesh1_type::index_type(t), alps::gf::momentum_index(k), s1,s2) +=
            factor * G(GF_RS::mesh1_type::index_type(t), alps::gf::real_space_index(i), alps::gf::real_space_index(j), s1, s2);
      }
    }
  }
};



struct Backward {
  template<class GF_RS, class GF_KS> static void add(GF_RS & G_new, const GF_KS& G,
      const int t, const int i, const int j, const int k, const int flavor,
      double factor) {
    G_new(GF_RS::mesh1_type::index_type(t), alps::gf::real_space_index(i), alps::gf::real_space_index(j), alps::gf::index(flavor)) +=
        factor * G(GF_RS::mesh1_type::index_type(t), alps::gf::momentum_index(k), alps::gf::index(flavor));
  }
  template<class GF_RS, class GF_KS>
  static void add(GF_RS& G_new, const GF_KS& G,
      const int t, const int i, const int j, const int k,
      double factor) {
    for(alps::gf::index s1(0); s1<2;++s1){
      for(alps::gf::index s2(0); s2<2;++s2){
        G_new(GF_RS::mesh1_type::index_type(t), alps::gf::real_space_index(i), alps::gf::real_space_index(j), s1,s2) +=
            factor * G(GF_RS::mesh1_type::index_type(t), alps::gf::momentum_index(k), s1, s2);
      }
    }
  }
};

cluster_matsubara_realspace_gf NormalStateClusterTransformer::transform_into_real_space(const cluster_matsubara_kspace_gf& G) const{
  throw std::logic_error("not implemented");
}
cluster_matsubara_kspace_gf NormalStateClusterTransformer::transform_into_k_space(const cluster_matsubara_realspace_gf& G) const{ throw std::logic_error("not implemented");};
cluster_itime_realspace_gf NormalStateClusterTransformer::transform_into_real_space(const cluster_itime_kspace_gf& G) const{ throw std::logic_error("not implemented");};
cluster_itime_kspace_gf NormalStateClusterTransformer::transform_into_k_space(const cluster_itime_realspace_gf& G) const{ throw std::logic_error("not implemented");};

/*cluster_nambu_matsubara_kspace_gf NambuClusterTransformer::initial_G0_nambu(const alps::params& parms) {
  std::cout<<"reading initial (NAMBU) G0_omega from: "<<parms["G0_OMEGA_INPUT"]<<std::endl;
  int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  int n_orbital=boost::lexical_cast<int>(parms["FLAVORS"]);
  int n_site=boost::lexical_cast<int>(parms["SITES"]);
  double beta =boost::lexical_cast<double>(parms["BETA"]);
  if(n_orbital!=2) throw std::invalid_argument("dunno how to handle nambu in >2 orbitals!");
  nambu_matsubara_green_function_t G0_omega(n_matsubara, n_site, n_orbital, diagonal);
  std::string g0path=parms["G0OMEGA_INPUT"];
  G0_omega.read(g0path.c_str());
  return G0_omega;
}*/


//Don't forget to delete the transformer afterwards!
/*ClusterTransformer* ClusterTransformer::generate_transformer(const alps::params& parms)
{
  ClusterTransformer* trans_ptr;
  bool antiferromagnetic = parms["ANTIFERROMAGNET"]|false;
  std::string cluster_loop=parms["CLUSTER_LOOP"];
  if(cluster_loop=="DCA"){
    if (antiferromagnetic) {
      trans_ptr=new AntiFerroDCATransformer(parms);
    }
    else {
      trans_ptr=new SimpleDCATransformer(parms);
    }
  }else if(cluster_loop=="CDMFT"){
    std::cout<<"using CDMFT transformer"<<std::endl;
    if(static_cast<int>(parms["SITES"])!=4) throw std::runtime_error("only four-site CDMFT is implemented so far.");
    trans_ptr=new LadderCDMFTTransformer(parms);
  }else
    throw std::logic_error("transformer pointer could not be initalized: unknown cluster loop.");
  return trans_ptr;
}
 */



/*std::complex<double>* ClusterTransformer::exp_iksite(const int sign) const
{
  std::complex<double> *m = new std::complex<double>[n_cluster_sites()*n_cluster_sites()];
  for (int k=0; k<n_cluster_sites(); ++k) {
    for (int i=0; i<n_cluster_sites(); ++i) {
      double phase = 0.;
      for (std::size_t d=0; d<dimension(); ++d){
        phase += cluster_momentum(k,d)*cluster_coordinate(i,d);
      }
      phase *= sign;
      m[k+n_cluster_sites()*i]=std::complex<double>(cos(phase), sin(phase));
    }
  }
  return m;
}


void ClusterTransformer::write_hdf5(alps::hdf5::archive &ar, const std::string &p) const{
  ClusterHelper::write_hdf5(ar,p);
}*/




