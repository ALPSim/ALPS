/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2007-2008 by Sebastian Fuchs <fuchs@comp-phys.org>
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


#include "cluster.h"


matsubara_green_function_t ClusterTransformer::initial_G0(const alps::Parameters& parms) {
  std::cerr<<"reading initial G0_omega"<<std::endl;
  int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  int n_orbital=boost::lexical_cast<int>(parms["FLAVORS"]);
  int n_site=boost::lexical_cast<int>(parms["SITES"]);
  matsubara_green_function_t G0_omega(n_matsubara, n_site, n_orbital);
  G0_omega.read(parms["G0OMEGA_INPUT"].c_str());
  return G0_omega;
}



std::complex<double>* DCATransformer::exp_iksite(const int sign) const
{
  std::complex<double> *m = new std::complex<double>[n_cluster_sites()*n_cluster_sites()];
  for (site_t k=0; k<n_cluster_sites(); ++k) {
    for (site_t i=0; i<n_cluster_sites(); ++i) {
      double phase = 0.;
      for (int d=0; d<dimension(); ++d) 
	phase += cluster_momentum(k,d)*cluster_coordinate(i,d);
      phase *= sign;
      m[k+n_cluster_sites()*i]=std::complex<double>(cos(phase), sin(phase));
    }
  }
  return m;
}



matsubara_green_function_t SimpleDCATransformer::operator()(const matsubara_green_function_t &G_omega, 
							    const matsubara_green_function_t &G0_omega, 
							    const double mu, const double h, const double beta)
{
  if(h!=0){std::cerr<<"don't know how to handle magn. fields: h="<<h<<" exiting"<<std::endl; abort(); }
  unsigned int nflavor=G0_omega.nflavor();
  unsigned int nfreq=G0_omega.nfreq();
  assert(G0_omega.nsite()==n_cluster_sites());
  matsubara_green_function_t selfenergy(G0_omega);
  matsubara_green_function_t G0_omega_new(G0_omega);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
	selfenergy(freq, kc, kc, flavor)  = 1./G0_omega(freq, kc, kc, flavor);
	selfenergy(freq, kc, kc, flavor) -= 1./G_omega(freq, kc, kc, flavor);
	std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
	G0_omega_new(freq, kc, kc, flavor) = 0.;
	for (site_t kl=0; kl<n_cluster_cells(); ++kl) { 
	  vector_type k(dimension());
	  for (unsigned int d=0; d<dimension(); ++d) 
	    k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
	  G0_omega_new(freq, kc, kc, flavor) 
	    += 1./(iomega_n+mu-epsilon(k)-selfenergy(freq, kc, kc, flavor));
	}
	G0_omega_new(freq, kc, kc, flavor) /= n_cluster_cells();
	G0_omega_new(freq, kc, kc, flavor) = 
	  1. / ( 1./G0_omega_new(freq, kc, kc, flavor) + selfenergy(freq, kc, kc, flavor) );
      }
    }
  }
  return G0_omega_new;
}



matsubara_green_function_t SimpleDCATransformer::initial_G0(const alps::Parameters& parms)
{
  if (parms.defined("G0OMEGA_INPUT")){
    std::ifstream in_file(parms["G0OMEGA_INPUT"].c_str());
    if(in_file.good()){
      return ClusterTransformer::initial_G0(parms);
    }
    else
      std::cerr<<"WARNING: Input file could not be found, proceeding with free solution."<<std::endl;
  }
  { //proceed with computation of new G0 if we can't find input G0
    std::cerr<<"calculating initial G0_omega"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    int n_tau=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    double beta = static_cast<double>(parms["BETA"]);
    double mu = static_cast<double>(parms["MU"]);
    double h = static_cast<double>(parms["H"]);
    //bool paramagnet = static_cast<bool>(parms["PARAMAGNET"]);
    matsubara_green_function_t G_omega(n_matsubara,n_cluster_sites(),n_flavor);
    matsubara_green_function_t G0_omega(n_matsubara,n_cluster_sites(),n_flavor);  
    G0_omega.clear();
    for (spin_t flavor=0; flavor<n_flavor; ++flavor) {
      for (site_t i=0; i<n_cluster_sites(); ++i) {
	for (frequency_t freq=0; freq<n_matsubara; ++freq) {
	  G0_omega(freq, i, i, flavor) = 1.;
	}
      }
    }
    G_omega = G0_omega; 
    //evoke transformation with G_omega==G0_omega (selfenergy=0)
    matsubara_green_function_t G0_omega_new = operator()(G_omega, G0_omega, mu, h, beta);
    
    itime_green_function_t G0_tau(n_tau+1, n_cluster_sites(), n_flavor);
    boost::shared_ptr<FourierTransformer> fourier_ptr;
    G0_omega_new.write("G0_omega_input");
    FourierTransformer::generate_transformer(parms, fourier_ptr, k_space_moments(parms["MU"]));
    fourier_ptr->backward_ft(G0_tau, G0_omega_new);
    G0_tau.write("G0_tau_input");
    matsubara_green_function_t G0_omega_real = transform_into_real_space(G0_omega_new);
    G0_omega_real.write("G0_omega_input_realspace");
    itime_green_function_t G0_tau_real = transform_into_real_space(G0_tau);
    G0_tau_real.write("G0_tau_input_realspace");
    return G0_omega_new;
  }
}



matsubara_green_function_t SimpleDCATransformer::transform_fourpoint(const matsubara_green_function_t& G_omega_fourpoint,
								     const matsubara_green_function_t& G_omega,
								     const matsubara_green_function_t& G0_omega, 
								     const double mu, const double h, 
								     const double beta) const
{
  std::cerr << "transforming cluster green function\n";
  assert(h==0); 
  unsigned int nz = G0_omega.nflavor();
  unsigned int n_matsubara = G0_omega.nfreq(); 
  assert(G0_omega.nsite()==n_cluster_sites());
  matsubara_green_function_t G0_omega_fourpoint_inv(n_matsubara*n_matsubara, n_cluster_sites(), nz);
  G0_omega_fourpoint_inv.clear();
  //matsubara_green_function_t G0_omega_fourpoint(n_matsubara*n_matsubara, n_cluster_sites(), nz);
  //G0_omega_fourpoint.clear();
  for(spin_t z=0; z<nz; ++z) {
    int zbar = z%nz==0 ? z+1 : z-1;
    for(site_t k=0;k<n_cluster_sites();++k) {
      int kpq = k%2==0 ? k+1 : k-1;
      for (frequency_t om=0; om<n_matsubara; ++om) {
	G0_omega_fourpoint_inv(om+n_matsubara*om,k,k,z) = 1. / (G_omega(om,kpq,kpq,z)*G_omega(om,k,k,zbar));
	//G0_omega_fourpoint(om+n_matsubara*om,k,k,z) = G_omega(om,kpq,kpq,z)*G_omega(om,k,k,zbar);
      } 
    }
  }
  /*std::ofstream greenfourpstr("G0FourPoint");
  for (int om1 =0; om1<n_matsubara; ++om1) {
    for (int om2 =0; om2<n_matsubara; ++om2) {
      greenfourpstr << om1 << " " << om2 << " ";
      for (int k1=0; k1<n_cluster_sites(); ++k1) {
	for (int k2=0; k2<n_cluster_sites(); ++k2) {
	  for (int z=0; z<nz; ++z) {
	    greenfourpstr << G0_omega_fourpoint(om2+n_matsubara*om1, k1, k2, z).real() << " " 
			  << G0_omega_fourpoint(om2+n_matsubara*om1, k1, k2, z).imag() << " ";
	  }
	}
      }
      greenfourpstr << std::endl;
    }
    }*/
  matsubara_green_function_t G_omega_fourpoint_inv = fourpoint_green_function_inverse(G_omega_fourpoint, n_matsubara);
  matsubara_green_function_t two_particle_vertex(n_matsubara*n_matsubara, n_cluster_sites(), nz);
  two_particle_vertex = fourpoint_green_function_subtract(G0_omega_fourpoint_inv, G_omega_fourpoint_inv);
  matsubara_green_function_t selfenergy(n_matsubara, n_cluster_sites(), nz);
  selfenergy.clear();
  for (spin_t z=0; z<nz; ++z) {
    for (site_t k=0; k<n_cluster_sites(); ++k) {
      for (frequency_t om=0; om<n_matsubara; ++om) {
	selfenergy(om, k, k, z) = 1./G0_omega(om, k, k, z);
	selfenergy(om, k, k, z) -= 1./G_omega(om, k, k, z);
      }
    }
  }
  matsubara_green_function_t G_omega_fourpoint_new(n_matsubara*n_matsubara, n_cluster_sites(), nz);
  G_omega_fourpoint_new.clear();
  //matsubara_green_function_t Gtest(n_matsubara*n_matsubara, n_cluster_sites(), nz);
  //Gtest.clear();
  for (spin_t z=0; z<nz; ++z) {
    spin_t zbar = z%nz==0 ? z+1 : z-1;
    for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
      int kcpq = kc%2==0 ? kc+1 : kc-1;
      for (frequency_t om=0; om<n_matsubara; ++om) {
	std::complex<double> iomega_n(0, (2.*om+1.)*M_PI/beta);	
	for (site_t kl=0; kl<n_cluster_cells(); ++kl) { 
	  vector_type k(dimension());
	  vector_type k_plus_q(dimension());
  	  for (unsigned int d=0; d<dimension(); ++d) { 
	    k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
	    k_plus_q[d] = lattice_momentum(kl,d) + cluster_momentum(kcpq,d);
	  }
	  G_omega_fourpoint_new(om+n_matsubara*om, kc, kc, z) 
	    += (1./(iomega_n+mu-epsilon(k_plus_q)-selfenergy(om, kcpq, kcpq, z)))
	    * (1./(iomega_n+mu-epsilon(k)-selfenergy(om, kc, kc, zbar)));
	}
	G_omega_fourpoint_new(om+n_matsubara*om, kc, kc, z) /= n_cluster_cells();
	//Gtest(om+n_matsubara*om, kc, kc, z) = G_omega_fourpoint_new(om+n_matsubara*om, kc, kc, z);
	G_omega_fourpoint_new(om+n_matsubara*om, kc, kc, z) = 1. / G_omega_fourpoint_new(om+n_matsubara*om, kc, kc, z);
      }
    }
  }
  /*std::ofstream greenteststr("GtestFourPoint");
  for (int om1 =0; om1<n_matsubara; ++om1) {
    for (int om2 =0; om2<n_matsubara; ++om2) {
      greenteststr << om1 << " " << om2 << " ";
      for (int k1=0; k1<n_cluster_sites(); ++k1) {
	for (int k2=0; k2<n_cluster_sites(); ++k2) {
	  for (int z=0; z<nz; ++z) {
	    greenteststr << Gtest(om2+n_matsubara*om1, k1, k2, z).real() << " " 
			  << Gtest(om2+n_matsubara*om1, k1, k2, z).imag() << " ";
	  }
	}
      }
      greenteststr << std::endl;
    }
    }*/
  G_omega_fourpoint_new = fourpoint_green_function_subtract(G_omega_fourpoint_new, two_particle_vertex);
  G_omega_fourpoint_new = fourpoint_green_function_inverse(G_omega_fourpoint_new, n_matsubara);
  std::cerr << "lattice green function calculated\n";
  //exit(1);
  return G_omega_fourpoint_new;
}



matsubara_green_function_t AntiFerroDCATransformer::operator()(const matsubara_green_function_t &G_omega, 
							       const matsubara_green_function_t &G0_omega, 
							       const double mu, const double h, const double beta)
{
  unsigned int nflavor=G0_omega.nflavor();
  unsigned int nfreq=G0_omega.nfreq();
  assert(G0_omega.nsite()==n_cluster_sites());
  assert(nflavor>1 && nflavor%2 == 0);
  matsubara_green_function_t G0_omega_new(nfreq,n_cluster_sites(),nflavor);
  G0_omega_new.clear();
  matsubara_green_function_t G0_omega_inv = green_function_inverse(G0_omega);
  matsubara_green_function_t G_omega_inv = green_function_inverse(G_omega);
  matsubara_green_function_t selfenergy = green_function_subtract(G0_omega_inv, G_omega_inv);
  for (site_t kl=0; kl<n_cluster_cells(); ++kl) { 
    matsubara_green_function_t G_kl(G0_omega);
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      for (site_t kc=0; kc<n_cluster_sites(); kc+=2) {
	vector_type k0(dimension());
	for (unsigned int d=0; d<dimension(); ++d) 
	  k0[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
	vector_type k1(dimension());
	for (unsigned int d=0; d<dimension(); ++d) 
	  k1[d] = lattice_momentum(kl,d) + cluster_momentum(kc+1,d);
	for (spin_t flavor=0; flavor<nflavor; flavor+=2) { 
	  G_kl(freq,kc,kc,flavor) = iomega_n+mu-epsilon(k0);
	  G_kl(freq,kc,kc,flavor+1) = iomega_n+mu-epsilon(k0);
	  G_kl(freq,kc+1,kc+1,flavor) = iomega_n+mu-epsilon(k1);
	  G_kl(freq,kc+1,kc+1,flavor+1) = iomega_n+mu-epsilon(k1);
	  G_kl(freq,kc+1,kc,flavor) = h/2;
	  G_kl(freq,kc+1,kc,flavor+1) = -h/2;
	  G_kl(freq,kc,kc+1,flavor) = h/2;
	  G_kl(freq,kc,kc+1,flavor+1) = -h/2;
	}
      }
    }
    G_kl = green_function_subtract(G_kl, selfenergy);
    G_kl = green_function_inverse(G_kl);
    G0_omega_new = green_function_add(G0_omega_new, G_kl);
  }
  G0_omega_new = green_function_divide(G0_omega_new, std::complex<double>(n_cluster_cells()));
  G0_omega_new = green_function_inverse(G0_omega_new);
  G0_omega_new = green_function_add(G0_omega_new, selfenergy);
  G0_omega_new = green_function_inverse(G0_omega_new);
  return G0_omega_new;
}



matsubara_green_function_t AntiFerroDCATransformer::initial_G0(const alps::Parameters& parms)
{
  if (parms.defined("G0OMEGA_INPUT")) 
    return ClusterTransformer::initial_G0(parms);
  else {
    std::cerr<<"calculating initial G0_omega"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    int n_tau=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    double beta = static_cast<double>(parms["BETA"]);
    double mu = static_cast<double>(parms["MU"]);
    double h = static_cast<double>(parms["H"]);
    double h_init = parms.value_or_default("H_INIT", 0.);
    matsubara_green_function_t G_omega(n_matsubara,n_cluster_sites(),n_flavor);
    matsubara_green_function_t G0_omega(n_matsubara,n_cluster_sites(),n_flavor);  
    G0_omega.clear();
    for (spin_t flavor=0; flavor<n_flavor; ++flavor) {
      for (site_t i=0; i<n_cluster_sites(); ++i) {
	for (frequency_t freq=0; freq<n_matsubara; ++freq) {
	  G0_omega(freq, i, i, flavor) = 1.;
	}
      }
    }
    G_omega = G0_omega; 
    //evoke transformation with G_omega==G0_omega (selfenergy=0)
    matsubara_green_function_t G0_omega_new = operator()(G_omega, G0_omega, mu, h_init+h, beta);
    boost::shared_ptr<FourierTransformer> fourier_ptr;
    FourierTransformer::generate_transformer(parms, fourier_ptr);
    G0_omega_new.write("G0_omega_input");
    transform_into_real_space(G0_omega_new).write("G0_omega_input_realspace");
    itime_green_function_t G0_tau(n_tau, n_cluster_sites(), n_flavor);
    fourier_ptr->backward_ft(G0_tau, G0_omega_new);
    G0_tau.write("G0_tau_input");
    transform_into_real_space(G0_tau).write("G0_tau_input_realspace");
    return G0_omega_new;
  }
}
  
hifreq_moments SimpleDCATransformer::k_space_moments(const double &mu) const{
  std::cout<<"computing moments for k-space fourier transform"<<std::endl;
  hifreq_moments moments;
  moments.c0.resize(2);
  moments.c1.resize(2);
  moments.c2.resize(2);
  moments.c3.resize(2);
  for(int s=0;s<2;++s){
    moments.c0[s].resize(n_cluster_sites());
    moments.c1[s].resize(n_cluster_sites());
    moments.c2[s].resize(n_cluster_sites());
    moments.c3[s].resize(n_cluster_sites());
    for(int i=0;i<n_cluster_sites();++i){
      moments.c0[s][i].resize(n_cluster_sites(),0.);
      moments.c1[s][i].resize(n_cluster_sites(),0.);
      moments.c2[s][i].resize(n_cluster_sites(),0.);
      moments.c3[s][i].resize(n_cluster_sites(),0.);
    }
  }

  std::vector<double> epsav  (n_cluster_sites()); //average of the band structure on this patch
  std::vector<double> epssqav(n_cluster_sites());
  compute_epsav(epsav, epssqav);

  for(int s=0;s<2;++s){ 
  for(int i=0;i<n_cluster_sites();++i){ //simple DCA has just the diagonal occupied
    moments.c1[s][i][i]=1.; //this is the commutator term
    moments.c2[s][i][i]=-(mu-epsav[i]); //this is the second order term
    moments.c3[s][i][i]=mu*mu-2*mu*epsav[i]+epssqav[i]; //this is the third order term
  }
  }
  return moments;
}
void SimpleDCATransformer::compute_epsav(std::vector<double> &epsav, std::vector<double> &epssqav) const{
  epsav  .resize(n_cluster_sites());
  epssqav.resize(n_cluster_sites());
  for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
    for (site_t kl=0; kl<n_cluster_cells(); ++kl) {
      vector_type k(dimension());
      for (unsigned int d=0; d<dimension(); ++d)
        k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
      epsav  [kc]+= epsilon(k);
      epssqav[kc]+= epsilon(k)*epsilon(k);
    }
    epsav  [kc] /= n_cluster_cells();
    epssqav[kc] /= n_cluster_cells();
  }
}
