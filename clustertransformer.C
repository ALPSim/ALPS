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


bool adjust_chemical_potential(const double &density, double &mu, bool &damping_mu, 
			       bool &computing_mu, int &mu_it, double &n_old, double &n_target, 
			       double &mu_old, double &mu_new, double &mu_orig);



matsubara_green_function_t ClusterTransformer::initial_G0(const alps::Parameters& parms) {
  std::cout<<"reading initial G0_omega"<<std::endl;
  int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  int n_orbital=boost::lexical_cast<int>(parms["FLAVORS"]);
  int n_site=boost::lexical_cast<int>(parms["SITES"]);
  matsubara_green_function_t G0_omega(n_matsubara, n_site, n_orbital);
  G0_omega.read(parms["G0OMEGA_INPUT"].c_str());
  return G0_omega;
}



void DCATransformer::generate_transformer(const alps::Parameters& parms, boost::shared_ptr<DCATransformer> &trans_ptr)
{
  bool antiferromagnetic = parms.value_or_default("ANTIFERROMAGNET", false);
  if (antiferromagnetic) {
    std::cout<<"using antiferromagnetic DCA transformer"<<std::endl;
    trans_ptr.reset(new AntiFerroDCATransformer(parms));
  }
  else {
    std::cout<<"using simple DCA transformer"<<std::endl;
    trans_ptr.reset(new SimpleDCATransformer(parms));
  }
}



DCATransformer::DCATransformer(const alps::Parameters& params) 
  : ClusterTransformer(params), 
    t_(params.defined("W") ? 1. : static_cast<double>(params["t"])), 
    t2_(params.value_or_default("t'", params.value_or_default("t2", 0.))),
    t3_(params.value_or_default("t3", 0.)),
    t4_(params.value_or_default("t4", 0.)),
    bandwidth_(params.value_or_default("W", 0.)),
    bandwidth_scale_factor_(1.),
    half_filling_(0.),
    epsav_(n_cluster_sites()), 
    epssqav_(n_cluster_sites()),
    dispersion_moments_(2, 0.)
{
  if (params.defined("t") && params.defined("W"))
    boost::throw_exception(std::invalid_argument("both t and W are defined"));
  //calculating dispersion moments, bandwidth and half filling
  std::vector<double> energies;
  for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
    for (site_t kl=0; kl<n_lattice_momenta(); ++kl) {
      vector_type k(dimension());
      for (unsigned int d=0; d<dimension(); ++d)
	k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
      double weight = lattice_momentum(kl, dimension());
      double E = epsilon(k);
      epsav_  [kc]+= weight*E;
      epssqav_[kc]+= weight*E*E;
      energies.push_back(E); 
    }    
    epsav_  [kc] /= n_cluster_cells();
    epssqav_[kc] /= n_cluster_cells();
    dispersion_moments_[0] += epsav_[kc];
    dispersion_moments_[1] += epssqav_[kc];
  }
  dispersion_moments_[0] /= n_cluster_sites();
  dispersion_moments_[1] /= n_cluster_sites();
  double hf = 0.;
  double bw = 0.;
  if (t2_==0. && t3_==0. && t4_==0.) {
    bw = dimension()*4*t_;
  }
  else {
    std::sort(energies.begin(), energies.end());
    double min = energies.front();
    double max = energies.back();
    bw = max - min;
    if (energies.size()%2) 
      hf = energies[energies.size()/2];
    else 
      hf = (energies[energies.size()/2-1] + energies[energies.size()/2])/2. ;
  }
  if (params.value_or_default("SHIFT_TO_HALF_FILLING", false)) {
    half_filling_ = hf;
    std::cout << "half-filling " << hf << " shifted to zero\n"; 
  }
  else 
    std::cout << "half-filling: " << hf << "\n"; 
  if (params.defined("W")) {
    bandwidth_scale_factor_ = bandwidth_/bw;
    std::cout << "bandwidth scale factor: " << bandwidth_scale_factor_ << "\n";
  }
  else {
    bandwidth_ = bw;
  }
  
  std::cout << "bandwidth: " << bandwidth_ << "\n";
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

int* DCATransformer::momentum_translation() const
{
	int *m = new int[n_cluster_sites()*n_cluster_sites()*n_cluster_sites()];
	std::valarray<double> phase(dimension());
	int res;
	for(site_t i=0; i<n_cluster_sites();++i){
		for(site_t j=0;j<n_cluster_sites();++j){
			for(site_t k=0;k<n_cluster_sites();++k){
				for(int d=0;d<dimension();++d){
					phase[d]=cluster_momentum(i,d)-cluster_momentum(j,d)+cluster_momentum(k,d);
				}
				for(site_t p=0;p<n_cluster_sites();++p){
					for(int d=0;d<dimension();++d){
						double tmp=(phase[d]-cluster_momentum(p,d))/(2.*M_PI);
						if(tmp==(int)tmp)
							res=1;
						else{
							res=0;
							break;
						}
					}
					if(res)
						m[i*n_cluster_sites()*n_cluster_sites()+j*n_cluster_sites()+k]=p;
				}
			}
		}
	}
	return m;
}

matsubara_green_function_t SimpleDCATransformer::operator()(const matsubara_green_function_t &G_omega, 
                                                            const matsubara_green_function_t &G0_omega, 
                                                            double &mu, const double h, const double beta)
{
  if(h!=0){std::cout<<"don't know how to handle magn. fields: h="<<h<<" exiting"<<std::endl; abort(); }
  unsigned int nflavor=G0_omega.nflavor();
  unsigned int nfreq=G0_omega.nfreq();
  assert(G0_omega.nsite()==n_cluster_sites());
  matsubara_green_function_t selfenergy(G0_omega);
  matsubara_green_function_t G0_omega_new(G0_omega);
  matsubara_green_function_t G_omega_lattice(G0_omega);
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
        selfenergy(freq, kc, kc, flavor)  = 1./G0_omega(freq, kc, kc, flavor);
        selfenergy(freq, kc, kc, flavor) -= 1./G_omega(freq, kc, kc, flavor);
      }
    }
  }
  
  bool converged=compute_filling_?false:true;
  bool damping_mu=false;
  bool computing_mu=true;
  int mu_it=1;
  double n_old;
  double mu_old;
  double mu_new;
  double mu_orig=mu;
  do{ //find the right particle density
    for (spin_t flavor=0; flavor<nflavor; ++flavor) {
      for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
	for (frequency_t freq=0; freq<nfreq; ++freq) {
	  G_omega_lattice(freq, kc, kc, flavor) = 0.;
          std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
	  for (site_t kl=0; kl<n_lattice_momenta(); ++kl) { 
	    vector_type k(dimension());
	    for (unsigned int d=0; d<dimension(); ++d) 
	      k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
	    double weight = lattice_momentum(kl,dimension());
	    G_omega_lattice(freq, kc, kc, flavor) += weight/(iomega_n+mu-epsilon(k)-selfenergy(freq, kc, kc, flavor));
	  }
	  G_omega_lattice(freq, kc, kc, flavor) /= n_cluster_cells();
	}
      }
    }
    std::cerr << "done" << std::endl;
    //compute density
    double density=0;
    for(spin_t flavor=0;flavor<nflavor;++flavor){
      for(site_t kc=0;kc<n_cluster_sites();++kc){
        for(frequency_t freq=0;freq<nfreq;++freq){
          //we could make this more accurate by using hi freq moments
          density+=2.*G_omega_lattice(freq,kc,kc,flavor).real();
        }
      }
    }
    density=(density/beta+0.5*nflavor*n_cluster_sites())/n_cluster_sites(); //density such that 'half filling' corresponds to 1
    if (compute_filling_) {
      std::cout<<"mu: "<<mu<<" filling "<<density<<" it: "<<mu_it<<" converged: "<<converged<<std::endl;
      //compute new chemical potential
      converged=adjust_chemical_potential(density, mu, damping_mu, computing_mu, mu_it, n_old, n_target_, mu_old, mu_new, mu_orig);
      std::cout<<"converged now: "<<converged<<" true is: "<<true<<" compute filling is: "<<compute_filling_<<std::endl;
    }
  } while(compute_filling_ && (!converged));
  (*const_cast<alps::Parameters*>(params_ptr_))["MU"]=mu; //write the chemical potential back to the parameter file
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
        G0_omega_new(freq, kc, kc, flavor) = 
	  1. / ( 1./G_omega_lattice(freq, kc, kc, flavor) + selfenergy(freq, kc, kc, flavor) );
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
      std::cout<<"WARNING: Input file could not be found, proceeding with free solution."<<std::endl;
  }
  { //proceed with computation of new G0 if we can't find input G0
    std::cout<<"calculating initial G0_omega"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    int n_tau=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    double beta = static_cast<double>(parms["BETA"]);
    double mu = static_cast<double>(parms["MU"]);
    //double h = static_cast<double>(parms["H"]);
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
    matsubara_green_function_t G0_omega_new = operator()(G_omega, G0_omega, mu, 0, beta);
    itime_green_function_t G0_tau(n_tau+1, n_cluster_sites(), n_flavor);
    boost::shared_ptr<FourierTransformer> fourier_ptr;
    G0_omega_new.write("G0_omega_input");
    FourierTransformer::generate_transformer(parms, fourier_ptr);
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
                                                                     const double beta, const int n_matsubara) const
{
  std::cout << "transforming cluster green function\n";
  assert(h==0); 
  unsigned int nz = G0_omega.nflavor();
  assert(G0_omega.nsite()==n_cluster_sites());
  matsubara_green_function_t G0_omega_fourpoint_inv(4*n_matsubara*n_matsubara, n_cluster_sites(), nz);
  G0_omega_fourpoint_inv.clear();
  for(spin_t z=0; z<nz; ++z) {
    int zbar = z%nz==0 ? z+1 : z-1;
    for(site_t k=0;k<n_cluster_sites();++k) {
      int kpq = k%2==0 ? k+1 : k-1;
      for (int om=0; om<2*n_matsubara; ++om) {
	std::complex<double> g_kpq_z;
	std::complex<double> g_k_zbar;
	if (om<n_matsubara) {
	  g_kpq_z = G_omega(om,kpq,kpq,z);
	  g_k_zbar = G_omega(om,k,k,zbar);
	} 
	else {
	  g_kpq_z = std::conj(G_omega(om-n_matsubara,kpq,kpq,z));
	  g_k_zbar = std::conj(G_omega(om-n_matsubara,k,k,zbar));
	}   
        G0_omega_fourpoint_inv(om+2*n_matsubara*om,k,k,z) = 1. / (g_kpq_z*g_k_zbar);
      } 
    }
  }
  matsubara_green_function_t G_omega_fourpoint_inv = fourpoint_green_function_inverse(G_omega_fourpoint, 2*n_matsubara);
  matsubara_green_function_t two_particle_vertex(4*n_matsubara*n_matsubara, n_cluster_sites(), nz);
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
  matsubara_green_function_t G_omega_fourpoint_new(4*n_matsubara*n_matsubara, n_cluster_sites(), nz);
  G_omega_fourpoint_new.clear();
  for (spin_t z=0; z<nz; ++z) {
    spin_t zbar = z%nz==0 ? z+1 : z-1;
    for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
      int kcpq = kc%2==0 ? kc+1 : kc-1;
      for (int om=0; om<2*n_matsubara; ++om) {
        std::complex<double> iomega_n;
	if(om<n_matsubara) 
	  iomega_n = std::complex<double>(0, (2.*om+1.)*M_PI/beta);
	else
	  iomega_n = std::complex<double>(0, -(2.*(om-n_matsubara)+1.)*M_PI/beta);
        for (site_t kl=0; kl<n_lattice_momenta(); ++kl) { 
	  vector_type k(dimension());
	  vector_type k_plus_q(dimension());
	  for (unsigned int d=0; d<dimension(); ++d) { 
	    k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
	    k_plus_q[d] = lattice_momentum(kl,d) + cluster_momentum(kcpq,d);
	  }
	  double weight = lattice_momentum(kl, dimension());
	  std::complex<double> s_kcpq_z;
	  std::complex<double> s_kc_zbar;
	  if(om<n_matsubara) {
	    s_kcpq_z = selfenergy(om, kcpq, kcpq, z);
	    s_kc_zbar = selfenergy(om, kc, kc, zbar);
	  }
	  else {
	    s_kcpq_z = std::conj(selfenergy(om-n_matsubara, kcpq, kcpq, z));
	    s_kc_zbar = std::conj(selfenergy(om-n_matsubara, kc, kc, zbar));
	  }
	  G_omega_fourpoint_new(om+2*n_matsubara*om, kc, kc, z) 
	    += weight*(1./(iomega_n+mu-epsilon(k_plus_q)-s_kcpq_z))
	    * (1./(iomega_n+mu-epsilon(k)-s_kc_zbar));
	}
	G_omega_fourpoint_new(om+2*n_matsubara*om, kc, kc, z) /= n_cluster_cells();
	G_omega_fourpoint_new(om+2*n_matsubara*om, kc, kc, z) = 1. / G_omega_fourpoint_new(om+2*n_matsubara*om, kc, kc, z);
      }
    }
  }
  G_omega_fourpoint_new = fourpoint_green_function_subtract(G_omega_fourpoint_new, two_particle_vertex);
  G_omega_fourpoint_new = fourpoint_green_function_inverse(G_omega_fourpoint_new, 2*n_matsubara);
  std::cout << "lattice green function calculated\n";
  return G_omega_fourpoint_new;
}



std::map<double, std::complex<double> > 
SimpleDCATransformer::local_green_function_from_selfenergies(const std::map<double, std::vector<std::complex<double> > >& selfenergies,
							     const double mu, const double delta)
{
  typedef std::map<double, std::vector<std::complex<double> > > selfenergy_t;
  std::map<double, std::complex<double> > local_G;
  for (selfenergy_t::const_iterator it=selfenergies.begin(); it!=selfenergies.end(); ++it) {
    for (int kc=0; kc<n_cluster_sites(); ++kc) { 
      for (int kl=0; kl<n_lattice_momenta(); ++kl) { 
	vector_type k(dimension());
	for (unsigned int d=0; d<dimension(); ++d) 
	  k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
	double weight = lattice_momentum(kl, dimension());
	std::complex<double> idelta(0., delta);
	if (local_G.find(it->first)==local_G.end()) 
	  local_G[it->first] = weight / (it->first + mu - epsilon(k) - it->second[kc] + idelta);
	else
	  local_G[it->first] += weight / (it->first + mu - epsilon(k) - it->second[kc] + idelta);
      }
    }
    local_G[it->first] /= n_cluster_cells()*n_cluster_sites();
  }
  return local_G;
}



alps::RealObsevaluator SimpleDCATransformer::kinetic_energy(const std::vector<std::vector<std::vector<alps::RealObsevaluator> > >& selfenergies, 
							    const double beta, const double mu, const bool error_propagation)
{
  const unsigned int nzone = selfenergies.size();
  const unsigned int nfreq = selfenergies[0][0].size()/2;
  if (error_propagation) {
    std::vector<alps::RealObsevaluator> energy_vec(nfreq);
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      double om = (2.*freq+1.)*M_PI/beta;
      for (int z=0; z<nzone; ++z) {
	for (int kc=0; kc<n_cluster_sites(); ++kc) { 
	  alps::RealObsevaluator Sigma_re = selfenergies[z][kc][2*freq];
	  alps::RealObsevaluator obsi = selfenergies[z][kc][2*freq+1];
	  Sigma_re -= mu;
	  obsi -= om;
	  alps::RealObsevaluator obsi_sq = obsi * obsi;
	  for (int kl=0; kl<n_lattice_momenta(); ++kl) { 
	    vector_type k(dimension());
	    for (unsigned int d=0; d<dimension(); ++d) 
	      k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
	    double weight = lattice_momentum(kl, dimension());
	    double eps = epsilon(k);
	    alps::RealObsevaluator obsr = Sigma_re;
	    obsr += eps;
	    obsr *= -1;
	    alps::RealObsevaluator obs = obsr;
	    obs *= obsr;
	    obs += obsi_sq;
	    obsr /= obs;
	    obsr *= eps;
	    // obsi /= obs;
	    // obsi *= eps;
	    if (z==0 && kc==0 && kl==0) {
	      energy_vec[freq] = weight*obsr;
	      // energy_vec[2*freq+1] = weight*obsi;
	    }
	    else {
	      energy_vec[freq] += weight*obsr;
	      // energy_vec[2*freq+1] += weight*obsi;
	    }
	  }
	}
      }
      energy_vec[freq] /= n_cluster_cells()*n_cluster_sites();
      // energy_vec[2*freq+1] /= n_cluster_cells()*n_cluster_sites();
    }
    std::cout << "dispmom: " << dispersion_moments(0) << "\t" << dispersion_moments(1) << "\n";
    std::ofstream test1str("test1.dat");
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      double om = (2.*freq+1.)*M_PI/beta;
      test1str << om << "\t" << energy_vec[freq].mean() << "\t"  << "\n";
    }
    test1str.close();
    alps::RealObsevaluator kin = energy_vec[0] + dispersion_moments(1)*nzone*beta*beta/M_PI/M_PI;
    for (frequency_t freq=1; freq<nfreq; ++freq) {
      double om = (2.*freq+1.)*M_PI/beta;
      kin += energy_vec[freq] + dispersion_moments(1)*nzone/om/om;
    }
    kin *= 2./beta;
    kin -= 0.25*nzone*dispersion_moments(1)*beta;
    return kin; 
  }
  else {
    std::vector<std::complex<double> > energy_test(nfreq, std::complex<double>(0., 0.));
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      std::complex<double> iom(0., (2.*freq+1.)*M_PI/beta);
      for (int z=0; z<nzone; ++z) {
	for (int kc=0; kc<n_cluster_sites(); ++kc) { 
	  std::complex<double> Sigma(selfenergies[z][kc][2*freq].mean(), selfenergies[z][kc][2*freq+1].mean());
	  for (int kl=0; kl<n_cluster_cells(); ++kl) { 
	    vector_type k(dimension());
	    for (unsigned int d=0; d<dimension(); ++d) 
	      k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
	    double weight = lattice_momentum(kl, dimension());
	    double eps = epsilon(k);
	    energy_test[freq] += weight*eps / ( iom+mu-eps-Sigma);
	  }
	}
      }
      energy_test[freq] /= n_cluster_cells()*n_cluster_sites();
    }
    std::ofstream test2str("test2.dat");
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      double om = (2.*freq+1.)*M_PI/beta;
      test2str << om << "\t" << energy_test[freq].real() << "\t"  << energy_test[freq].imag() << "\n";
    }
    test2str.close();
    double kin = 0.;
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      double om = (2.*freq+1.)*M_PI/beta;
      kin += energy_test[freq].real() + dispersion_moments(1)*nzone/om/om;
    }
    kin *= 2./beta;
    kin -= 0.25*nzone*dispersion_moments(1)*beta;
    alps::SimpleRealObservable kin_obs("kinetic energy");
    kin_obs.reset(true);
    kin_obs << kin;
    alps::RealObsevaluator kin_eval(kin_obs);
    return kin_eval; 
  }
}

matsubara_green_function_t AntiFerroDCATransformer::operator()(const matsubara_green_function_t &G_omega, 
                                                               const matsubara_green_function_t &G0_omega, 
                                                               double &mu, const double h, const double beta)
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
    for (int sign=-1; sign<=1; sign+=2) {
      matsubara_green_function_t G_kl(G0_omega);
      for (frequency_t freq=0; freq<nfreq; ++freq) {
	std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
	for (site_t kc=0; kc<n_cluster_sites(); kc+=2) {
	  vector_type k0(dimension());
	  for (unsigned int d=0; d<dimension(); ++d) 
	    k0[d] = sign*lattice_momentum(kl,d) + cluster_momentum(kc,d);
	  vector_type k1(dimension());
	  for (unsigned int d=0; d<dimension(); ++d) 
	    k1[d] = sign*lattice_momentum(kl,d) + cluster_momentum(kc+1,d);
	  for (spin_t flavor=0; flavor<nflavor; flavor+=2) { 
	    G_kl(freq,kc,kc,flavor) = iomega_n+mu-epsilon(k0);
	    G_kl(freq,kc,kc,flavor+1) = iomega_n+mu-epsilon(k0);
	    G_kl(freq,kc+1,kc+1,flavor) = iomega_n+mu-epsilon(k1);
	    G_kl(freq,kc+1,kc+1,flavor+1) = iomega_n+mu-epsilon(k1);
	    G_kl(freq,kc+1,kc,flavor) = -h/2;
	    G_kl(freq,kc+1,kc,flavor+1) = h/2;
	    G_kl(freq,kc,kc+1,flavor) = -h/2;
	    G_kl(freq,kc,kc+1,flavor+1) = h/2;
	  }
	}
      }
      G_kl = green_function_subtract(G_kl, selfenergy);
      G_kl = green_function_inverse(G_kl);
      G0_omega_new = green_function_add(G0_omega_new, G_kl);
    }
  }
  G0_omega_new = green_function_divide(G0_omega_new, std::complex<double>(2*n_cluster_cells()));
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
    std::cout<<"calculating initial G0_omega"<<std::endl;
    int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    int n_tau=boost::lexical_cast<int>(parms["N"]);
    int n_flavor=parms.value_or_default("FLAVORS", 2);
    double beta = static_cast<double>(parms["BETA"]);
    double mu = static_cast<double>(parms["MU"]);
    double h = static_cast<double>(parms["H"]);
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





std::map<double, ClusterHelper::Matrix<std::complex<double> >::Type> 
AntiFerroDCATransformer::local_green_function_from_selfenergies(const std::map<double, 
								std::vector<ClusterHelper::Matrix<std::complex<double> >::Type> >& /*selfenergies*/,
								const double /*mu*/, const double /*delta*/)
{
  boost::throw_exception(std::invalid_argument("to be implemented"));
  std::map<double, ClusterHelper::Matrix<std::complex<double> >::Type> dummy;
  return dummy;
}


alps::RealObsevaluator AntiFerroDCATransformer::kinetic_energy(const std::vector<std::vector<std::vector<alps::RealObsevaluator> > >& /*selfenergies*/,
							       const double /*beta*/, const double /*mu*/, const bool /*error_propagation*/)
{
  boost::throw_exception(std::invalid_argument("to be implemented"));
  return alps::RealObsevaluator();
}

				       



bool adjust_chemical_potential(const double &n_tot,
    double &mu,          //chemical potential, will be adjusted
    bool &damping_mu,    //do we use dampening?
    bool &computing_mu,  //are we computing the chemical potential?
    int &mu_it,          //counter for the iterations
    double &n_old,       //old density
    double &n_target,    //target density
    double &mu_old,      //old chemical pot
    double &mu_new,      //new chemical pot
    double &mu_orig      //original chemical potential
    )
{
  //static bool damping_mu=false;
  //static bool computing_mu=true;
  //static int mu_it=0;
  //double n_target=parms["N_TARGET"];
  //static double n_old, mu_old, mu_new, mu_orig=mu;
  //once we are converged.
  std::cout<<"running switch for: mu: "<<mu<<" mu old: "<<mu_old<<" mu orig: "
	   <<mu_orig<<" n_tot: "<<n_tot<<" n old: "<<n_old<<" n_target: "<<n_target<<std::endl;
  if(!damping_mu && (std::abs(n_tot-n_target)<1.e-5 || mu_it==10)){
    damping_mu=true;
    mu_new=0.7*mu_new +0.3*mu_orig; //that should damp oscillations a bit by weakening corrections
    std::cout<<"converged n/mu loop with mu "<<mu<<" applying damping: orig mu: "
	     <<mu_orig<<" new mu: "<<mu_new<<std::endl;
    n_old=n_tot;
    mu_old=mu;
    mu=mu_new;
  }
  //damping applied, now find the mu for the next iteration
  else if(damping_mu){
    n_old=n_tot;
    mu_old=mu;
    mu=mu_new;
    computing_mu=false;
  }
  //not yet converged, find the next chemical potential.
  else{
    computing_mu=true;
    if(mu_it==0){ //first simulation
      mu_new = mu; //let's try again, with the newly recomputed dc term
    }else if(mu_it==1){ //first simulation
      //that's quite a big change... radically tune mu such that we can do bisection.
      mu_new=n_tot<n_target?mu-0.1:mu+0.1;
      std::cout<<"first iteration. mu was: "<<mu<<" new mu is: "<<mu_new<<std::endl;
    }else if(n_tot < 0){
      //oops.. that means that we shifted the chemical potential up by too much
      std::cout<<"oops, get negative number of interacting particles!"<<std::endl;
      mu_new = (mu + mu_orig)/2; //this should bring us out of the interval where we have a negative particle number
    }else if((n_old< n_target && n_tot > n_target) || 
	     (n_old > n_target && n_tot < n_target)){ //target density is somewhere in the middle
      //bisection
      mu_new = (mu - mu_old)/(n_tot - n_old)*(n_target-n_old) + mu_old;
      std::cout<<"bisection: old mus were"<<mu_old<<" "<<mu<<" new mu is: "<<mu_new<<std::endl;
    }else{
      //we dropped out of the bisection interval :-( that implies a change of mu > 1.
      //let's limit it to one and be done with it.
      mu_new = (mu - mu_old)/(n_tot - n_old)*(n_target-n_old) + mu_old;
      std::cout<<"outside interval: new mu is"<<mu_new<<std::endl;
    }
    std::cout<<"new chemical potential is: "<<mu_new<<std::endl;
    n_old=n_tot;
    mu_old=mu;
    mu=mu_new;
    mu_it++;
  }
  return computing_mu?false:true;
}
