/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2007-2010 by Emanuel Gull <gull@phys.columbia.edu>
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
#include "alps/alea.h"
#include <boost/numeric/ublas/io.hpp>

bool adjust_chemical_potential(const double &density, double &mu, bool &damping_mu, 
                               bool &computing_mu, int &mu_it, double &n_old, double &n_target, 
                               double &mu_old, double &mu_new, double &mu_orig);



matsubara_green_function_t ClusterTransformer::initial_G0(const alps::Parameters& parms) {
  std::cout<<"reading initial G0_omega from: "<<parms["G0OMEGA_INPUT"]<<std::endl;
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
#ifndef USE_MPI
    std::cout<<"creating simple DCA transformer"<<std::endl;
#else
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank==0) std::cout<<"creating simple DCA transformer"<<std::endl;
#endif
    trans_ptr.reset(new SimpleDCATransformer(parms));
  }
}

void CDMFTTransformer::generate_transformer(const alps::Parameters& parms, boost::shared_ptr<CDMFTTransformer> &trans_ptr)
{
  bool antiferromagnetic = parms.value_or_default("ANTIFERROMAGNET", false);
  if (antiferromagnetic) {
    throw std::logic_error("Anti Ferro CDMFT is not implemented.");
  }
  else {
    std::cout<<"using simple CDMFT transformer"<<std::endl;
    if(static_cast<int>(parms["SITES"])!=4) throw std::runtime_error("only four-site CDMFT is implemented so far.");
    trans_ptr.reset(new FoursiteCDMFTTransformer(parms));
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
  //std::cout << "dispersion moments: " << dispersion_moments(0) << " " 
  //<< dispersion_moments(1) << std::endl;
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
  else{ 
    ; //std::cout << "half-filling: " << hf << "\n"; 
  }
  if (params.defined("W")) {
    bandwidth_scale_factor_ = bandwidth_/bw;
    std::cout << "bandwidth scale factor: " << bandwidth_scale_factor_ << "\n";
  }
  else {
    bandwidth_ = bw;
  }
#ifndef USE_MPI 
  std::cout << "bandwidth: " << bandwidth_ << "\n";
#else
  int rank;   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0) std::cout << "bandwidth: " << bandwidth_ << "\n";
#endif
}


std::complex<double>* DCATransformer::exp_iksite(const int sign) const
{
  std::complex<double> *m = new std::complex<double>[n_cluster_sites()*n_cluster_sites()];
  for (site_t k=0; k<n_cluster_sites(); ++k) {
    for (site_t i=0; i<n_cluster_sites(); ++i) {
      double phase = 0.;
      for (std::size_t d=0; d<dimension(); ++d) 
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
  int res=-1;
  for(site_t i=0; i<n_cluster_sites();++i){
    for(site_t j=0;j<n_cluster_sites();++j){
      for(site_t k=0;k<n_cluster_sites();++k){
        for(std::size_t d=0;d<dimension();++d){
          phase[d]=cluster_momentum(i,d)-cluster_momentum(j,d)+cluster_momentum(k,d);
        }
        for(site_t p=0;p<n_cluster_sites();++p){
          for(std::size_t d=0;d<dimension();++d){
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
  matsubara_green_function_t selfenergy(G_omega);
  return this->operator()(G_omega, G0_omega, selfenergy, mu, h, beta, false);
}

matsubara_green_function_t SimpleDCATransformer::operator()(const matsubara_green_function_t &G_omega, 
                                                            const matsubara_green_function_t &G0_omega, 
                                                            matsubara_green_function_t &selfenergy,
                                                            double &mu, const double h, const double beta, const bool use_selfenergy)
{
  if(h!=0) {
    throw std::invalid_argument("don't know how to handle magn. fields (use the AFM transformer for that!), exiting"); 
  }
  unsigned int nflavor=G0_omega.nflavor();
  unsigned int nfreq=G0_omega.nfreq();
  assert(G0_omega.nsite()==n_cluster_sites());
  matsubara_green_function_t G0_omega_new(G0_omega);
  matsubara_green_function_t G_omega_lattice(G0_omega);
  if(!use_selfenergy){ //compute selfenergy out of G and G0 first
    for (spin_t flavor=0; flavor<nflavor; ++flavor) {
      for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
        for (frequency_t freq=0; freq<nfreq; ++freq) {
          selfenergy(freq, kc, kc, flavor)  = 1./G0_omega(freq, kc, kc, flavor);
          selfenergy(freq, kc, kc, flavor) -= 1./G_omega(freq, kc, kc, flavor);
        }
      }
    }
  }
  bool converged=compute_filling_?false:true;
  bool damping_mu=false;
  bool computing_mu=true;
  int mu_it=1;

  int current_it = static_cast<int>( (*const_cast<alps::Parameters*>(params_ptr_)).value_or_default("CURRENT_IT",0));

  double n_old;
  double mu_old=mu;
  double mu_new=mu;
  double mu_orig=mu;
  std::vector<std::vector<double> >epsilon_kckl(n_cluster_sites(), std::vector<double>(n_lattice_momenta()));
  std::vector<double> weight_kl(n_lattice_momenta());
  for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
    for (site_t kl=0; kl<n_lattice_momenta(); ++kl) { 
      vector_type k(dimension());
      for (unsigned int d=0; d<dimension(); ++d)
        k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
      epsilon_kckl[kc][kl]=epsilon(k);
    }
  }for (site_t kl=0; kl<n_lattice_momenta(); ++kl) { 
    weight_kl[kl]= lattice_momentum(kl,dimension());
  }
  do{ //find the right particle density
    for (spin_t flavor=0; flavor<nflavor; ++flavor) {
      for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
        for (frequency_t freq=0; freq<nfreq; ++freq) {
          G_omega_lattice(freq, kc, kc, flavor) = 0.;
          std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
          for (site_t kl=0; kl<n_lattice_momenta(); ++kl) { 
            std::complex<double> zeta = iomega_n + mu + (flavor%2 ? -h : h);
            G_omega_lattice(freq, kc, kc, flavor) += weight_kl[kl]/(zeta-epsilon_kckl[kc][kl]-selfenergy(freq, kc, kc, flavor));
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
    density=density/beta/n_cluster_sites()+(0.5*n_cluster_sites())/n_cluster_sites()*nflavor; //density such that 'half filling' corresponds to 1
    if (compute_filling_ && current_it>1) {
      std::cout<<"mu: "<<mu<<" filling "<<density<<" it: "<<mu_it<<" converged: "<<converged<<std::endl;
      //compute new chemical potential
      converged=adjust_chemical_potential(density, mu, damping_mu, computing_mu, mu_it, n_old, n_target_, mu_old, mu_new, mu_orig);
      std::cout<<"converged now: "<<converged<<" true is: "<<true<<" compute filling is: "<<compute_filling_<<std::endl;
    }

	if (current_it<2){converged=true;}

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
  //proceed with computation of new G0 if we can't find input G0
  std::cout<<"calculating initial G0_omega"<<std::endl;
  unsigned int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  unsigned int n_tau=boost::lexical_cast<int>(parms["N"]);
  unsigned int n_flavor=parms.value_or_default("FLAVORS", 2);
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  //double h = static_cast<double>(parms["H"]);
  //bool paramagnet = static_cast<bool>(parms["PARAMAGNET"]);
  matsubara_green_function_t G0_omega_new(n_matsubara,n_cluster_sites(),n_flavor);
  G0_omega_new.clear();
  if (parms.value_or_default("INSULATING", false)) {
    std::cout <<"calculating insulating G0_omega";
    for (frequency_t freq=0; freq<n_matsubara; ++freq) {
      const std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      const std::complex<double> G0 = 1./(iomega_n+mu-dispersion_moments(0));
      for (spin_t flavor=0; flavor<n_flavor; ++flavor) {
        for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
          G0_omega_new(freq, kc, kc, flavor) = G0;
        }
      }
    }
  }
  else if (parms.value_or_default("ISOLATED", false)) {
    std::cout << "calculating G0_omega for isolated cluster"<<std::endl;
    for (frequency_t freq=0; freq<n_matsubara; ++freq) {
      for (spin_t flavor=0; flavor<n_flavor; ++flavor) {
        const std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
        for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
          vector_type k(dimension());
          for (unsigned int d=0; d<dimension(); ++d) 
            k[d] = cluster_momentum(kc,d);
          const std::complex<double> G0 = 1./(iomega_n+mu-epsilon(k));
          G0_omega_new(freq, kc, kc, flavor) = G0;
        }
      }
    }
  }
  else {
    matsubara_green_function_t G0_omega(n_matsubara,n_cluster_sites(),n_flavor);  
    matsubara_green_function_t G_omega(n_matsubara,n_cluster_sites(),n_flavor);
    G0_omega.clear();
    for (spin_t flavor=0; flavor<n_flavor; ++flavor) {
      for (site_t i=0; i<n_cluster_sites(); ++i) {
        for (frequency_t freq=0; freq<n_matsubara; ++freq) {
          G0_omega(freq, i, i, flavor) = 1.;
        }
      }
    }
    G_omega = G0_omega; 
    //invoke transformation with G_omega==G0_omega (selfenergy=0)
    G0_omega_new = operator()(G_omega, G0_omega, mu, 0, beta);
  }
  itime_green_function_t G0_tau(n_tau+1, n_cluster_sites(), n_flavor);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  G0_omega_new.write("G0_omega_input");
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega_new);
  G0_tau.write("G0_tau_input");
  matsubara_green_function_t G0_omega_real = transform_into_real_space(G0_omega_new);
  G0_omega_real.write("G0_omega_input_realspace");
  itime_green_function_t G0_tau_real = transform_into_real_space(G0_tau);
  G0_tau_real.write("G0_tau_input_realspace");
  //std::ofstream g0_omega_input_hr("g0_omega_input_hr.dat"); g0_omega_input_hr<<std::make_pair(G0_omega_new,beta);
  //std::ofstream g0_omega_input_rs_hr("g0_omega_input_rs_hr.dat"); g0_omega_input_rs_hr<<std::make_pair(G0_omega_real,beta);
  return G0_omega_new;
}




matsubara_green_function_t SimpleDCATransformer::transform_fourpoint(const matsubara_green_function_t& G_omega_fourpoint,
                                                                     const matsubara_green_function_t& G_omega,
                                                                     const matsubara_green_function_t& G0_omega, 
                                                                     const double mu, const double h, 
                                                                     const double beta, const unsigned int n_matsubara) const
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
      for (std::size_t om=0; om<2*n_matsubara; ++om) {
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
      for (std::size_t om=0; om<2*n_matsubara; ++om) {
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
    for (unsigned int kc=0; kc<n_cluster_sites(); ++kc) { 
      for (unsigned int kl=0; kl<n_lattice_momenta(); ++kl) { 
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
                                                            const double beta, const double mu, 
                                                            const bool error_propagation, const bool paramagnet)
{
  const unsigned int nzone = selfenergies.size();
  const unsigned int nfreq = selfenergies[0][0].size()/2;
  if (error_propagation) {
    std::cerr << "calculating the kinetic energy using error propagation\n"
    << "CAUTION: this may take a while\n";
    int del_z = paramagnet ? 2 : 1;
    std::vector<alps::RealObsevaluator> energy_vec(nfreq);
    std::vector<alps::RealObsevaluator> energy_kc(nfreq);
    for (unsigned int z=0; z<nzone; z+=del_z) {
      for (site_class_set_type::const_iterator it_kc=site_classes_k_space().first; 
           it_kc!=site_classes_k_space().second; ++it_kc) {
        const int multiplier = it_kc->size();
        const site_t kc = *((*it_kc).begin());
        vector_type kc_vec(dimension());
        for (unsigned int d=0; d<dimension(); ++d) 
          kc_vec[d] = cluster_momentum(kc,d);
        for (frequency_t freq=0; freq<nfreq; ++freq) {
          std::cerr << ".";
          const double om = (2.*freq+1.)*M_PI/beta;
          alps::RealObsevaluator obsi = selfenergies[z][kc][2*freq+1];
          obsi -= om;
          alps::RealObsevaluator obsi_sq = obsi * obsi;
          alps::RealObsevaluator mSigma_re = -selfenergies[z][kc][2*freq];
          mSigma_re += mu;      
          // kl == 0
          alps::RealObsevaluator obsr0 = mSigma_re;
          vector_type k0 = kc_vec;
          for (unsigned int d=0; d<dimension(); ++d) 
            k0[d] += lattice_momentum(0,d);
          const double weight0 = lattice_momentum(0, dimension());
          const double eps0 = epsilon(k0);
          obsr0 -= eps0;
          alps::RealObsevaluator obs0 = obsr0;
          obs0 *= obsr0;
          obs0 += obsi_sq;
          obsr0 /= obs0;
          obsr0 *= eps0;
          energy_kc[freq] = weight0*obsr0;
          for (site_t kl=1; kl<n_lattice_momenta(); ++kl) { 
            alps::RealObsevaluator obsr = mSigma_re;
            vector_type k = kc_vec;
            for (unsigned int d=0; d<dimension(); ++d) 
              k[d] += lattice_momentum(kl,d);
            const double weight = lattice_momentum(kl, dimension());
            const double eps = epsilon(k);
            obsr -= eps;
            alps::RealObsevaluator obs = obsr;
            obs *= obsr;
            obs += obsi_sq;
            obsr /= obs;
            obsr *= eps;
            energy_kc[freq] += weight*obsr;
          }
          energy_kc[freq] *= multiplier;
          if (it_kc==site_classes_k_space().first && z==0) 
            energy_vec[freq] = energy_kc[freq];
          else 
            energy_vec[freq] += energy_kc[freq];
        }
        std::cerr << "\n";
      }
    }
    for (frequency_t freq=0; freq<nfreq; ++freq) 
      energy_vec[freq] *= (double)del_z/(double)n_cluster_sites()/(double)n_cluster_cells();
    // frequency sum
    alps::RealObsevaluator kin = energy_vec[0] 
    + dispersion_moments(1)*nzone*beta*beta/M_PI/M_PI;
    for (frequency_t freq=1; freq<nfreq; ++freq) {
      double om = (2.*freq+1.)*M_PI/beta;
      kin += energy_vec[freq]  + dispersion_moments(1)*nzone/om/om;
    }
    kin *= 2./beta;
    kin -= 0.25*nzone*dispersion_moments(1)*beta;
    return kin; 
  }
  else {
    
    std::ofstream selfout("selftest.dat");
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      selfout << freq << " ";
      for (unsigned int kc=0; kc<n_cluster_sites(); ++kc) {
        for (unsigned int z=0; z<nzone; ++z) {
          selfout << selfenergies[z][kc][2*freq].mean() << " " <<  selfenergies[z][kc][2*freq+1].mean() << " ";
        }
      }
      selfout << "\n";
    }
    selfout.close();
    
    
    std::vector<std::complex<double> > energy_vec(nfreq, std::complex<double>(0., 0.));
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      const std::complex<double> iom(0., (2.*freq+1.)*M_PI/beta);
      for (unsigned int z=0; z<nzone; ++z) {
        for (unsigned int kc=0; kc<n_cluster_sites(); ++kc) {
          const std::complex<double> Sigma(selfenergies[z][kc][2*freq].mean(), selfenergies[z][kc][2*freq+1].mean());
          for (unsigned int kl=0; kl<n_lattice_momenta(); ++kl) { 
            vector_type k(dimension());
            for (unsigned int d=0; d<dimension(); ++d) 
              k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
            const double weight = lattice_momentum(kl, dimension());
            const double eps = epsilon(k);
            energy_vec[freq] += weight*eps / ( iom + mu -eps -Sigma );
          }
        }
      }
      energy_vec[freq] /= n_cluster_cells()*n_cluster_sites();
    }
    std::ofstream estr("etest.dat");
    double kin = 0.;
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      double om = (2.*freq+1.)*M_PI/beta;
      estr << om << " " << energy_vec[freq].real() << " "  << energy_vec[freq].imag() << "\n";
      kin += energy_vec[freq].real() + dispersion_moments(1)*nzone/om/om;
    }
    estr.close();
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
    unsigned int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
    unsigned int n_tau=boost::lexical_cast<int>(parms["N"]);
    unsigned int n_flavor=parms.value_or_default("FLAVORS", 2);
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
    fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega_new);
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
                                                               const double /*beta*/, const double /*mu*/, 
                                                               const bool /*error_propagation*/, const bool /*paramagnet*/)
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////FUNCTIONS FOR NAMBU SPACE//////////////////////////////////////////////////////////////////////////////////
void SimpleDCATransformer::construct_green_function(double mu, double beta, double eta1, double eta2, const std::vector<std::vector<double> > &epsilon_kckl, const std::vector<std::vector<double> > &dwave_prefactor_kckl, const std::vector<double> &weight_kl, const nambu_matsubara_green_function_t &selfenergy,nambu_matsubara_green_function_t &G_omega_lattice_kav){
  
  unsigned int nflavor=G_omega_lattice_kav.nflavor();
  unsigned int nfreq=G_omega_lattice_kav.nfreq();
  for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      G_omega_lattice_kav(freq, kc, kc, 0,0) = 0.;
      G_omega_lattice_kav(freq, kc, kc, 0,1) = 0.;
      G_omega_lattice_kav(freq, kc, kc, 1,0) = 0.;
      G_omega_lattice_kav(freq, kc, kc, 1,1) = 0.;
      std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      for (site_t kl=0; kl<n_lattice_momenta(); ++kl) { 
        //double prefactor=1.;//(k[0]==0 && k[1]==0)?1.:0.; //cos kx - cos ky
        double prefactor=dwave_prefactor_kckl[kc][kl];
        //this is the self consistency formula!!!
        //note the pairing field eta + i eta2 on the 0,1 component, eta - i eta2 on the 1,0 component.
        std::complex<double> G_omega_latt_inv[4]=
        {iomega_n + mu -epsilon_kckl[kc][kl]-selfenergy(freq, kc, kc, 0,0),
          -selfenergy(freq, kc, kc, 0,1) - prefactor*std::complex<double>( eta1, eta2),
          -selfenergy(freq, kc, kc, 1,0) - prefactor*std::complex<double>( eta1,-eta2),
          iomega_n - mu +epsilon_kckl[kc][kl]-selfenergy(freq, kc, kc, 1,1)};
        
        std::complex<double> G0omega_latt_inv_det=G_omega_latt_inv[0]*G_omega_latt_inv[3]-G_omega_latt_inv[1]*G_omega_latt_inv[2];
        G_omega_lattice_kav(freq, kc, kc, 0,0) += weight_kl[kl]*G_omega_latt_inv[3]/G0omega_latt_inv_det;
        G_omega_lattice_kav(freq, kc, kc, 0,1) -= weight_kl[kl]*G_omega_latt_inv[1]/G0omega_latt_inv_det;
        G_omega_lattice_kav(freq, kc, kc, 1,0) -= weight_kl[kl]*G_omega_latt_inv[2]/G0omega_latt_inv_det;
        G_omega_lattice_kav(freq, kc, kc, 1,1) += weight_kl[kl]*G_omega_latt_inv[0]/G0omega_latt_inv_det;
      }
      G_omega_lattice_kav(freq, kc, kc, 0,0) /= n_cluster_cells();
      G_omega_lattice_kav(freq, kc, kc, 0,1) /= n_cluster_cells();
      G_omega_lattice_kav(freq, kc, kc, 1,0) /= n_cluster_cells();
      G_omega_lattice_kav(freq, kc, kc, 1,1) /= n_cluster_cells();
    }
  }
}
double SimpleDCATransformer::compute_density(const nambu_matsubara_green_function_t &G_omega_lattice_kav, double beta){
  double density=0;
  double density_up_part=0;
  double density_dn_hole=0;
  //std::cout<<"density part and hole: "<<density_up_part<<" "<<density_dn_hole<<" "<<std::endl;
  for(site_t kc=0;kc<n_cluster_sites();++kc){
    for(frequency_t freq=0;freq<G_omega_lattice_kav.nfreq();++freq){
      //we could make this more accurate by using hi freq moments
      density_up_part+=2./beta*G_omega_lattice_kav(freq,kc,kc,0,0).real(); //particles
      density_dn_hole-=2./beta*G_omega_lattice_kav(freq,kc,kc,1,1).real(); //holes
    }
  }
  //compute density
  density_up_part=(density_up_part+n_cluster_sites()/2.0)/n_cluster_sites(); //density such that 'half filling' corresponds to 1
  density_dn_hole=(density_dn_hole+n_cluster_sites()/2.0)/n_cluster_sites(); //density such that 'half filling' corresponds to 1
  density=density_dn_hole+density_up_part;
  //std::cout<<"density total: "<<density<<std::endl;

  return density;
}
nambu_matsubara_green_function_t SimpleDCATransformer::operator()(const nambu_matsubara_green_function_t &G_omega, 
                                                                  const nambu_matsubara_green_function_t &G0_omega, 
                                                                  nambu_matsubara_green_function_t &selfenergy, 
                                                                  double &mu, const double h, const double eta1, const double eta2, const double beta, bool use_selfenergy)
{
  if(h!=0 ) { 
    throw std::invalid_argument("don't know how to handle magn. fields, exiting"); 
  }
  if(dimension()!=2 && (eta1 !=0 || eta2 !=0)) throw std::invalid_argument("not clear what to do with dwave symmetry away from dimension 2");
  unsigned int nflavor=G0_omega.nflavor();
  unsigned int nfreq=G0_omega.nfreq();
  nambu_matsubara_green_function_t G0_omega_new(G0_omega);
  nambu_matsubara_green_function_t G_omega_lattice_kav(G0_omega);
  
  std::vector<std::vector<double> >epsilon_kckl        (n_cluster_sites(), std::vector<double>(n_lattice_momenta()));
  std::vector<std::vector<double> >dwave_prefactor_kckl(n_cluster_sites(), std::vector<double>(n_lattice_momenta()));
  std::vector<double> weight_kl(n_lattice_momenta());
  for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
    for (site_t kl=0; kl<n_lattice_momenta(); ++kl) { 
      vector_type k(dimension());
      for (unsigned int d=0; d<dimension(); ++d)
        k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
      epsilon_kckl[kc][kl]=epsilon(k);
      dwave_prefactor_kckl[kc][kl]= cos(k[0]) - cos(k[1]); //cos kx - cos ky
    }
  }for (site_t kl=0; kl<n_lattice_momenta(); ++kl) { 
    weight_kl[kl]= lattice_momentum(kl,dimension());
  }

  if(!use_selfenergy){ //compute selfenergy out of G and G0 first
    selfenergy=G0_omega;
    //compute the lattice = cluster self energy
    for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
      for (frequency_t freq=0; freq<nfreq; ++freq) {
        std::complex<double> g0_block[4]={G0_omega(freq, kc, kc, 0,0),G0_omega(freq, kc, kc, 0,1),G0_omega(freq, kc, kc, 1,0),G0_omega(freq, kc, kc, 1,1)};
        std::complex<double>  g_block[4]={ G_omega(freq, kc, kc, 0,0), G_omega(freq, kc, kc, 0,1), G_omega(freq, kc, kc, 1,0), G_omega(freq, kc, kc, 1,1)};
        std::complex<double>      det_g0=g0_block[0]*g0_block[3]-g0_block[1]*g0_block[2];
        std::complex<double>      det_g = g_block[0]* g_block[3]- g_block[1]* g_block[2];
        
        selfenergy(freq,kc,kc,0,0)= g0_block[3]/det_g0 - g_block[3]/det_g;
        selfenergy(freq,kc,kc,0,1)=-g0_block[1]/det_g0 + g_block[1]/det_g;
        selfenergy(freq,kc,kc,1,0)=-g0_block[2]/det_g0 + g_block[2]/det_g;
        selfenergy(freq,kc,kc,1,1)= g0_block[0]/det_g0 - g_block[0]/det_g;
        
        //      std::cout<<"self energy is: "<<selfenergy(freq,kc,kc,0,0)<<" "<<selfenergy(freq,kc,kc,0,1)<<" "<<selfenergy(freq,kc,kc,1,0)<<" "<<selfenergy(freq,kc,kc,1,1)<<std::endl;
      }
    }
  }
  if(compute_filling_){
    double mu_min=-10.;
    double mu_max=10;
    double mu_upper=mu_max;
    double mu_lower=mu_min;
    double tol = 1.e-8;
    double density_upper=0., density_lower=0., density_now;
    
    std::cout<<"entering dyson lattice chem pot search."<<std::endl;
    //find upper density
    construct_green_function( mu_upper,  beta,  eta1,  eta2, epsilon_kckl, dwave_prefactor_kckl, weight_kl, selfenergy,G_omega_lattice_kav);
    density_upper=compute_density(G_omega_lattice_kav, beta);
    //find lower density
    construct_green_function( mu_lower,  beta,  eta1,  eta2, epsilon_kckl, dwave_prefactor_kckl, weight_kl, selfenergy,G_omega_lattice_kav);
    density_lower=compute_density(G_omega_lattice_kav, beta);
    
    //perform bisection to adjust chem. pot.
    do{
      mu=mu_lower+(mu_upper-mu_lower)/2.;
      construct_green_function( mu,  beta,  eta1,  eta2, epsilon_kckl, dwave_prefactor_kckl, weight_kl, selfenergy,G_omega_lattice_kav);
      density_now=compute_density(G_omega_lattice_kav, beta);
      if(density_now>n_target_){
        density_upper=density_now;
        mu_upper=mu;
        //std::cout<<"bisection found upper bound: "<<mu<<" "<<density_now<<std::endl;
      }else{
        density_lower=density_now;
        mu_lower=mu;
        //std::cout<<"bisection found lower bound: "<<mu<<" "<<density_now<<std::endl;
      }
    }while(std::abs(density_upper-density_lower)>tol);
    
    std::cout<<"bisection converged for mu: "<<mu<<" "<<density_lower<<" "<<density_upper<<" "<<density_now<<std::endl;
    (*const_cast<alps::Parameters*>(params_ptr_))["MU"]=mu; //write the chemical potential back to the parameter file
  }else{
    construct_green_function( mu,  beta,  eta1,  eta2, epsilon_kckl, dwave_prefactor_kckl, weight_kl, selfenergy,G_omega_lattice_kav);
  }
  
  
  
  
  //compute the bath Green's function G0
  for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      std::complex<double> glatt_block[4]={ G_omega_lattice_kav(freq, kc, kc, 0,0), G_omega_lattice_kav(freq, kc, kc, 0,1), G_omega_lattice_kav(freq, kc, kc, 1,0), G_omega_lattice_kav(freq, kc, kc, 1,1)};
      std::complex<double> det_glatt     = glatt_block[0]* glatt_block[3]- glatt_block[1]* glatt_block[2];
      std::complex<double> glatt_inv[4]  ={glatt_block[3]/det_glatt,-glatt_block[1]/det_glatt,-glatt_block[2]/det_glatt,glatt_block[0]/det_glatt};
      std::complex<double> g0_inv[4]     ={glatt_inv[0]+selfenergy(freq,kc,kc,0,0), glatt_inv[1]+selfenergy(freq,kc,kc,0,1),glatt_inv[2]+selfenergy(freq,kc,kc,1,0),glatt_inv[3]+selfenergy(freq,kc,kc,1,1)};
      std::complex<double> det_g0inv     =g0_inv[0]*g0_inv[3]-g0_inv[1]*g0_inv[2];
      //      std::cout<<" two dets: "<<det_glatt<<" "<<det_g0inv<<std::endl;
      G0_omega_new(freq, kc, kc, 0,0) = g0_inv[3]/det_g0inv;
      G0_omega_new(freq, kc, kc, 0,1) =-g0_inv[1]/det_g0inv;
      G0_omega_new(freq, kc, kc, 1,0) =-g0_inv[2]/det_g0inv;
      G0_omega_new(freq, kc, kc, 1,1) = g0_inv[0]/det_g0inv;
      //      std::cout<<"G0_omega_new: "<<G0_omega_new(freq,kc,kc,0,0)<<" "<<G0_omega_new(freq,kc,kc,0,1)<<" "<<G0_omega_new(freq,kc,kc,1,0)<<" "<<G0_omega_new(freq,kc,kc,1,1)<<std::endl;
    }
  }
  return G0_omega_new;
}

nambu_matsubara_green_function_t ClusterTransformer::initial_G0_nambu(const alps::Parameters& parms) {
  std::cout<<"reading initial (NAMBU) G0_omega from: "<<parms["G0OMEGA_INPUT"]<<std::endl;
  int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  int n_orbital=boost::lexical_cast<int>(parms["FLAVORS"]);
  int n_site=boost::lexical_cast<int>(parms["SITES"]);
  double beta =boost::lexical_cast<double>(parms["BETA"]);
  if(n_orbital!=2) throw std::invalid_argument("dunno how to handle nambu in >2 orbitals!");
  nambu_matsubara_green_function_t G0_omega(n_matsubara, n_site, n_orbital, diagonal);
  G0_omega.read(parms["G0OMEGA_INPUT"].c_str());
  return G0_omega;
}

nambu_matsubara_green_function_t SimpleDCATransformer::initial_G0_nambu(const alps::Parameters& parms)
{
  if (parms.defined("G0OMEGA_INPUT")){
    std::ifstream in_file(parms["G0OMEGA_INPUT"].c_str());
    if(in_file.good()){
      return ClusterTransformer::initial_G0_nambu(parms);
    }
    else
      std::cout<<"WARNING: Input file "<<parms["G0OMEGA_INPUT"]<<" defined but could not be found, proceeding with free solution."<<std::endl;
  }
  //proceed with computation of new G0 if we can't find input G0
  std::cout<<"calculating initial G0_omega"<<std::endl;
  unsigned int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  unsigned int n_tau=boost::lexical_cast<int>(parms["N"]);
  unsigned int n_flavor=parms.value_or_default("FLAVORS", 2);
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  
  double eta1_init = static_cast<double>(parms["ETA1_INIT"]);
  double eta2_init = static_cast<double>(parms["ETA2_INIT"]);
  double h_init = static_cast<double>(parms["H_INIT"]);
  
  nambu_matsubara_green_function_t G0_omega_new(n_matsubara,n_cluster_sites(),n_flavor, diagonal);
  nambu_matsubara_green_function_t selfenergy(n_matsubara,n_cluster_sites(),n_flavor, diagonal);
  G0_omega_new.clear();
  if (parms.defined("INSULATING")) {
    std::cout <<"calculating insulating G0_omega starting solution";
    for (frequency_t freq=0; freq<n_matsubara; ++freq) {
      const std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      const std::complex<double> G0 = 1./(iomega_n+mu-dispersion_moments(0));
      for (spin_t flavor=0; flavor<n_flavor; ++flavor) {
        for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
          G0_omega_new(freq, kc, kc, 0,0) = G0;
          G0_omega_new(freq, kc, kc, 1,1) = G0;
          G0_omega_new(freq, kc, kc, 0,1) = 0.;
          G0_omega_new(freq, kc, kc, 1,0) = 0.;
        }
      }
    }
  }else if (parms.value_or_default("ISOLATED", false)){
    for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
      for (frequency_t freq=0; freq<n_matsubara; ++freq) {
        G0_omega_new(freq, kc, kc, 0,0) = 0.;
        G0_omega_new(freq, kc, kc, 0,1) = 0.;
        G0_omega_new(freq, kc, kc, 1,0) = 0.;
        G0_omega_new(freq, kc, kc, 1,1) = 0.;
        std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
        
        vector_type k(dimension());
        for (unsigned int d=0; d<dimension(); ++d) 
          k[d] = cluster_momentum(kc,d);
        double dwave_prefactor= cos(k[0]) - cos(k[1]); //cos kx - cos ky
        double prefactor=dwave_prefactor;
        
        //note the pairing field eta + i eta2 on the 0,1 component, eta - i eta2 on the 1,0 component.
        std::complex<double> G0_omega_latt_inv[4]=
        {iomega_n + mu -epsilon(k),
          - prefactor*std::complex<double>( eta1_init, eta2_init),
          - prefactor*std::complex<double>( eta1_init,-eta2_init),
          iomega_n - mu +epsilon(k)};
        
        std::complex<double> G0omega_latt_inv_det=G0_omega_latt_inv[0]*G0_omega_latt_inv[3]-G0_omega_latt_inv[1]*G0_omega_latt_inv[2];
        G0_omega_new(freq, kc, kc, 0,0) =  G0_omega_latt_inv[3]/G0omega_latt_inv_det;
        G0_omega_new(freq, kc, kc, 0,1) =- G0_omega_latt_inv[1]/G0omega_latt_inv_det;
        G0_omega_new(freq, kc, kc, 1,0) =- G0_omega_latt_inv[2]/G0omega_latt_inv_det;
        G0_omega_new(freq, kc, kc, 1,1) =  G0_omega_latt_inv[0]/G0omega_latt_inv_det;
      }
    }
  }
  else {
    nambu_matsubara_green_function_t G0_omega(n_matsubara,n_cluster_sites(),n_flavor, diagonal);  
    nambu_matsubara_green_function_t G_omega(n_matsubara,n_cluster_sites(),n_flavor, diagonal);
    G0_omega.clear();
    for (site_t i=0; i<n_cluster_sites(); ++i) {
      for (frequency_t freq=0; freq<n_matsubara; ++freq) {
        G0_omega(freq, i, i, 0,0) = 1.;
        G0_omega(freq, i, i, 0,1) = 0.;
        G0_omega(freq, i, i, 1,0) = 0.;
        G0_omega(freq, i, i, 1,1) = 1.;
      }
    }
    G_omega = G0_omega; 
    //invoke transformation with G_omega==G0_omega (selfenergy=0)
    std::cout<<"invoking free transform in initial."<<std::endl;
    G0_omega_new = operator()(G_omega, G0_omega, selfenergy, mu, h_init, eta1_init, eta2_init, beta, false);
  }
  nambu_itime_green_function_t G0_tau(n_tau+1, n_cluster_sites(), n_flavor, diagonal);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  G0_omega_new.write("G0_omega_input");
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega_new);
  G0_tau.write("G0_tau_input");
  //nambu_matsubara_green_function_t G0_omega_real = transform_into_real_space(G0_omega_new);
  //G0_omega_real.write("G0_omega_input_realspace");
  //nambu_itime_green_function_t G0_tau_real = transform_into_real_space(G0_tau);
  //G0_tau_real.write("G0_tau_input_realspace");
  //std::ofstream g0_omega_input_hr("g0_omega_input_hr.dat"); g0_omega_input_hr<<std::make_pair(G0_omega_new,beta);
  //std::ofstream g0_omega_input_rs_hr("g0_omega_input_rs_hr.dat"); g0_omega_input_rs_hr<<std::make_pair(G0_omega_real,beta);
  return G0_omega_new;
}
nambu_matsubara_green_function_t AntiFerroDCATransformer::operator()(const nambu_matsubara_green_function_t &G_omega, 
                                                                     const nambu_matsubara_green_function_t &G0_omega, 
                                                                     nambu_matsubara_green_function_t &selfenergy, 
                                                                     double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy){
  throw std::invalid_argument("antiferro & nambu DCA has not yet been implemented!");
}
nambu_green_function<double> SimpleDCATransformer::spin_symmetrizer(const nambu_green_function<double>& G, const bool paramagnet) const
{
  nambu_green_function<double> G_new(G);
  if (paramagnet) {
    for (site_t k1=0; k1<n_cluster_sites(); ++k1) {
      for (site_t k2=0; k2<n_cluster_sites(); ++k2) {
        for(itime_index_t n=0; n<G_new.nfreq(); ++n) {
          itime_index_t other_n=G_new.nfreq()-1-n;
          G_new(n, k1, k2, 0,0)= 0.5*(G(n, k1, k2, 0,0)+G(other_n, k1, k2, 1,1));
          G_new(other_n, k1, k2, 1,1)= G_new(n, k1, k2, 0,0);
          G_new(n, k1, k2, 0,1)= 0.5*(G(n, k1, k2, 0,1)-G(other_n, k1, k2, 1,0));
          G_new(other_n, k1, k2, 1,0)=-G_new(n, k1, k2, 0,1);
        }
        //std::cout<<"symmetrizing: "<<G(0, k1, k2, 0,0)<<" "<<G(G_new.nfreq()-1, k1, k2, 1,1)<<" and: "<<G(0, k1, k2, 0,1)<<" "<<G(G_new.nfreq()-1, k1, k2, 1,0)<<std::endl;
        //std::cout<<"new values  : "<<G_new(0, k1, k2, 0,0)<<" "<<G_new(G_new.nfreq()-1, k1, k2, 1,1)<<" and: "<<G_new(0, k1, k2, 0,1)<<" "<<G_new(G_new.nfreq()-1, k1, k2, 1,0)<<std::endl;
      }
    }
  }else {
    throw std::logic_error("do Nambu & AFM!!");
  }
  return G_new;
}
nambu_green_function<std::complex<double> > SimpleDCATransformer::spin_symmetrizer(const nambu_green_function<std::complex<double> >& G, const bool paramagnet) const
{
  nambu_green_function<std::complex<double> > G_new(G);
  //std::ofstream nambu_before_file("beforespin.dat");
  //nambu_before_file<<std::make_pair(G, 10.)<<std::endl;;
  if (paramagnet) {
    for (site_t k1=0; k1<n_cluster_sites(); ++k1) {
      for (site_t k2=0; k2<n_cluster_sites(); ++k2) {
        for(frequency_t freq=0; freq<G_new.nfreq(); ++freq) {
          G_new(freq, k1, k2, 0,0).imag()= 0.5*(G(freq, k1, k2, 0,0).imag()+G_new(freq, k1, k2, 1,1).imag());
          G_new(freq, k1, k2, 0,0).real()= 0.5*(G(freq, k1, k2, 0,0).real()-G_new(freq, k1, k2, 1,1).real());
          G_new(freq, k1, k2, 1,1).real()= -G_new(freq, k1, k2, 0,0).real();
          G_new(freq, k1, k2, 1,1).imag()=  G_new(freq, k1, k2, 0,0).imag();
          //symmetrize the  F.
          G_new(freq, k1, k2, 0,1)= (G(freq, k1, k2, 0,1)+std::conj(G(freq, k1, k2, 1,0)))/2.;
          G_new(freq, k1, k2, 1,0)= std::conj(G_new(freq, k1, k2, 0,1));
        }
        //std::cout<<"symmetrizing: "<<G(0, k1, k2, 0,0)<<" "<<G(0, k1, k2, 1,1)<<" and: "<<G(0, k1, k2, 0,1)<<" "<<G(0, k1, k2, 1,0)<<std::endl;
        //std::cout<<"new values  : "<<G_new(0, k1, k2, 0,0)<<" "<<G_new(0, k1, k2, 1,1)<<" and: "<<G_new(0, k1, k2, 0,1)<<" "<<G_new(0, k1, k2, 1,0)<<std::endl;
      }
    }
  }
  //std::ofstream nambu_after_file("afterspin.dat");
  //nambu_after_file<<std::make_pair(G_new, 10.);
  return G_new;
}

/*
 *
 *
 *
 These are all foursite CDMFT things. relatively easily extendable to other four-site things...
 *
 *
 *
 */

//Auxiliary functions
//invert a matrix
void dense_matrix_inverse(complex_dense_matrix &M /*input*/, complex_dense_matrix &M_inv/*output*/) 
{
  M_inv = boost::numeric::ublas::identity_matrix<std::complex<double> >(M.size1());
  boost::numeric::ublas::vector<fortran_int_t> ipivot(M.size1());
  //memset(&(ipivot[0]), 0, sizeof(fortran_int_t)*M.size1());
  boost::numeric::bindings::lapack::gesv(M,ipivot, M_inv);
}

//write a gf to a matrix
template<class T> void gf_to_matrix(int omega, int flavor, const green_function<T> &gf, boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> &M){
  for(unsigned int i=0;i<gf.nsite();++i){
    for(unsigned int j=0;j<gf.nsite();++j){
      M(i,j)=gf(omega, i,j,flavor);
    }
  }
}
//write a matrix to a gf
template<class T> void matrix_to_gf(int omega, int flavor,  green_function<T> &gf, const boost::numeric::ublas::matrix<T,boost::numeric::ublas::column_major> &M){
  for(unsigned int i=0;i<gf.nsite();++i){
    for(unsigned int j=0;j<gf.nsite();++j){
      gf(omega, i,j,flavor)=M(i,j);
    }
  }
}

CDMFTTransformer::CDMFTTransformer(const alps::Parameters& params) 
: ClusterTransformer(params), 
t_(params.defined("W") ? 1. : static_cast<double>(params["t"])), 
t2_(params.value_or_default("t'", params.value_or_default("t2", 0.))),
t3_(params.value_or_default("t3", 0.)),
t4_(params.value_or_default("t4", 0.)),
bandwidth_(params.value_or_default("W", 0.)),
bandwidth_scale_factor_(1.),
half_filling_(0.),
L_(params["L"]),
dispersion_moments_(2, 0.)
{
  if (params.defined("t") && params.defined("W"))
    boost::throw_exception(std::invalid_argument("both t and W are defined"));
  //calculating dispersion moments, bandwidth and half filling
  std::cout<<"warning: the CDMFT transformer has NO bandwidths, dispersions, etc implemented.";
}
matsubara_green_function_t FoursiteCDMFTTransformer::operator()(const matsubara_green_function_t &G_omega, 
                                                                const matsubara_green_function_t &G0_omega, 
                                                                double &mu, const double h, const double beta)
{
  if(h!=0) { // && G0_omega.nsite()!=1 <-this not understood - looks like a bug to me.
    std::cerr<<"H is: "<<h<<std::endl;
    throw std::invalid_argument("don't know how to handle magn. fields (use the AFM transformer for that!), exiting"); 
  }
  unsigned int nflavor=G0_omega.nflavor();
  unsigned int nfreq=G0_omega.nfreq();
  if(G0_omega.nsite()!=n_cluster_sites()) throw std::runtime_error("problem with cluster sites");
  matsubara_green_function_t selfenergy(G0_omega);
  matsubara_green_function_t G0_omega_new(G0_omega);
  matsubara_green_function_t G_omega_lattice(G0_omega);
  
  std::cout<<"computing self energy."<<std::endl;
  //compute the Matsubara self energy
  complex_dense_matrix G0_mat(n_cluster_sites(),n_cluster_sites()), G_mat(n_cluster_sites(),n_cluster_sites()), Sigma_mat(n_cluster_sites(),n_cluster_sites());
  complex_dense_matrix G0_inv_mat(n_cluster_sites(),n_cluster_sites()), G_inv_mat(n_cluster_sites(),n_cluster_sites());
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (frequency_t omega=0; omega<nfreq; ++omega) {
      gf_to_matrix(omega, flavor, G0_omega, G0_mat);
      gf_to_matrix(omega, flavor, G_omega, G_mat);
      dense_matrix_inverse(G0_mat, G0_inv_mat);
      dense_matrix_inverse(G_mat, G_inv_mat);
      Sigma_mat=G0_inv_mat-G_inv_mat;
      matrix_to_gf(omega, flavor, selfenergy, Sigma_mat);
    }
  }
  std::cout<<"computing k summation."<<std::endl;
  
  selfenergy = transform_into_real_space(selfenergy); //the integration is done in real space
  
  //the k-summation
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (frequency_t freq=0; freq<nfreq; ++freq) {
      std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      complex_dense_matrix diag=boost::numeric::ublas::identity_matrix<std::complex<double> >(n_cluster_sites())*(iomega_n+mu);
      complex_dense_matrix H0(n_cluster_sites(), n_cluster_sites());
      //complex_dense_matrix H0K(n_cluster_sites(), n_cluster_sites());
      complex_dense_matrix tmp(n_cluster_sites(), n_cluster_sites());
      complex_dense_matrix Ginv(n_cluster_sites(), n_cluster_sites());
      G_mat=boost::numeric::ublas::zero_matrix<std::complex<double> >(n_cluster_sites(),n_cluster_sites());
      H0   =boost::numeric::ublas::zero_matrix<std::complex<double> >(n_cluster_sites(),n_cluster_sites());
      tmp  =boost::numeric::ublas::zero_matrix<std::complex<double> >(n_cluster_sites(),n_cluster_sites());
      gf_to_matrix(freq, flavor, selfenergy, Sigma_mat);
      for(int nx=0;nx<L_;++nx){
        double kx=-M_PI+2*M_PI*nx/L_;
        for(int ny=0;ny<L_;++ny){
          //transformation of the H0 - matrix to k-space
          double ky=-M_PI+2*M_PI*ny/L_;
          set_Tij_matrix(kx,ky,H0);
          //Ginv = iomega +mu - epsilon - Sigma
          Ginv=diag+H0+Sigma_mat;
          dense_matrix_inverse(Ginv,tmp);
          G_mat+=tmp/(L_*L_);
        }
      }
      matrix_to_gf(freq, flavor, G_omega_lattice, G_mat);
    }
  }     
  std::cout<<"lattice GF is: "<<std::endl;
  std::ofstream lgf_file("lgf.dat");
  std::cout<<"in real space: "<<std::endl;
  lgf_file<<std::make_pair(G_omega_lattice,beta)<<std::endl;
  std::ofstream lgfks_file("lgfks.dat");
  std::cout<<"in k space: "<<std::endl;
  G_omega_lattice = transform_into_k_space(G_omega_lattice);
  lgfks_file<<std::make_pair(G_omega_lattice,beta)<<std::endl;
  
  std::cout<<"computing G0 from sigma"<<std::endl;
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (frequency_t omega=0; omega<nfreq; ++omega) {
      gf_to_matrix(omega, flavor, selfenergy, Sigma_mat);
      gf_to_matrix(omega, flavor, G_omega_lattice, G_mat);
      dense_matrix_inverse(G_mat, G_inv_mat);
      G0_inv_mat=Sigma_mat+G_inv_mat;
      dense_matrix_inverse(G0_inv_mat, G0_mat);
      matrix_to_gf(omega, flavor, G0_omega_new, G0_mat);
    }
  }
  std::cout<<"returning from CDMFT self consistency"<<std::endl;
  std::ofstream g0newks_file("g0newks.dat");
  g0newks_file<<std::make_pair(G0_omega_new,beta)<<std::endl;
  return G0_omega_new;
}



matsubara_green_function_t FoursiteCDMFTTransformer::initial_G0(const alps::Parameters& parms)
{
  if (parms.defined("G0OMEGA_INPUT")){
    std::ifstream in_file(parms["G0OMEGA_INPUT"].c_str());
    if(in_file.good()){
      return ClusterTransformer::initial_G0(parms);
    }
    else
      std::cout<<"WARNING: Input file could not be found, proceeding with free solution."<<std::endl;
  }
  //proceed with computation of new G0 if we can't find input G0
  std::cout<<"calculating initial G0_omega"<<std::endl;
  unsigned int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  unsigned int n_tau=boost::lexical_cast<int>(parms["N"]);
  unsigned int n_flavor=parms.value_or_default("FLAVORS", 2);
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  //double h = static_cast<double>(parms["H"]);
  //bool paramagnet = static_cast<bool>(parms["PARAMAGNET"]);
  matsubara_green_function_t G0_omega_new(n_matsubara,n_cluster_sites(),n_flavor);
  G0_omega_new.clear();
  if (parms.value_or_default("INSULATING", false)) {
    std::cout <<"calculating insulating G0_omega";
    throw std::runtime_error("insulating g0 in foursite transformer not implemented.");
  }
  else if (parms.value_or_default("ISOLATED", false)) {
    std::cout << "calculating G0_omega for isolated cluster"<<std::endl;
    throw std::runtime_error("isolated g0 in foursite transformer not implemented.");
  }
  else {
    matsubara_green_function_t G0_omega(n_matsubara,n_cluster_sites(),n_flavor);  
    matsubara_green_function_t G_omega(n_matsubara,n_cluster_sites(),n_flavor);
    G0_omega.clear();
    for (spin_t flavor=0; flavor<n_flavor; ++flavor) {
      for (site_t i=0; i<n_cluster_sites(); ++i) {
        for (frequency_t freq=0; freq<n_matsubara; ++freq) {
          G0_omega(freq, i, i, flavor) = 1.;
        }
      }
    }
    G_omega = G0_omega; 
    //invoke transformation with G_omega==G0_omega (selfenergy=0)
    G0_omega_new = operator()(G_omega, G0_omega, mu, 0, beta);
  }
  itime_green_function_t G0_tau(n_tau+1, n_cluster_sites(), n_flavor);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  G0_omega_new.write("G0_omega_input");
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega_new);
  G0_tau.write("G0_tau_input");
  matsubara_green_function_t G0_omega_real = transform_into_real_space(G0_omega_new);
  G0_omega_real.write("G0_omega_input_realspace");
  itime_green_function_t G0_tau_real = transform_into_real_space(G0_tau);
  G0_tau_real.write("G0_tau_input_realspace");
  std::ofstream g0_omega_input_hr("g0_omega_input_hr.dat"); g0_omega_input_hr<<std::make_pair(G0_omega_new,beta);
  std::ofstream g0_omega_input_rs_hr("g0_omega_input_rs_hr.dat"); g0_omega_input_rs_hr<<std::make_pair(G0_omega_real,beta);
  return G0_omega_new;
}




matsubara_green_function_t FoursiteCDMFTTransformer::transform_fourpoint(const matsubara_green_function_t& G_omega_fourpoint,
                                                                         const matsubara_green_function_t& G_omega,
                                                                         const matsubara_green_function_t& G0_omega, 
                                                                         const double mu, const double h, 
                                                                         const double beta, const unsigned int n_matsubara) const
{
  throw std::runtime_error("transform_fourpoint in foursite transformer not implemented.");
  
  return G_omega_fourpoint;
}



std::map<double, std::complex<double> > 
FoursiteCDMFTTransformer::local_green_function_from_selfenergies(const std::map<double, std::vector<std::complex<double> > >& selfenergies,
                                                                 const double mu, const double delta)
{
  throw std::runtime_error("local green to selfenergies in foursite transformer not implemented.");
  return std::map<double, std::complex<double> >();
}



alps::RealObsevaluator FoursiteCDMFTTransformer::kinetic_energy(const std::vector<std::vector<std::vector<alps::RealObsevaluator> > >& selfenergies, 
                                                                const double beta, const double mu, 
                                                                const bool error_propagation, const bool paramagnet)
{
  throw std::runtime_error("foursite CDMFT transformer not implemented.");
}
nambu_matsubara_green_function_t FoursiteCDMFTTransformer::initial_G0_nambu(alps::Parameters const&){
  throw std::runtime_error("initial nambu_matsubara_green_function_t CDMFT transformer not implemented.");
}
nambu_matsubara_green_function_t FoursiteCDMFTTransformer::operator()(const nambu_matsubara_green_function_t &G_omega, 
                                                                      const nambu_matsubara_green_function_t &G0_omega, 
                                                                      nambu_matsubara_green_function_t &selfenergy, 
                                                                      double &mu, const double h, const double eta1, const double eta2, const double beta, const bool use_selfenergy){
  throw std::runtime_error("nambu operator() in CDMFT transformer not implemented.");
}
/*nambu_green_function<double> FoursiteCDMFTTransformer::transform_into_k_space(nambu_green_function<double> const& gf) const{
 throw std::runtime_error("transform_into_k_space CDMFT transformer not implemented.");
 }
 nambu_green_function<std::complex<double> > FoursiteCDMFTTransformer::transform_into_k_space(nambu_green_function<std::complex<double> > const& gf) const{
 throw std::runtime_error("transform_into_k_space CDMFT transformer not implemented.");
 }*/
matsubara_green_function_t FoursiteCDMFTTransformer::transform_into_real_space(const matsubara_green_function_t& G) const{
  spin_t nflavor=G.nflavor();
  frequency_t nfreq=G.nfreq();
  matsubara_green_function_t G_transformed(G);
  complex_dense_matrix U=U_;
  complex_dense_matrix tmp=U_;
  complex_dense_matrix G_mat(n_cluster_sites(),n_cluster_sites()), G_mat_new(n_cluster_sites(),n_cluster_sites());
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (frequency_t omega=0; omega<nfreq; ++omega) {
      gf_to_matrix(omega, flavor, G, G_mat);
      tmp=boost::numeric::ublas::prod(G_mat,U);
      G_mat_new=boost::numeric::ublas::prod(U,tmp);
      matrix_to_gf(omega, flavor, G_transformed, G_mat_new);
    }
  }
  return G_transformed;
}
matsubara_green_function_t FoursiteCDMFTTransformer::transform_into_k_space(const matsubara_green_function_t& G) const{
  return transform_into_real_space(G); //this is a special foursite case.
}
itime_green_function_t FoursiteCDMFTTransformer::transform_into_real_space(const itime_green_function_t& G) const{
  spin_t nflavor=G.nflavor();
  itime_index_t ntime=G.ntime();
  itime_green_function_t G_transformed(G);
  dense_matrix U(U_.size1(), U_.size1());
  for(unsigned int i=0;i<U_.size1();++i){ for (unsigned int j=0;j<U_.size1();++j){ U(i,j)=U_(i,j).real();} }
  dense_matrix tmp=U;
  dense_matrix G_mat(n_cluster_sites(),n_cluster_sites()), G_mat_new(n_cluster_sites(),n_cluster_sites());
  for (spin_t flavor=0; flavor<nflavor; ++flavor) {
    for (itime_index_t tau=0; tau<ntime; ++tau) {
      gf_to_matrix(tau, flavor, G, G_mat);
      tmp=boost::numeric::ublas::prod(G_mat,U);
      G_mat_new=boost::numeric::ublas::prod(U,tmp);
      matrix_to_gf(tau, flavor, G_transformed, G_mat_new);
    }
  }
  return G_transformed;
}  
itime_green_function_t FoursiteCDMFTTransformer::transform_into_k_space(const itime_green_function_t& G) const{
  return transform_into_real_space(G);
}
nambu_matsubara_green_function_t FoursiteCDMFTTransformer::transform_into_real_space(const nambu_matsubara_green_function_t& G) const{
  throw std::runtime_error("this has not yet been implemented.");
}
nambu_matsubara_green_function_t FoursiteCDMFTTransformer::transform_into_k_space(const nambu_matsubara_green_function_t& G) const{
  throw std::runtime_error("this has not yet been implemented.");
}
nambu_itime_green_function_t FoursiteCDMFTTransformer::transform_into_real_space(const nambu_itime_green_function_t& G) const{
  throw std::runtime_error("this has not yet been implemented.");
}
nambu_itime_green_function_t FoursiteCDMFTTransformer::transform_into_k_space(const nambu_itime_green_function_t& G) const{
  throw std::runtime_error("this has not yet been implemented.");
}
void FoursiteCDMFTTransformer::set_Tij_matrix(double kx, double ky, complex_dense_matrix &H0){
  H0(0,0)=0;  H0(1,1)=0; H0(2,2)=0; H0(3,3)=0;
  double coskx=cos(kx);
  double cosky=cos(ky);
  
  std::complex<double> exp_ikx=std::exp(std::complex<double>(0, kx));
  std::complex<double> exp_mikx=std::conj(exp_ikx);
  std::complex<double> exp_iky=std::exp(std::complex<double>(0, ky));
  std::complex<double> exp_miky=std::conj(exp_iky);
  //these are the factors of 't'
  H0(0,1)=-2.*t_*coskx*exp_ikx;
  H0(0,3)=-2.*t_*cosky*exp_miky;
  H0(2,1)=-2.*t_*cosky*exp_iky;
  H0(2,3)=-2.*t_*coskx*exp_mikx;
  H0(1,0)=-2.*t_*coskx*exp_mikx;
  H0(1,2)=-2.*t_*cosky*exp_miky;
  H0(3,0)=-2.*t_*cosky*exp_iky;
  H0(3,2)=-2.*t_*coskx*exp_ikx;
  //these are the factors of 'tprime'
  H0(0,2)=-4.*std::exp(std::complex<double>(0, ( kx - ky)))*tprime_*coskx*cosky;
  H0(1,3)=-4.*std::exp(std::complex<double>(0, (-kx - ky)))*tprime_*coskx*cosky;
  H0(0,2)=-4.*std::exp(std::complex<double>(0, (-kx + ky)))*tprime_*coskx*cosky;
  H0(0,2)=-4.*std::exp(std::complex<double>(0, ( kx + ky)))*tprime_*coskx*cosky;
}
std::vector<double> FoursiteCDMFTTransformer::epsav() const{
  return epsav_;
}
std::vector<double> FoursiteCDMFTTransformer::epssqav() const{
  return epssqav_;
}


