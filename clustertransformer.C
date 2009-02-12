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


DCATransformer::DCATransformer(const alps::Parameters& params) 
  : ClusterTransformer(params), 
    t_(boost::lexical_cast<double>(params["t"])), 
    t2_(params.value_or_default("t'", params.value_or_default("t2", 0.))),
    t3_(params.value_or_default("t3", 0.)),
    t4_(params.value_or_default("t4", 0.)),
    epsav_(n_cluster_sites()), 
    epssqav_(n_cluster_sites())
{
  //calculating dispersion moments, bandwidth and half filling
  std::vector<double> energies;
  for (site_t kc=0; kc<n_cluster_sites(); ++kc) {
    for (site_t kl=0; kl<n_cluster_cells(); ++kl) {
      vector_type k(dimension());
      for (unsigned int d=0; d<dimension(); ++d)
        k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
      double E = epsilon(k);
      epsav_  [kc]+= E;
      epssqav_[kc]+= E*E;
      energies.push_back(E); 
    }    
    epsav_  [kc] /= n_cluster_cells();
    epssqav_[kc] /= n_cluster_cells();
  }
  if (t2_==0. && t3_==0. && t4_==0.) {
    bandwidth_ = dimension()*4*t_;
    half_filling_ = 0;
  }
  else {
    std::sort(energies.begin(), energies.end());
    double min = energies.front();
    double max = energies.back();
    bandwidth_ = max - min;
    if (energies.size()%2) 
      half_filling_ = energies[energies.size()/2];
    else 
      half_filling_ = (energies[energies.size()/2-1] + energies[energies.size()/2])/2. ;
  }
  std::cout << "Bandwidth: " << bandwidth_ << ", half filling at: " << half_filling_ << "\n";
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
	        for (site_t kl=0; kl<n_cluster_cells(); ++kl) { 
	          vector_type k(dimension());
	          for (unsigned int d=0; d<dimension(); ++d) 
	            k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
	          G_omega_lattice(freq, kc, kc, flavor) += 1./(iomega_n+mu-epsilon(k)-selfenergy(freq, kc, kc, flavor));
          }
          G_omega_lattice(freq, kc, kc, flavor) /= n_cluster_cells();
        }
      }
    }
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
    FourierTransformer::generate_transformer(parms, fourier_ptr, epsav(), epssqav());
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
  std::cout << "transforming cluster green function\n";
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
  std::cout << "lattice green function calculated\n";
  //exit(1);
  return G_omega_fourpoint_new;
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
    std::cout<<"calculating initial G0_omega"<<std::endl;
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
    FourierTransformer::generate_transformer(parms, fourier_ptr, epsav(), epssqav());
    G0_omega_new.write("G0_omega_input");
    transform_into_real_space(G0_omega_new).write("G0_omega_input_realspace");
    itime_green_function_t G0_tau(n_tau, n_cluster_sites(), n_flavor);
    fourier_ptr->backward_ft(G0_tau, G0_omega_new);
    G0_tau.write("G0_tau_input");
    transform_into_real_space(G0_tau).write("G0_tau_input_realspace");
    return G0_omega_new;
  }
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
