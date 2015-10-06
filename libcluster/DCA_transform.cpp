/*****************************************************************************
 *
 * Copyright (C) 2007-2015 by Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#include "alps/numeric/vector_functions.hpp"
#include "alps/hdf5/vector.hpp"

#include "cluster.h"
//#include "cluster_fourier_transform.h"

#include "DCA_transform.hxx"
#include "SimpleDCA_transform.hxx"
#include "CDMFT_transform.h"

bool adjust_chemical_potential(const double &density, double &mu, bool &damping_mu,
    bool &computing_mu, int &mu_it, double &n_old, double &n_target,
    double &mu_old, double &mu_new, double &mu_orig);



matsubara_green_function_t ClusterTransformer::initial_G0(const alps::params& parms, shape_t shape) {
  std::cout<<"reading initial G0_omega from: "<<parms["G0OMEGA_INPUT"]<<std::endl;
  int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  int n_orbital=boost::lexical_cast<int>(parms["FLAVORS"]);
  int n_site=boost::lexical_cast<int>(parms["SITES"]);
// TODO: This needs to be checked.  Is it always ok to assume that the G0 is diagonal?
  matsubara_green_function_t G0_omega(n_matsubara, n_site, n_orbital, shape);
//  matsubara_green_function_t G0_omega(n_matsubara, n_site, n_orbital, diagonal);
  std::string g0path=parms["G0OMEGA_INPUT"];
  G0_omega.read(g0path.c_str());
  return G0_omega;
}

matsubara_green_function_t DCATransformer::initial_G0(const alps::params& parms, shape_t shape) {
  return ClusterTransformer::initial_G0(parms, shape);
}

//Don't forget to delete the transformer afterwards!
ClusterTransformer* ClusterTransformer::generate_transformer(const alps::params& parms)
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


DCATransformer::DCATransformer(const alps::params& params)
: ClusterTransformer(params),
t_(params.defined("W") ? 1. : static_cast<double>(params["t"])),
t2_(params["t'"]|params["t2"]|0.),
t3_(params["t3"]|0.),
t4_(params["t4"]|0.),
t_flavor_(params["t_flavor"]|1.0),
bandwidth_(params["W"]|0.),
bandwidth_scale_factor_(1.),
half_filling_(0.),
epsav_(n_cluster_sites()),
epssqav_(n_cluster_sites()),
dispersion_moments_(2, 0.),
params_(params)
{
  if (params.defined("t") && params.defined("W"))
    boost::throw_exception(std::invalid_argument("both t and W are defined"));
  //calculating dispersion moments, bandwidth and half filling
  std::vector<double> energies;
  for (int kc=0; kc<n_cluster_sites(); ++kc) {
    for (int kl=0; kl<n_lattice_momenta(); ++kl) {
      std::vector<double> k(dimension());
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
  if (params["SHIFT_TO_HALF_FILLING"]|false) {
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
}


std::complex<double>* ClusterTransformer::exp_iksite(const int sign) const
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

int* DCATransformer::momentum_translation() const
{
  int *m = new int[n_cluster_sites()*n_cluster_sites()*n_cluster_sites()];
  std::valarray<double> phase(dimension());
  int res=-1;
  for(int i=0; i<n_cluster_sites();++i){
    for(int j=0;j<n_cluster_sites();++j){
      for(int k=0;k<n_cluster_sites();++k){
        for(std::size_t d=0;d<dimension();++d){
          phase[d]=cluster_momentum(i,d)-cluster_momentum(j,d)+cluster_momentum(k,d);
        }
        for(int p=0;p<n_cluster_sites();++p){
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

void DCATransformer::print_momenta(const char * filename, const matsubara_green_function_t & G) const{
  std::ofstream momenta_file;
  momenta_file.open(filename);
  momenta_file<<std::setprecision(14);

  int dim = this->dimension();
  shape_t gf_shape = G.shape();

  if (gf_shape == diagonal) {
    for (unsigned int k = 0; k < G.nsite(); ++k) {
      for (unsigned int d = 0; d < dim; ++d) {
        momenta_file << cluster_momentum(k, d) << "\t";
      }
      momenta_file << "\n";
    }
  } else if (gf_shape == blockdiagonal) {
    for (unsigned int k = 0; k < G.nsite(); k += 2)
      for (unsigned int p1 = 0; p1 <= 1; p1++)
        for (unsigned int p2 = 0; p2 <= 1; p2++) {
          for (unsigned int d = 0; d < dim; ++d) {
            momenta_file << cluster_momentum(k + p1, d) << "\t";
          }
          momenta_file << "\n";
          for (unsigned int d = 0; d < dim; ++d) {
            momenta_file << cluster_momentum(k + p2, d) << "\t";
          }
          momenta_file << "\n";
        }

  } else if (gf_shape == nondiagonal) {
    for (unsigned int s1 = 0; s1 < G.nsite(); ++s1)
      for (unsigned int s2 = 0; s2 < G.nsite(); ++s2) {
        for (unsigned int d = 0; d < dim; ++d) {
          momenta_file << cluster_momentum(s1, d) << "\t";
        }
        momenta_file << "\n";
        for (unsigned int d = 0; d < dim; ++d) {
          momenta_file << cluster_momentum(s2, d) << "\t";
        }
        momenta_file <<"\n";
      }
  }

  momenta_file.close();
}

std::string DCATransformer::generate_header_kspace(const matsubara_green_function_t & G) const
{
  std::string header="Freq\t\t";
  int dim = this->dimension();
  shape_t gf_shape = G.shape();

  if (gf_shape==diagonal) {
    for(unsigned int k=0; k<G.nsite(); ++k){
      std::ostringstream temp;
      temp.precision(3);
      temp << "(";
      for(unsigned int d=0; d<dim-1; ++d){
        temp<<cluster_momentum(k,d)<<",";
      }
      temp<<cluster_momentum(k,dim-1)<<")";
      for(unsigned int f=0;f<G.nflavor();++f){
        header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
      }
    }
  }else if(gf_shape==blockdiagonal){
    for (unsigned int k=0; k<G.nsite(); k+=2)
      for (unsigned int p1=0; p1<=1; p1++)
        for (unsigned int p2=0; p2<=1; p2++){
          std::ostringstream temp;
          temp.precision(3);
          temp << "(";
          for(unsigned int d=0; d<dim-1; ++d){
            temp<<cluster_momentum(k+p1,d)<<",";
          }
          temp<<cluster_momentum(k+p1,dim-1)<<")";
          temp << "(";
          for(unsigned int d=0; d<dim-1; ++d){
            temp<<cluster_momentum(k+p2,d)<<",";
          }
          temp<<cluster_momentum(k+p2,dim-1)<<")";
          for(unsigned int f=0;f<G.nflavor();++f){
            header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
          }
        }

  }else if(gf_shape==nondiagonal){
    for (unsigned int s1=0; s1<G.nsite(); ++s1)
      for (unsigned int s2=0; s2<G.nsite(); ++s2){

        std::ostringstream temp;
        temp.precision(3);
        temp << "(";
        for(unsigned int d=0; d<dim-1; ++d){
          temp<<cluster_momentum(s1,d)<<",";
        }
        temp<<cluster_momentum(s1,dim-1)<<")";
        temp << "(";
        for(unsigned int d=0; d<dim-1; ++d){
          temp<<cluster_momentum(s2,d)<<",";
        }
        temp<<cluster_momentum(s2,dim-1)<<")";
        for(unsigned int f=0;f<G.nflavor();++f){
          header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
        }
      }
  }

  return header;
}

std::string DCATransformer::generate_header_realspace(const matsubara_green_function_t & G) const
{
  std::string header="Freq\t\t";
  int dim = this->dimension();
  shape_t gf_shape = G.shape();

  if (gf_shape==diagonal) {
    for(unsigned int k=0; k<G.nsite(); ++k){
      std::ostringstream temp;
      temp.precision(3);
      temp << "(";
      for(unsigned int d=0; d<dim-1; ++d){
        temp<<cluster_coordinate(k,d)<<",";
      }
      temp<<cluster_coordinate(k,dim-1)<<")";
      for(unsigned int f=0;f<G.nflavor();++f){
        header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
      }
    }
  }else if(gf_shape==blockdiagonal){
    for (unsigned int k=0; k<G.nsite(); k+=2)
      for (unsigned int p1=0; p1<=1; p1++)
        for (unsigned int p2=0; p2<=1; p2++){
          std::ostringstream temp;
          temp.precision(3);
          temp << "(";
          for(unsigned int d=0; d<dim-1; ++d){
            temp<<cluster_coordinate(k+p1,d)<<",";
          }
          temp<<cluster_coordinate(k+p1,dim-1)<<")";
          temp << "(";
          for(unsigned int d=0; d<dim-1; ++d){
            temp<<cluster_coordinate(k+p2,d)<<",";
          }
          temp<<cluster_coordinate(k+p2,dim-1)<<")";
          for(unsigned int f=0;f<G.nflavor();++f){
            header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
          }
        }

  }else if(gf_shape==nondiagonal){
    for (unsigned int s1=0; s1<G.nsite(); ++s1)
      for (unsigned int s2=0; s2<G.nsite(); ++s2){

        std::ostringstream temp;
        temp.precision(3);
        temp << "(";
        for(unsigned int d=0; d<dim-1; ++d){
          temp<<cluster_coordinate(s1,d)<<",";
        }
        temp<<cluster_coordinate(s1,dim-1)<<")";
        temp << "(";
        for(unsigned int d=0; d<dim-1; ++d){
          temp<<cluster_coordinate(s2,d)<<",";
        }
        temp<<cluster_coordinate(s2,dim-1)<<")";
        for(unsigned int f=0;f<G.nflavor();++f){
          header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
        }
      }
  }

  return header;
}

std::string DCATransformer::generate_header_kspace(const itime_green_function_t & G) const
{
  std::string header="Tau\t";
  int dim = this->dimension();
  shape_t gf_shape = G.shape();

  if (gf_shape==diagonal) {
    for(unsigned int k=0; k<G.nsite(); ++k){
      std::ostringstream temp;
      temp.precision(3);
      temp << "(";
      for(unsigned int d=0; d<dim-1; ++d){
        temp<<cluster_momentum(k,d)<<",";
      }
      temp<<cluster_momentum(k,dim-1)<<")";
      for(unsigned int f=0;f<G.nflavor();++f){
        header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
      }
    }
  }else if(gf_shape==blockdiagonal){
    for (unsigned int k=0; k<G.nsite(); k+=2)
      for (unsigned int p1=0; p1<=1; p1++)
        for (unsigned int p2=0; p2<=1; p2++){
          std::ostringstream temp;
          temp.precision(3);
          temp << "(";
          for(unsigned int d=0; d<dim-1; ++d){
            temp<<cluster_momentum(k+p1,d)<<",";
          }
          temp<<cluster_momentum(k+p1,dim-1)<<")";
          temp << "(";
          for(unsigned int d=0; d<dim-1; ++d){
            temp<<cluster_momentum(k+p2,d)<<",";
          }
          temp<<cluster_momentum(k+p2,dim-1)<<")";
          for(unsigned int f=0;f<G.nflavor();++f){
            header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
          }
        }

  }else if(gf_shape==nondiagonal){
    for (unsigned int s1=0; s1<G.nsite(); ++s1)
      for (unsigned int s2=0; s2<G.nsite(); ++s2){
        std::ostringstream temp;
        temp.precision(3);
        temp << "(";
        for(unsigned int d=0; d<dim-1; ++d){
          temp<<cluster_momentum(s1,d)<<",";
        }
        temp<<cluster_momentum(s1,dim-1)<<")";
        temp << "(";
        for(unsigned int d=0; d<dim-1; ++d){
          temp<<cluster_momentum(s2,d)<<",";
        }
        temp<<cluster_momentum(s2,dim-1)<<")";
        for(unsigned int f=0;f<G.nflavor();++f){
          header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
        }
      }
  }

  return header;
}

std::string DCATransformer::generate_header_realspace(const itime_green_function_t & G) const{
  std::string header="Tau\t";
  int dim = this->dimension();
  shape_t gf_shape = G.shape();

  if (gf_shape==diagonal) {
    for(unsigned int k=0; k<G.nsite(); ++k){
      std::ostringstream temp;
      temp.precision(3);
      temp << "(";
      for(unsigned int d=0; d<dim-1; ++d){
        temp<<cluster_coordinate(k,d)<<",";
      }
      temp<<cluster_coordinate(k,dim-1)<<")";
      for(unsigned int f=0;f<G.nflavor();++f){
        header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
      }
    }
  }else if(gf_shape==blockdiagonal){
    for (unsigned int k=0; k<G.nsite(); k+=2)
      for (unsigned int p1=0; p1<=1; p1++)
        for (unsigned int p2=0; p2<=1; p2++){
          std::ostringstream temp;
          temp.precision(3);
          temp << "(";
          for(unsigned int d=0; d<dim-1; ++d){
            temp<<cluster_coordinate(k+p1,d)<<",";
          }
          temp<<cluster_coordinate(k+p1,dim-1)<<")";
          temp << "(";
          for(unsigned int d=0; d<dim-1; ++d){
            temp<<cluster_coordinate(k+p2,d)<<",";
          }
          temp<<cluster_coordinate(k+p2,dim-1)<<")";
          for(unsigned int f=0;f<G.nflavor();++f){
            header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
          }
        }

  }else if(gf_shape==nondiagonal){
    for (unsigned int s1=0; s1<G.nsite(); ++s1)
      for (unsigned int s2=0; s2<G.nsite(); ++s2){
        std::ostringstream temp;
        temp.precision(3);
        temp << "(";
        for(unsigned int d=0; d<dim-1; ++d){
          temp<<cluster_coordinate(s1,d)<<",";
        }
        temp<<cluster_coordinate(s1,dim-1)<<")";
        temp << "(";
        for(unsigned int d=0; d<dim-1; ++d){
          temp<<cluster_coordinate(s2,d)<<",";
        }
        temp<<cluster_coordinate(s2,dim-1)<<")";
        for(unsigned int f=0;f<G.nflavor();++f){
          header.append(temp.str()+"f"+boost::lexical_cast<std::string>(f)+"\t");
        }
      }
  }

  return header;
}

matsubara_green_function_t SimpleDCATransformer::operator()(const matsubara_green_function_t &G_omega,
    const matsubara_green_function_t &G0_omega,
    double &mu, const double h, const double beta)
{
  matsubara_green_function_t selfenergy(G_omega);
  return this->operator()(G_omega, G0_omega, selfenergy, mu, h, beta, false);
}

void SimpleDCATransformer::compute_sigma_from_G0_and_G(const matsubara_green_function_t& G0_omega,
    const matsubara_green_function_t& G_omega,
    matsubara_green_function_t& selfenergy) const {
  //compute selfenergy out of G and G0 first
  unsigned int nflavor=G0_omega.nflavor();
  unsigned int nfreq=G0_omega.nfreq();

  for (int flavor = 0; flavor < nflavor; ++flavor) {
    for (int kc = 0; kc < n_cluster_sites(); ++kc) {
      for (int freq = 0; freq < nfreq; ++freq) {
        selfenergy(freq, kc, kc, flavor) = 1. / G0_omega(freq, kc, kc, flavor);
        selfenergy(freq, kc, kc, flavor) -= 1. / G_omega(freq, kc, kc, flavor);
      }
    }
  }
}
void SimpleDCATransformer::compute_integration_weight_and_dispersion(bool use_phase,std::vector<std::vector<std::vector<double> > > &epsilon_kckl, std::vector<double> &weight_kl, unsigned int nflavor )const{
  std::vector<double> phi=use_phase?phase_vector():std::vector<double>();
  for (int kc=0; kc<n_cluster_sites(); ++kc) {
    for (int kl=0; kl<n_lattice_momenta(); ++kl) {
      std::vector<double> k(dimension());
      for (unsigned int d=0; d<dimension(); ++d)
        k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
for (unsigned int flavor=0; flavor<nflavor; ++flavor){
        epsilon_kckl[kc][kl][flavor]=use_phase?epsilon_phi(k,phi):epsilon_flavor(k,flavor);
	}
     }
  }


  for (int kl=0; kl<n_lattice_momenta(); ++kl) {
    weight_kl[kl]= lattice_momentum(kl,dimension());
  }
}
void SimpleDCATransformer::compute_lattice_green_function(matsubara_green_function_t &G_omega_lattice, const matsubara_green_function_t &selfenergy, const double &mu, const double &beta, const double &h, std::vector<std::vector<std::vector<double> > > &epsilon_kckl, std::vector<double> &weight_kl) const{
  unsigned int nflavor=selfenergy.nflavor();
  unsigned int nfreq=selfenergy.nfreq();
  for (int flavor=0; flavor<nflavor; ++flavor) {
    for (int kc=0; kc<n_cluster_sites(); ++kc) {
      for (int freq=0; freq<nfreq; ++freq) {
        G_omega_lattice(freq, kc, kc, flavor) = 0.;
        std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
        for (int kl=0; kl<n_lattice_momenta(); ++kl) {
          std::complex<double> zeta = iomega_n + mu + (flavor%2 ? -h : h);
          G_omega_lattice(freq, kc, kc, flavor) += weight_kl[kl]/(zeta-epsilon_kckl[kc][kl][flavor]-selfenergy(freq, kc, kc, flavor));
        }
        G_omega_lattice(freq, kc, kc, flavor) /= n_cluster_cells();
      }
    }
  }
}
///compute the density, given a Green's function
///TODO: we will need to make this more accurate by adding the high frequency behavior; at least the c2 is needed here.
double SimpleDCATransformer::compute_density(const matsubara_green_function_t &G_omega_lattice, const double &beta){
  unsigned int nflavor=G_omega_lattice.nflavor();
  unsigned int nfreq=G_omega_lattice.nfreq();
  double density=0;
  for(int flavor=0;flavor<nflavor;++flavor){
    for(int kc=0;kc<n_cluster_sites();++kc){
      for(int freq=0;freq<nfreq;++freq){
        //we need to make this more accurate by using hi freq moments
        density+=2.*G_omega_lattice(freq,kc,kc,flavor).real();
      }
    }
  }
  density=density/beta/n_cluster_sites()+(0.5*n_cluster_sites())/n_cluster_sites()*nflavor; //density such that 'half filling' corresponds to 1
  return density;
}

void SimpleDCATransformer::update_parameter_file(int current_it, double& mu) {
  current_it = current_it + 1;
  (*const_cast<alps::params*>(&params_))["CURRENT_IT"] = current_it;
  (*const_cast<alps::params*>(&params_))["MU"] = mu;
  std::cout << "Updating mu in basename file with value " << mu << std::endl;
  std::string basenamepath = params_["BASENAME"];
  std::ofstream os(basenamepath.c_str());
  os << params_;
  os.close();
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
    compute_sigma_from_G0_and_G(G0_omega, G_omega, selfenergy);
  }
  bool converged=compute_filling_?false:true;
  bool damping_mu=false;
  bool computing_mu=true;
  int mu_it=1;


  int current_it = params_["CURRENT_IT"]|0;
  bool use_phase = params_["PHASE"]|false;

 // std::vector<std::vector<double> >epsilon_kckl(n_cluster_sites(), std::vector<double>(n_lattice_momenta()));

 std::vector<std::vector<std::vector<double> > > epsilon_kckl(n_cluster_sites(), std::vector<std::vector<double> >(n_lattice_momenta(), std::vector<double>(nflavor)) );
  std::vector<double> weight_kl(n_lattice_momenta());

  compute_integration_weight_and_dispersion(use_phase, epsilon_kckl, weight_kl, nflavor);

  if(!compute_filling_ || current_it < 2)
    compute_lattice_green_function(G_omega_lattice, selfenergy, mu, beta, h, epsilon_kckl, weight_kl);
  else{
    //find the right particle density
    do{

      //--> this will need to be reworked. Here the density computation should be inside the mu calculation.

      double n_old;
      double mu_old=mu;
      double mu_new=mu;
      double mu_orig=mu;

      //compute density
      double density=compute_density(G_omega_lattice, beta);

      std::cout<<"mu: "<<mu<<" filling "<<density<<" it: "<<mu_it<<" converged: "<<converged<<std::endl;
        //compute new chemical potential
        converged=adjust_chemical_potential(density, mu, damping_mu, computing_mu, mu_it, n_old, n_target_, mu_old, mu_new, mu_orig);

        std::cout<<"converged now: "<<converged<<" true is: "<<true<<" compute filling is: "<<compute_filling_<<std::endl;
    } while(!converged);
    update_parameter_file(current_it, mu);
  }
  for (int flavor=0; flavor<nflavor; ++flavor) {
    for (int kc=0; kc<n_cluster_sites(); ++kc) {
      for (int freq=0; freq<nfreq; ++freq) {
        G0_omega_new(freq, kc, kc, flavor) =
            1. / ( 1./G_omega_lattice(freq, kc, kc, flavor) + selfenergy(freq, kc, kc, flavor) );
      }
    }
  }


  return G0_omega_new;
}



matsubara_green_function_t SimpleDCATransformer::initial_G0(const alps::params& parms, shape_t shape)
{
  if (parms.defined("G0OMEGA_INPUT")){
    std::string g0path=parms["G0OMEGA_INPUT"];
    std::ifstream in_file(g0path.c_str());
    if(in_file.good()){
      return DCATransformer::initial_G0(parms, shape);
    }
    else
      std::cout<<"WARNING: Input file could not be found, proceeding with free solution."<<std::endl;
  }
  //proceed with computation of new G0 if we can't find input G0
  std::cout<<"calculating initial G0_omega"<<std::endl;
  unsigned int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  unsigned int n_tau=boost::lexical_cast<int>(parms["N"]);
  unsigned int n_flavor=parms["FLAVORS"]|2;
  double beta = static_cast<double>(parms["BETA"]);
  double mu = static_cast<double>(parms["MU"]);
  //double h = static_cast<double>(parms["H"]);
  //bool paramagnet = static_cast<bool>(parms["PARAMAGNET"]);
  matsubara_green_function_t G0_omega_new(n_matsubara,n_cluster_sites(),n_flavor, diagonal);
  G0_omega_new.clear();
  if (parms["INSULATING"]|false) {
    std::cout <<"calculating insulating G0_omega";
    for (int freq=0; freq<n_matsubara; ++freq) {
      const std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      const std::complex<double> G0 = 1./(iomega_n+mu-dispersion_moments(0));
      for (int flavor=0; flavor<n_flavor; ++flavor) {
        for (int kc=0; kc<n_cluster_sites(); ++kc) {
          G0_omega_new(freq, kc, kc, flavor) = G0;
        }
      }
    }
  }
  else if (parms["ISOLATED"]|false) {
    std::cout << "calculating G0_omega for isolated cluster"<<std::endl;
    bool phase=parms["PHASE"]|false;
    for (int freq=0; freq<n_matsubara; ++freq) {
      for (int flavor=0; flavor<n_flavor; ++flavor) {
        const std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
        for (int kc=0; kc<n_cluster_sites(); ++kc) {
          std::vector<double> k(dimension());
          std::vector<double> phi(dimension(),0.);
          if(phase) phi=phase_vector();
          for (unsigned int d=0; d<dimension(); ++d)
            k[d] = cluster_momentum(kc,d);
          const std::complex<double> G0 = 1./(iomega_n+mu-epsilon_phi(k,phi));
          G0_omega_new(freq, kc, kc, flavor) = G0;
        }
      }
    }
  }
  else {
    matsubara_green_function_t G0_omega(n_matsubara,n_cluster_sites(),n_flavor,  diagonal);
    matsubara_green_function_t G_omega(n_matsubara,n_cluster_sites(),n_flavor,  diagonal);
    G0_omega.clear();
    for (int flavor=0; flavor<n_flavor; ++flavor) {
      for (int i=0; i<n_cluster_sites(); ++i) {
        for (int freq=0; freq<n_matsubara; ++freq) {
          G0_omega(freq, i, i, flavor) = 1.;
        }
      }
    }
    G_omega = G0_omega;
    G0_omega_new = operator()(G_omega, G0_omega, mu, 0, beta);
  }
  itime_green_function_t G0_tau(n_tau+1, n_cluster_sites(), n_flavor,  diagonal);
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  G0_omega_new.write("G0_omega_input");
  FourierTransformer::generate_cluster_transformer(parms, fourier_ptr);
  fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega_new);
  G0_tau.write("G0_tau_input");
  if(!parms.defined("WRITE_SPARSE")){
    matsubara_green_function_t G0_omega_real = transform_into_real_space(G0_omega_new);
    G0_omega_real.write("G0_omega_input_realspace");
    itime_green_function_t G0_tau_real = transform_into_real_space(G0_tau);
    G0_tau_real.write("G0_tau_input_realspace");
  }
  return G0_omega_new;
}
//fourier transform epsilon_bar to real space using the same convention as in the real-space FT of the GF
std::vector<double> DCATransformer::tmatrix() const{
  unsigned int nsite=n_cluster_sites();
  std::vector<double> tij(nsite*nsite, 0.);
  for (int i=0; i<nsite; ++i) {
    for (int j=0; j<nsite; ++j) {
      for (int k=0; k<nsite; ++k) {
        double phase = 0.;
        for (unsigned int d=0; d<dimension(); ++d)
          phase += cluster_momentum(k,d)*(cluster_coordinate(i,d)-cluster_coordinate(j,d));
        //std::cout<<"i j k: "<<i<<" "<<j<<" "<<k<<" "<<" epsav: "<<epsav_[k]<<" phase: "<<phase<<std::endl;
        tij[i*nsite+j]+=epsav_[k]*cos(phase)/double(nsite);
      }
    }
  }
  return tij;

}
//write all the essential information to file: averaged moments etc.
void DCATransformer::write_hdf5(alps::hdf5::archive &ar, const std::string &p) const{
  std::cout<<"writing DCA transformer information to: "<<p<<" in file: "<<ar.get_filename()<<std::endl;
  std::cout<<"using row major format"<<std::endl;
  std::cout<<" name is: "<<(p+"/epsav")<<std::endl;
  ar<<alps::make_pvp(p+"/epsav", epsav_);
  ar<<alps::make_pvp(p+"/epssqav", epssqav_);
  ar<<alps::make_pvp(p+"/dispersion_moments", dispersion_moments_);
  ar<<alps::make_pvp(p+"/tmatrix", tmatrix());
  ClusterTransformer::write_hdf5(ar,p);
}
void ClusterTransformer::write_hdf5(alps::hdf5::archive &ar, const std::string &p) const{
  ClusterHelper::write_hdf5(ar,p);
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
  double density_tolerance=1.e-10;
  int max_mu_it =10;

  std::cout<<"Adjusting chemical potential: current chem. pot: "<<mu<<" previous chem. pot.: "<<mu_old<<" mu orig: "
      <<mu_orig<<" current density: "<<n_tot<<" previous density: "<<n_old<<" target density: "<<n_target<<std::endl;
  if(!damping_mu && (std::abs(n_tot-n_target)<density_tolerance || mu_it==10)){
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
      //regula falsi
      //use c_k = b_k -f(b_k)*(b_k - a_k)/(f(b_k) - f(a_k))
      //mu_new = mu_old -n_old*(mu_old - mu)/(n_tot-n_old);
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

void SimpleDCATransformer::construct_green_function(double mu, double beta, double eta1, double eta2, const std::vector<std::vector<std::vector<double> > > &epsilon_kckl, const std::vector<std::vector<double> > &dwave_prefactor_kckl, const std::vector<double> &weight_kl, const nambu_matsubara_green_function_t &selfenergy,nambu_matsubara_green_function_t &G_omega_lattice_kav){
  
  unsigned int nflavor=G_omega_lattice_kav.nflavor();
  unsigned int nfreq=G_omega_lattice_kav.nfreq();
  for (int kc=0; kc<n_cluster_sites(); ++kc) {
    for (int freq=0; freq<nfreq; ++freq) {
      G_omega_lattice_kav(freq, kc, kc, 0,0) = 0.;
      G_omega_lattice_kav(freq, kc, kc, 0,1) = 0.;
      G_omega_lattice_kav(freq, kc, kc, 1,0) = 0.;
      G_omega_lattice_kav(freq, kc, kc, 1,1) = 0.;
      std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      for (int kl=0; kl<n_lattice_momenta(); ++kl) {
        //double prefactor=1.;//(k[0]==0 && k[1]==0)?1.:0.; //cos kx - cos ky
        double prefactor=dwave_prefactor_kckl[kc][kl];
        //this is the self consistency formula!!!
        //note the pairing field eta + i eta2 on the 0,1 component, eta - i eta2 on the 1,0 component.
        std::complex<double> G_omega_latt_inv[4]=
        {iomega_n + mu -epsilon_kckl[kc][kl][0]-selfenergy(freq, kc, kc, 0,0),
          -selfenergy(freq, kc, kc, 0,1) - prefactor*std::complex<double>( eta1, eta2),
          -selfenergy(freq, kc, kc, 1,0) - prefactor*std::complex<double>( eta1,-eta2),
          iomega_n - mu +epsilon_kckl[kc][kl][1]-selfenergy(freq, kc, kc, 1,1)};
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
  for(int kc=0;kc<n_cluster_sites();++kc){
    for(int freq=0;freq<G_omega_lattice_kav.nfreq();++freq){
      //we could make this more accurate by using hi freq moments
      density_up_part+=2./beta*G_omega_lattice_kav(freq,kc,kc,0,0).real(); //particles
      density_dn_hole-=2./beta*G_omega_lattice_kav(freq,kc,kc,1,1).real(); //holes
    }
  }
  //compute density
  density_up_part=(density_up_part+n_cluster_sites()/2.0)/n_cluster_sites(); //density such that 'half filling' corresponds to 1
  density_dn_hole=(density_dn_hole+n_cluster_sites()/2.0)/n_cluster_sites(); //density such that 'half filling' corresponds to 1
  density=density_dn_hole+density_up_part;

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

  //std::vector<std::vector<double> >epsilon_kckl        (n_cluster_sites(), std::vector<double>(n_lattice_momenta()));

  std::vector<std::vector<std::vector<double> > > epsilon_kckl(n_cluster_sites(), std::vector<std::vector<double> >(n_lattice_momenta(), std::vector<double>(nflavor)) );
  std::vector<std::vector<double> >dwave_prefactor_kckl(n_cluster_sites(), std::vector<double>(n_lattice_momenta()));
  std::vector<double> weight_kl(n_lattice_momenta());
  for (int kc=0; kc<n_cluster_sites(); ++kc) {
    for (int kl=0; kl<n_lattice_momenta(); ++kl) {
      std::vector<double> k(dimension());
      for (unsigned int d=0; d<dimension(); ++d)
        k[d] = lattice_momentum(kl,d) + cluster_momentum(kc,d);
     for (unsigned int flavor=0; flavor<nflavor; ++flavor){
        epsilon_kckl[kc][kl][flavor]=epsilon_flavor(k,flavor);
	}
      dwave_prefactor_kckl[kc][kl]= cos(k[0]) - cos(k[1]); //cos kx - cos ky
    }
  }for (int kl=0; kl<n_lattice_momenta(); ++kl) {
    weight_kl[kl]= lattice_momentum(kl,dimension());
  }

  if(!use_selfenergy){ //compute selfenergy out of G and G0 first
    selfenergy=G0_omega;
    //compute the lattice = cluster self energy
    for (int kc=0; kc<n_cluster_sites(); ++kc) {
      for (int freq=0; freq<nfreq; ++freq) {
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
    throw std::runtime_error("fix writing back of chemical potential if needed");
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
    (*const_cast<alps::params*>(&params_))["MU"]=mu;
    //write the chemical potential back to the parameter file
    /// \todo fix this if needed
  }else{
    construct_green_function( mu,  beta,  eta1,  eta2, epsilon_kckl, dwave_prefactor_kckl, weight_kl, selfenergy,G_omega_lattice_kav);
  }




  //compute the bath Green's function G0
  for (int kc=0; kc<n_cluster_sites(); ++kc) {
    for (int freq=0; freq<nfreq; ++freq) {
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

nambu_matsubara_green_function_t ClusterTransformer::initial_G0_nambu(const alps::params& parms) {
  std::cout<<"reading initial (NAMBU) G0_omega from: "<<parms["G0OMEGA_INPUT"]<<std::endl;
  int n_matsubara=boost::lexical_cast<int>(parms["NMATSUBARA"]);
  int n_orbital=boost::lexical_cast<int>(parms["FLAVORS"]);
  int n_site=boost::lexical_cast<int>(parms["SITES"]);
  double beta =boost::lexical_cast<double>(parms["BETA"]);
  if(n_orbital!=2) throw std::invalid_argument("dunno how to handle nambu in >2 orbitals!");
  nambu_matsubara_green_function_t G0_omega(n_matsubara, n_site, n_orbital, diagonal);
  std::string g0path=parms["G0OMEGA_INPUT"];
  G0_omega.read(g0path.c_str());
  return G0_omega;
}

nambu_matsubara_green_function_t SimpleDCATransformer::initial_G0_nambu(const alps::params& parms)
{
  if (parms.defined("G0OMEGA_INPUT")){
    std::string g0path=parms["G0OMEGA_INPUT"];
    std::ifstream in_file(g0path.c_str());
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
  unsigned int n_flavor=parms["FLAVORS"]|2;
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
    for (int freq=0; freq<n_matsubara; ++freq) {
      const std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);
      const std::complex<double> G0 = 1./(iomega_n+mu-dispersion_moments(0));
      for (int flavor=0; flavor<n_flavor; ++flavor) {
        for (int kc=0; kc<n_cluster_sites(); ++kc) {
          G0_omega_new(freq, kc, kc, 0,0) = G0;
          G0_omega_new(freq, kc, kc, 1,1) = G0;
          G0_omega_new(freq, kc, kc, 0,1) = 0.;
          G0_omega_new(freq, kc, kc, 1,0) = 0.;
        }
      }
    }
  }else if (parms["ISOLATED"]|false){
    for (int kc=0; kc<n_cluster_sites(); ++kc) {
      for (int freq=0; freq<n_matsubara; ++freq) {
        G0_omega_new(freq, kc, kc, 0,0) = 0.;
        G0_omega_new(freq, kc, kc, 0,1) = 0.;
        G0_omega_new(freq, kc, kc, 1,0) = 0.;
        G0_omega_new(freq, kc, kc, 1,1) = 0.;
        std::complex<double> iomega_n(0, (2.*freq+1.)*M_PI/beta);

        std::vector<double> k(dimension());
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
    for (int i=0; i<n_cluster_sites(); ++i) {
      for (int freq=0; freq<n_matsubara; ++freq) {
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
  FourierTransformer::generate_cluster_transformer(parms, fourier_ptr);
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
    for (int k1=0; k1<n_cluster_sites(); ++k1) {
      for (int k2=0; k2<n_cluster_sites(); ++k2) {
        for(int n=0; n<G_new.nfreq(); ++n) {
          int other_n=G_new.nfreq()-1-n;
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
    for (int k1=0; k1<n_cluster_sites(); ++k1) {
      for (int k2=0; k2<n_cluster_sites(); ++k2) {
        for(int freq=0; freq<G_new.nfreq(); ++freq) {
          //CLANG really hates individual assignements to real() and imag()
          G_new(freq, k1, k2, 0,0)=std::complex<double>(0.5*(G(freq, k1, k2, 0,0).real()-G_new(freq, k1, k2, 1,1).real()),0.5*(G(freq, k1, k2, 0,0).imag()+G_new(freq, k1, k2, 1,1).imag()));
          G_new(freq, k1, k2, 1,1)=std::complex<double>(-G_new(freq, k1, k2, 0,0).real(),G_new(freq, k1, k2, 0,0).imag());
          //G_new(freq, k1, k2, 0,0).imag()= 0.5*(G(freq, k1, k2, 0,0).imag()+G_new(freq, k1, k2, 1,1).imag());
          //G_new(freq, k1, k2, 0,0).real()= 0.5*(G(freq, k1, k2, 0,0).real()-G_new(freq, k1, k2, 1,1).real());
          //G_new(freq, k1, k2, 1,1).real()= -G_new(freq, k1, k2, 0,0).real();
          //G_new(freq, k1, k2, 1,1).imag()=  G_new(freq, k1, k2, 0,0).imag();
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

double DCATransformer::epsilon_phi(const std::vector<double> &k, const std::vector<double> &phi) const{
  using namespace alps::numeric;
  return epsilon(k+phi);
}
