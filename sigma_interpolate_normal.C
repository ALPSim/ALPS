#include<iostream>
#include<fstream>
#include<stdexcept>
#include<vector>
#include<complex>
#include<cmath>
#include<cstdlib>
#include"nambu_green_function.h"
#include<alps/parameter.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "cluster.h"
#include <boost/program_options.hpp>


int main(int argc, char **argv){
  namespace po = boost::program_options;
  std::string pfilename;
  std::string latticefilename;
  std::string clustername;
  int new_cluster_sites=-1;
  // Declare the supported options.
  double mult_anomalous;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "show this help")
    ("beta", po::value<double>(), "new inverse temperature")
    ("mu", po::value<double>(), "new chemical pot")
    ("n_matsubara", po::value<int>(), "new number of matsubara frequencies")
    ("pfile", po::value<std::string>(&pfilename)->default_value("nambu.param"), "parameter file")
    ("lattice", po::value<std::string>(&latticefilename)->default_value("cluster.xml"), "lattice file cluster.xml")
    ("cluster", po::value<std::string>(&clustername), "cluster geometry")
    ("sites", po::value<int>(&new_cluster_sites), "number of cluster sites for new cluster")
    ("iteration", po::value<int>(), "iteration")
          ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    std::cout<<desc;
    return 1;
  }
  if(!vm.count("pfile")) throw std::runtime_error("you need to specify the parameter file as pfile.");
  std::ifstream parm_file(pfilename.c_str());
  if(!parm_file.is_open()) throw std::runtime_error("parm file not open.");
  alps::Parameters parms(parm_file);
  double beta=parms["BETA"];
  double mu=parms["MU"];
  double h=parms["H"];
  double nsites_old=parms["SITES"];
  double nsites_new=vm.count("sites")?new_cluster_sites:(int)(parms["SITES"]);
  int nfreq = parms["NMATSUBARA"];
  int ntau= parms["N"];
  double beta_new,mu_new,eta1_new;
  int nfreq_new;
  if(vm.count("beta")) beta_new=vm["beta"].as<double>(); else beta_new=beta;
  if(vm.count("mu")) mu_new=vm["mu"].as<double>(); else mu_new=mu;
  if(vm.count("n_matsubara")) nfreq_new=vm["n_matsubara"].as<int>(); else nfreq_new=nfreq;
  bool paramagnet=true;
  
  std::cout<<"adjusting parameters: "<<std::endl;
  std::cout<<"new beta: "<<beta_new<<"(old: "<<beta<<")"<<std::endl;
  std::cout<<"new mu: "<<mu_new<<"(old: "<<mu<<")"<<std::endl;
 
  if(!vm.count("lattice")) throw std::runtime_error("you need to specify the cluster.xml file as lattice.");
  parms["LATTICE_LIBRARY"]=vm["lattice"].as<std::string>();
  
  matsubara_green_function_t selfenergy_orig(nfreq, nsites_old, 2, diagonal);
  matsubara_green_function_t selfenergy_kspace_interp(nfreq, nsites_new, 2, diagonal);
  matsubara_green_function_t selfenergy_new(nfreq_new, nsites_new, 2, diagonal);
  matsubara_green_function_t green(nfreq, nsites_new, 2, diagonal);
  matsubara_green_function_t green0(nfreq_new, nsites_new, 2, diagonal);
  itime_green_function_t green0_itime(ntau, nsites_new, 2, diagonal);
  
  if(!vm.count("iteration")) throw std::runtime_error("please specify iteration");
  int iteration=vm["iteration"].as<int>();
  std::stringstream selfenergy_name; selfenergy_name<<"selfenergy_"<<iteration;
  {
    std::ifstream sigma_file(selfenergy_name.str().c_str());
    if(!sigma_file.is_open()) throw std::runtime_error("self energy file not open.");
    double ignored;
    for(int i=0;i<nfreq;++i){
      sigma_file>>ignored; //omega_0 = (2*0+1)*pi/beta, beta=pi/omega_0
      for(int j=0;j<nsites_old;++j){
        double re00,im00,re01,im01,re10,im10,re11,im11;
        sigma_file>>re00>>im00>>re11>>im11;
        selfenergy_orig(i,j,j,0)=std::complex<double>(re00,im00);
        selfenergy_orig(i,j,j,1)=std::complex<double>(re11,im11);
      }
    }
  }
  if(vm.count("cluster")){
    boost::shared_ptr<DCATransformer> clusterhandler_old;
    boost::shared_ptr<DCATransformer> clusterhandler_new;
    
    DCATransformer::generate_transformer(parms, clusterhandler_old);
    if(new_cluster_sites==-1) throw std::logic_error("please specify new cluster sites!");
    parms["SITES"]=new_cluster_sites;
    parms["LATTICE"]=clustername;
    DCATransformer::generate_transformer(parms, clusterhandler_new);
    std::vector<std::complex<double> > values(nsites_old);
    std::vector<double> kvec(2);
    for(int i=0;i<nfreq;++i){
      for(int flavor=0;flavor<2;++flavor){
          for(int k=0;k<nsites_old;++k){
            values[k]=selfenergy_orig(i,k,k,flavor);
          }
          for(int k=0;k<nsites_new;++k){
            for (unsigned int d=0; d<2; ++d){
              kvec[d] = clusterhandler_new->cluster_momentum(k,d);
            }
            selfenergy_kspace_interp(i,k,k,flavor)=clusterhandler_old->interpolate(kvec,values, false);
          }
        }
    }
  }else{
    selfenergy_kspace_interp=selfenergy_orig;
  }
  
  std::vector<double> matsubara_freqs_old(nfreq);
  std::vector<double> matsubara_freqs_new(nfreq);
  for(int i=0;i<nfreq;++i){
    matsubara_freqs_old[i]=(2.*i+1)*M_PI/beta;
    matsubara_freqs_new[i]=(2.*i+1)*M_PI/beta_new;
  }
  //interpolate the self energy
  gsl_interp_accel *acc_real = gsl_interp_accel_alloc ();
  gsl_spline *spline_real = gsl_spline_alloc (gsl_interp_cspline, nfreq);
  gsl_interp_accel *acc_imag = gsl_interp_accel_alloc ();
  gsl_spline *spline_imag = gsl_spline_alloc (gsl_interp_cspline, nfreq);
  
  std::vector<double> sigma_real(nfreq);
  std::vector<double> sigma_imag(nfreq);
  for(int j=0;j<nsites_new;++j){
    for(int flavor=0;flavor<2;++flavor){
        for(int i=0;i<nfreq;++i){
          sigma_real[i]=selfenergy_kspace_interp(i,j,j,flavor).real();
          sigma_imag[i]=selfenergy_kspace_interp(i,j,j,flavor).imag();
        }
        gsl_spline_init (spline_real, &(matsubara_freqs_old[0]), &(sigma_real[0]), nfreq);
        gsl_spline_init (spline_imag, &(matsubara_freqs_old[0]), &(sigma_imag[0]), nfreq);
        
        double xi,yi;                     
        for (int i=0; i<nfreq_new; ++i){
          if(i<nfreq){
            selfenergy_new(i,j,j,flavor) = std::complex<double>(gsl_spline_eval (spline_real, matsubara_freqs_new[i], acc_real), gsl_spline_eval (spline_imag, matsubara_freqs_new[i], acc_imag));
          }else{
            selfenergy_new(i,j,j,flavor) = selfenergy_kspace_interp(nfreq-1,j,j,flavor); //TODO: append tail from high freq expansion.
          }
        }
    }
  }
  gsl_spline_free (spline_real);
  gsl_spline_free (spline_imag);
  gsl_interp_accel_free (acc_real);
  gsl_interp_accel_free (acc_imag);
  
  //compute the self consistency
  //std::cout<<"new self energy is: "<<std::make_pair(selfenergy_new,beta_new)<<std::endl;
  parms["BETA"]=beta_new;
  parms["MU"]=mu_new;

  SimpleDCATransformer clustertrans(parms);
  green0 = clustertrans(green, green0, selfenergy_new, mu_new, h, beta_new, true);  
  
  green0 = clustertrans.symmetrize(green0, paramagnet);
  green0 = clustertrans.symmetrize_k_space(green0);
  
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  fourier_ptr->frequency_to_time_ft(green0_itime, green0);
  std::ofstream green0_interpolated_tau_file("green0_interpolated_tau.dat");
  
  green0_interpolated_tau_file<<std::make_pair(green0_itime,beta_new)<<std::endl;
  std::stringstream filename_new; filename_new<<"G0_omega_input_beta"<<beta_new<<"_mu"<<mu_new;
  std::stringstream filename_sigma_new; filename_sigma_new<<"sigma_beta"<<beta_new;
  green0.write(filename_new.str().c_str());
  std::ofstream sigma_file(filename_sigma_new.str().c_str());
  sigma_file<<std::make_pair(selfenergy_new,beta)<<std::endl;
}

