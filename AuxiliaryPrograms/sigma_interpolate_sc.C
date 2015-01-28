#include<iostream>
#include<fstream>
#include<stdexcept>
#include<vector>
#include<complex>
#include<cmath>
#include<cstdlib>
#include"nambu_green_function.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "cluster.h"
#include "fouriertransform.h"
#include <boost/program_options.hpp>
#include "../clusterhelper.hxx"
#include "../DCA_transformer.hxx"


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
    ("eta", po::value<double>(), "new dwave field")
    ("n_matsubara", po::value<int>(), "new number of matsubara frequencies")
    ("pfile", po::value<std::string>(&pfilename)->default_value("nambu.param"), "parameter file")
    ("lattice", po::value<std::string>(&latticefilename)->default_value("cluster.xml"), "lattice file cluster.xml")
    ("cluster", po::value<std::string>(&clustername), "cluster geometry")
    ("sites", po::value<int>(&new_cluster_sites), "number of cluster sites for new cluster")
    ("iteration", po::value<std::string>(), "iteration, use \"all\" to average all available iterations.")
    ("mult-anomalous", po::value<double>(&mult_anomalous)->default_value(1.), "multiply anomalous self energy by a factor of x")
          ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
    std::cout<<desc;
    return 1;
  }
  if(!vm.count("pfile")) throw std::runtime_error("you need to specify the parameter file as pfile.");
  alps::params parms(pfilename);
  double beta=parms["BETA"];
  double mu=parms["MU"];
  double h=parms["H"];
  double eta1=parms["ETA1"];
  double eta2=parms["ETA2"];
  double nsites_old=parms["SITES"];
  double nsites_new=vm.count("sites")?new_cluster_sites:(int)(parms["SITES"]);
  int nfreq = parms["NMATSUBARA"];
  int ntau= parms["N"];
  double beta_new,mu_new,eta1_new;
  int nfreq_new;
  if(vm.count("beta")) beta_new=vm["beta"].as<double>(); else beta_new=beta;
  if(vm.count("mu")) mu_new=vm["mu"].as<double>(); else mu_new=mu;
  if(vm.count("eta")) eta1_new=vm["eta"].as<double>(); else eta1_new=eta1;
  if(vm.count("n_matsubara")) nfreq_new=vm["n_matsubara"].as<int>(); else nfreq_new=nfreq;
  bool paramagnet=true;
  
  std::cout<<"adjusting parameters: "<<std::endl;
  std::cout<<"new beta: "<<beta_new<<"(old: "<<beta<<")"<<std::endl;
  std::cout<<"new mu: "<<mu_new<<"(old: "<<mu<<")"<<std::endl;
  std::cout<<"new eta: "<<eta1_new<<"(old: "<<eta1<<")"<<std::endl;
 
  if(!vm.count("lattice")) throw std::runtime_error("you need to specify the cluster.xml file as lattice.");
  parms["LATTICE_LIBRARY"]=vm["lattice"].as<std::string>();
  
  nambu_matsubara_green_function_t selfenergy_orig(nfreq, nsites_old, 2, diagonal);
  nambu_matsubara_green_function_t selfenergy_kspace_interp(nfreq, nsites_new, 2, diagonal);
  nambu_matsubara_green_function_t selfenergy_new(nfreq_new, nsites_new, 2, diagonal);
  nambu_matsubara_green_function_t green(nfreq, nsites_new, 2, diagonal);
  nambu_matsubara_green_function_t green0(nfreq_new, nsites_new, 2, diagonal);
  nambu_itime_green_function_t green0_itime(ntau, nsites_new, 2, diagonal);
  
  if(!vm.count("iteration")) throw std::runtime_error("please specify iteration");
  if(vm["iteration"].as<std::string>()!= std::string("all")){
    std::stringstream selfenergy_name; 
    int iteration=vm["iteration"].as<int>();
    selfenergy_name<<"selfenergy_"<<iteration;
    std::ifstream sigma_file(selfenergy_name.str().c_str());
    if(!sigma_file.is_open()) throw std::logic_error("sigma file not open!");
    double ignored;
    for(int i=0;i<nfreq;++i){
      sigma_file>>ignored; //omega_0 = (2*0+1)*pi/beta, beta=pi/omega_0
      for(int j=0;j<nsites_old;++j){
        double re00,im00,re01,im01,re10,im10,re11,im11;
        sigma_file>>re00>>im00>>re01>>im01>>re10>>im10>>re11>>im11;
        selfenergy_orig(i,j,j,0,0)=std::complex<double>(re00,im00);
        selfenergy_orig(i,j,j,0,1)=std::complex<double>(re01,im01);
        selfenergy_orig(i,j,j,1,0)=std::complex<double>(re10,im10);
        selfenergy_orig(i,j,j,1,1)=std::complex<double>(re11,im11);
      }
    }
  }else{
    int iteration=1;
    for(;;++iteration){
      std::stringstream selfenergy_name; 
      selfenergy_name<<"selfenergy_"<<iteration;
      std::ifstream sigma_file(selfenergy_name.str().c_str());
      if(!sigma_file.is_open()) break;
      std::cout<<"averaging iteration : "<<iteration<<std::endl;;
      double ignored;
      for(int i=0;i<nfreq;++i){
        sigma_file>>ignored; //omega_0 = (2*0+1)*pi/beta, beta=pi/omega_0
        for(int j=0;j<nsites_old;++j){
          double re00,im00,re01,im01,re10,im10,re11,im11;
          sigma_file>>re00>>im00>>re01>>im01>>re10>>im10>>re11>>im11;
          selfenergy_orig(i,j,j,0,0)+=std::complex<double>(re00,im00);
          selfenergy_orig(i,j,j,0,1)+=std::complex<double>(re01,im01);
          selfenergy_orig(i,j,j,1,0)+=std::complex<double>(re10,im10);
          selfenergy_orig(i,j,j,1,1)+=std::complex<double>(re11,im11);
        }
      }
    }
    iteration--;
    for(int i=0;i<nfreq;++i){
      for(int j=0;j<nsites_old;++j){
        selfenergy_orig(i,j,j,0,0)/=iteration;
        selfenergy_orig(i,j,j,0,1)/=iteration;
        selfenergy_orig(i,j,j,1,0)/=iteration;
        selfenergy_orig(i,j,j,1,1)/=iteration;
      }
    } 
    std::cout<<"averaging: "<<iteration<<" iterations."<<std::endl;
    std::ofstream sigma_orig_av_file("sigma_orig_av.dat");
    sigma_orig_av_file<<std::make_pair(selfenergy_orig,beta)<<std::endl;
  }
  if(vm.count("cluster")){
    ClusterTransformer* clusterhandler_old=ClusterTransformer::generate_transformer(parms);
    
    if(new_cluster_sites==-1) throw std::logic_error("please specify new cluster sites!");
    parms["SITES"]=new_cluster_sites;
    parms["LATTICE"]=clustername;
    ClusterTransformer *clusterhandler_new=ClusterTransformer::generate_transformer(parms);
    std::vector<std::complex<double> > values(nsites_old);
    std::vector<double> kvec(2);
    for(int i=0;i<nfreq;++i){
      for(int p=0;p<2;++p){
        for(int q=0;q<2;++q){
          for(int k=0;k<nsites_old;++k){
            values[k]=selfenergy_orig(i,k,k,p,q);
          }
          for(int k=0;k<nsites_new;++k){
            for (unsigned int d=0; d<2; ++d){
              kvec[d] = clusterhandler_new->cluster_momentum(k,d);
            }
            selfenergy_kspace_interp(i,k,k,p,q)=clusterhandler_old->interpolate(kvec,values, false);
          }
        }
      }
    }
    delete clusterhandler_old;
    delete clusterhandler_new;
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
    for(int p=0;p<2;++p){
      for(int q=0;q<2;++q){
        for(int i=0;i<nfreq;++i){
          sigma_real[i]=selfenergy_kspace_interp(i,j,j,p,q).real();
          sigma_imag[i]=selfenergy_kspace_interp(i,j,j,p,q).imag();
        }
        gsl_spline_init (spline_real, &(matsubara_freqs_old[0]), &(sigma_real[0]), nfreq);
        gsl_spline_init (spline_imag, &(matsubara_freqs_old[0]), &(sigma_imag[0]), nfreq);
        
        double xi,yi;                     
        for (int i=0; i<nfreq_new; ++i){
          if(i<nfreq){
            selfenergy_new(i,j,j,p,q) = std::complex<double>(gsl_spline_eval (spline_real, matsubara_freqs_new[i], acc_real), gsl_spline_eval (spline_imag, matsubara_freqs_new[i], acc_imag));
          }else{
            selfenergy_new(i,j,j,p,q) = selfenergy_kspace_interp(nfreq-1,j,j,p,q); //TODO: append tail from high freq expansion.
          }
        }
        /*if(selfenergy_new(0,j,j,p,q).real()*selfenergy_new(nfreq-1,j,j,p,q).real()<0.){ //there is a sign change in the real part of the self energy
          std::cout<<"we have a sign change in the real part of the self energy: "<<selfenergy_new(0,j,j,p,q)<<" "<<selfenergy_new(nfreq-1,j,j,p,q)<<" k-vector: "<<j<<" index: "<<p<<" "<<q<<" product: "<<selfenergy_new(0,j,j,p,q).real()*selfenergy_new(nfreq-1,j,j,p,q).real()<<std::endl;
          if(selfenergy_new(0,j,j,p,q).real()*selfenergy_new(nfreq-1,j,j,p,q).real()<-1.e-3){ //if we have a 'significant' deviation.
            for(int i=0;i<nfreq_new;++i){
              selfenergy_new(i,j,j,p,q)=selfenergy_new(nfreq-1,j,j,p,q).real();
            }
          }
        }*/
        /*if(selfenergy_new(0,j,j,p,q).imag()*selfenergy_new(nfreq-1,j,j,p,q).imag()<0.){ //there is a sign change in the real part of the self energy
          std::cout<<"we have a sign change in the imag part of the self energy: "<<selfenergy_new(0,j,j,p,q)<<" "<<selfenergy_new(nfreq-1,j,j,p,q)<<" k-vector: "<<j<<" index: "<<p<<" "<<q<<" product: "<<selfenergy_new(0,j,j,p,q).real()*selfenergy_new(nfreq-1,j,j,p,q).real()<<std::endl;
        }*/
        if(p != q && mult_anomalous!=1){
          for(int i=0;i<nfreq;++i){
            selfenergy_new(i,j,j,p,q)*=mult_anomalous;
          }
        }
        /*        for(int i=0;i<nfreq;++i){
         sigma_real[i]=selfenergy(i,j,j,p,q).real();
         sigma_imag[i]=selfenergy(i,j,j,p,q).imag();
         }
         gsl_spline_init (spline_real, &(matsubara_freqs_old[0]), &(sigma_real[0]), nfreq);
         gsl_spline_init (spline_imag, &(matsubara_freqs_old[0]), &(sigma_imag[0]), nfreq);
         
         double xi,yi;                     
         double selfenergy_real_sign=selfenergy(0,j,j,p,q).real()>0?1:-1;
         double selfenergy_imag_sign=selfenergy(0,j,j,p,q).imag()>0?1:-1;
         for (int i=0; i<nfreq; ++i){
         selfenergy_new(i,j,j,p,q) = std::complex<double>(gsl_spline_eval (spline_real, matsubara_freqs_new[i], acc_real), gsl_spline_eval (spline_imag, matsubara_freqs_new[i], acc_imag));
         if(selfenergy_new(i,j,j,p,q).imag()>0 && selfenergy_imag_sign<0){
         std::cerr<<"problem with self energy sign after interpolation."<<i<<" "<<j<<" "<<p<<" "<<q<<std::endl;
         }else if(selfenergy_new(i,j,j,p,q).real()>0 && selfenergy_real_sign<0){
         std::cerr<<"problem with self energy real after interpolation."<<i<<" "<<j<<" "<<p<<" "<<q<<std::endl;
         }else if(selfenergy_new(i,j,j,p,q).real()<0 && selfenergy_real_sign>0){
         std::cerr<<"problem with self energy real after interpolation."<<i<<" "<<j<<" "<<p<<" "<<q<<std::endl;
         }
         }*/
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
  parms["ETA1"]=eta1_new;
  SimpleDCATransformer clustertrans(parms);
  green0 = clustertrans(green, green0, selfenergy_new, mu_new, h, eta1_new, eta2, beta_new, true);  
  
  green0 = clustertrans.symmetrize(green0, paramagnet);
  green0 = clustertrans.symmetrize_k_space(green0);
  
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  fourier_ptr->frequency_to_time_ft(green0_itime, green0);
  std::ofstream green0_interpolated_tau_file("green0_interpolated_tau.dat");
  
  green0_interpolated_tau_file<<std::make_pair(green0_itime,beta_new)<<std::endl;
  std::stringstream filename_new; filename_new<<"G0_omega_input_beta"<<beta_new<<"_mu"<<mu_new<<"_eta"<<eta1_new;
  std::stringstream filename_sigma_new; filename_sigma_new<<"sigma_beta"<<beta_new;
  green0.write(filename_new.str().c_str());
  std::ofstream sigma_file(filename_sigma_new.str().c_str());
  sigma_file<<std::make_pair(selfenergy_new,beta)<<std::endl;
}

