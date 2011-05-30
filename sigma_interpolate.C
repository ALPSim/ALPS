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


int main(int argc, char **argv){
  if(argc !=8) throw(std::invalid_argument("start up program with ./program_name parm_file lattice_file iteration new_beta new_mu new_eta1 new_n_matsubara"));
  std::ifstream parm_file(argv[1]);
  if(!parm_file.is_open()) throw std::runtime_error("parm file not open.");
  alps::Parameters parms(parm_file);
  double beta=parms["BETA"];
  double mu=parms["MU"];
  double h=parms["H"];
  double eta1=parms["ETA1"];
  double eta2=parms["ETA2"];
  double nsites=parms["SITES"];
  int nfreq = parms["NMATSUBARA"];
  int ntau= parms["N"];
  double beta_new=atof(argv[4]);
  double mu_new=atof(argv[5]);
  double eta1_new=atof(argv[6]);
  int nfreq_new=atoi(argv[7]);
  bool paramagnet=true;
  
  std::cout<<"adjusting parameters: "<<std::endl;
  std::cout<<"new beta: "<<beta_new<<"(old: "<<beta<<")"<<std::endl;
  std::cout<<"new mu: "<<mu_new<<"(old: "<<mu<<")"<<std::endl;
  std::cout<<"new eta: "<<eta1_new<<"(old: "<<eta1<<")"<<std::endl;
  
  parms["LATTICE_LIBRARY"]=std::string(argv[2]);
  
  
  nambu_matsubara_green_function_t selfenergy(nfreq, nsites, 2, diagonal);
  nambu_matsubara_green_function_t selfenergy_new(nfreq_new, nsites, 2, diagonal);
  nambu_matsubara_green_function_t green(nfreq, nsites, 2, diagonal);
  nambu_matsubara_green_function_t green0(nfreq_new, nsites, 2, diagonal);
  nambu_itime_green_function_t green0_itime(ntau, nsites, 2, diagonal);
  
  int iteration=atoi(argv[3]);
  std::stringstream selfenergy_name; selfenergy_name<<"selfenergy_"<<iteration;
  {
    std::ifstream sigma_file(selfenergy_name.str().c_str());
    if(!sigma_file.is_open()) throw std::runtime_error("self energy file not open.");
    double ignored;
    for(int i=0;i<nfreq;++i){
      sigma_file>>ignored; //omega_0 = (2*0+1)*pi/beta, beta=pi/omega_0
      for(int j=0;j<nsites;++j){
        double re00,im00,re01,im01,re10,im10,re11,im11;
        sigma_file>>re00>>im00>>re01>>im01>>re10>>im10>>re11>>im11;
        selfenergy(i,j,j,0,0)=std::complex<double>(re00,im00);
        selfenergy(i,j,j,0,1)=std::complex<double>(re01,im01);
        selfenergy(i,j,j,1,0)=std::complex<double>(re10,im10);
        selfenergy(i,j,j,1,1)=std::complex<double>(re11,im11);
      }
    }
  }
  //std::cout<<std::make_pair(selfenergy,beta)<<std::endl;
  
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
  for(int j=0;j<nsites;++j){
    for(int p=0;p<2;++p){
      for(int q=0;q<2;++q){
        for(int i=0;i<nfreq;++i){
          sigma_real[i]=selfenergy(i,j,j,p,q).real();
          sigma_imag[i]=selfenergy(i,j,j,p,q).imag();
        }
        gsl_spline_init (spline_real, &(matsubara_freqs_old[0]), &(sigma_real[0]), nfreq);
        gsl_spline_init (spline_imag, &(matsubara_freqs_old[0]), &(sigma_imag[0]), nfreq);
        
        double xi,yi;                     
        for (int i=0; i<nfreq_new; ++i){
          if(i<nfreq){
            selfenergy_new(i,j,j,p,q) = std::complex<double>(gsl_spline_eval (spline_real, matsubara_freqs_new[i], acc_real), gsl_spline_eval (spline_imag, matsubara_freqs_new[i], acc_imag));
          }else{
            selfenergy_new(i,j,j,p,q) = selfenergy(nfreq-1,j,j,p,q); //TODO: append tail from high freq expansion.
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

