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


///This program takes a self-energy from a file, an old parameter file that specifies for which parameters it has been generated, and a new parameter file that has the new parameters.
///It then generates a few
int main(int argc, char **argv){
  namespace po = boost::program_options;
  std::string pfilename_new,pfilename_old;
  std::string latticefilename;
  std::string clustername;
  std::string selfenergy_filename;
  int new_cluster_sites=-1;
  // Declare the supported options.
  double mult_anomalous;
  po::options_description desc("Allowed options");
  desc.add_options()
  ("help", "show this help")
  ("pfile_old", po::value<std::string>(&pfilename_old), "parameter file of the old parameter file")
  ("pfile_new", po::value<std::string>(&pfilename_new), "parameter file of the new parameter file")
  ("selfenergy", po::value<std::string>(&selfenergy_filename), "self-energy to be interpolated")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout<<desc;
    return 1;
  }
  //parse parameters: error check
  if(!vm.count("pfile_old")) throw std::runtime_error("you need to specify the old parameter file.");
  if(!vm.count("pfile_new")) throw std::runtime_error("you need to specify the new parameter file.");
  if(!vm.count("selfenergy")) throw std::runtime_error("you need to specify the selfenergy.");
  
  //parse parameters: open old parameter file
  std::ifstream parm_file_old(pfilename_old.c_str());
  if(!parm_file_old.is_open()) throw std::runtime_error("old parm file not open.");
  alps::Parameters parms_old(parm_file_old);
  
  //parse parameters: open new parameter file
  std::ifstream parm_file_new(pfilename_new.c_str());
  if(!parm_file_new.is_open()) throw std::runtime_error("new parm file not open.");
  alps::Parameters parms_new(parm_file_new);
  
  //general error check
  if(parms_old["CLUSTER_LOOP"]!="DCA") throw std::invalid_argument("this interpolation assumes DCA");
  double beta_old=parms_old["BETA"];
  double mu_old=parms_old["MU"];
  double h_old=parms_old["H"];
  double nsites_old=parms_old["SITES"];
  int nfreq_old = parms_new["NMATSUBARA"];
  
  double beta_new=parms_new["BETA"];
  double mu_new=parms_new["MU"];
  double h_new=parms_new["H"];
  int nfreq_new=parms_new["NMATSUBARA"];
  double nsites_new=parms_new["SITES"];
  
  bool paramagnet=!(parms_new.value_or_default("ANTIFERROMAGNET",false));
  shape_t gf_shape=paramagnet?diagonal:blockdiagonal;
  
  std::cout<<"main parameters: "<<std::endl;
  std::cout<<"new beta: "<<beta_new<<"(old: "<<beta_old<<")"<<std::endl;
  std::cout<<"new mu: "<<mu_new<<"(old: "<<mu_old<<")"<<std::endl;
  std::cout<<"new h: "<<h_new<<"(old: "<<h_old<<")"<<std::endl;
  
  //definition of Green's functions
  matsubara_green_function_t selfenergy_orig(nfreq_old, nsites_old, 2, gf_shape);
  matsubara_green_function_t selfenergy_kspace_interp(nfreq_old, nsites_new, 2, gf_shape);
  matsubara_green_function_t selfenergy_new(nfreq_new, nsites_new, 2, gf_shape);
  matsubara_green_function_t green(nfreq_old, nsites_new, 2, gf_shape);
  matsubara_green_function_t green0(nfreq_new, nsites_new, 2, gf_shape);
  
  //open and read in self-energy
  if(!vm.count("selfenergy")) throw std::runtime_error("please specify self-energy");
  {
    double ignored;
    std::ifstream sigma_file(selfenergy_filename.c_str());
    sigma_file>>selfenergy_orig;
    std::ofstream s_file("sigma_orig.dat"); s_file<<std::make_pair(selfenergy_orig, beta_old);
  }
  
  //interpolate cluster geometry in k-space
  if(std::string(parms_old["CLUSTER"]).compare(std::string(parms_new["CLUSTER"]))!=0){
    std::cout<<"...interpolating momentum points..."<<std::endl;
    DCATransformer* clustertransformer_old=(DCATransformer*)ClusterTransformer::generate_transformer(parms_old);
    DCATransformer* clustertransformer_new=(DCATransformer*)ClusterTransformer::generate_transformer(parms_new);

    for(int i=0;i<nfreq_old;++i){
      for(int flavor=0;flavor<2;++flavor){
        std::vector<std::complex<double> > values(nsites_old);
        for(int k=0;k<nsites_old;++k){
          values[k]=selfenergy_orig(i,k,k,flavor);
        }
        for(int k=0;k<nsites_new;++k){
          std::vector<double> kvec(2);
          for (unsigned int d=0; d<2; ++d){
            kvec[d] = clustertransformer_new->cluster_momentum(k,d);
          }
          selfenergy_kspace_interp(i,k,k,flavor)=clustertransformer_old->interpolate(kvec,values, false);
        }
      }
      if(!paramagnet){
        for(int flavor=0;flavor<2;++flavor){
          std::vector<std::complex<double> > values_kq(nsites_old);
          std::vector<std::complex<double> > values_qk(nsites_old);
          for(int k=0;k<nsites_old;k+=2){
            values_kq[k]=selfenergy_orig(i,k,k+1,flavor);
            values_qk[k]=selfenergy_orig(i,k+1,k,flavor);
          }
          for(int k=0;k<nsites_new;++k){
            std::vector<double> kvec(2);
            for (unsigned int d=0; d<2; ++d){
              kvec[d] = clustertransformer_new->cluster_momentum(k,d);
            }
            selfenergy_kspace_interp(i,k,k+1,flavor)=clustertransformer_old->interpolate(kvec,values_kq, false);
            selfenergy_kspace_interp(i,k+1,k,flavor)=clustertransformer_old->interpolate(kvec,values_qk, false);
          }
        }
      }
    }
    delete clustertransformer_old;
    delete clustertransformer_new;
  }else{
    selfenergy_kspace_interp=selfenergy_orig;
    std::ofstream s_file("sigma_k_int.dat"); s_file<<std::make_pair(selfenergy_kspace_interp, beta_old);
  }
  
  //interpolate for new temperature
  double nfreq_min=std::min(nfreq_old, nfreq_new);
  std::vector<double> matsubara_freqs_old(nfreq_min);
  std::vector<double> matsubara_freqs_new(nfreq_min);
  for(int i=0;i<nfreq_min;++i){
    matsubara_freqs_old[i]=(2.*i+1)*M_PI/beta_old;
    matsubara_freqs_new[i]=(2.*i+1)*M_PI/beta_new;
  }

  gsl_interp_accel *acc_real = gsl_interp_accel_alloc ();
  gsl_spline *spline_real = gsl_spline_alloc (gsl_interp_cspline, nfreq_min);
  gsl_interp_accel *acc_imag = gsl_interp_accel_alloc ();
  gsl_spline *spline_imag = gsl_spline_alloc (gsl_interp_cspline, nfreq_min);

  gsl_interp_accel *acc_real_kq = gsl_interp_accel_alloc ();
  gsl_spline *spline_real_kq = gsl_spline_alloc (gsl_interp_cspline, nfreq_min);
  gsl_interp_accel *acc_imag_kq = gsl_interp_accel_alloc ();
  gsl_spline *spline_imag_kq = gsl_spline_alloc (gsl_interp_cspline, nfreq_min);

  gsl_interp_accel *acc_real_qk = gsl_interp_accel_alloc ();
  gsl_spline *spline_real_qk = gsl_spline_alloc (gsl_interp_cspline, nfreq_min);
  gsl_interp_accel *acc_imag_qk = gsl_interp_accel_alloc ();
  gsl_spline *spline_imag_qk = gsl_spline_alloc (gsl_interp_cspline, nfreq_min);

  for(int j=0;j<nsites_new;++j){
    for(int flavor=0;flavor<2;++flavor){
      std::vector<double> sigma_real(nfreq_min);
      std::vector<double> sigma_imag(nfreq_min);
      for(int i=0;i<nfreq_min;++i){
        sigma_real[i]=selfenergy_kspace_interp(i,j,j,flavor).real();
        sigma_imag[i]=selfenergy_kspace_interp(i,j,j,flavor).imag();
      }
      gsl_spline_init (spline_real, &(matsubara_freqs_old[0]), &(sigma_real[0]), nfreq_min);
      gsl_spline_init (spline_imag, &(matsubara_freqs_old[0]), &(sigma_imag[0]), nfreq_min);
      
      double xi,yi;
      for (int i=0; i<nfreq_min; ++i){
        selfenergy_new(i,j,j,flavor) = std::complex<double>(gsl_spline_eval (spline_real, matsubara_freqs_new[i], acc_real), gsl_spline_eval (spline_imag, matsubara_freqs_new[i], acc_imag));
      }
      for(int i=nfreq_min;i<nfreq_new;++i){
        selfenergy_new(i,j,j,flavor) = selfenergy_kspace_interp(nfreq_old-1,j,j,flavor); //TODO: append tail from high freq expansion.
      }
    }
  }
  if(!paramagnet){
    for(int j=0;j<nsites_new;j+=2){
      for(int flavor=0;flavor<2;++flavor){
        std::vector<double> sigma_real_kq(nfreq_min);
        std::vector<double> sigma_imag_kq(nfreq_min);
        std::vector<double> sigma_real_qk(nfreq_min);
        std::vector<double> sigma_imag_qk(nfreq_min);
        for(int i=0;i<nfreq_min;++i){
          sigma_real_kq[i]=selfenergy_kspace_interp(i,j,j+1,flavor).real();
          sigma_imag_kq[i]=selfenergy_kspace_interp(i,j,j+1,flavor).imag();
          sigma_real_qk[i]=selfenergy_kspace_interp(i,j+1,j,flavor).real();
          sigma_imag_qk[i]=selfenergy_kspace_interp(i,j+1,j,flavor).imag();
        }
        gsl_spline_init (spline_real_kq, &(matsubara_freqs_old[0]), &(sigma_real_kq[0]), nfreq_min);
        gsl_spline_init (spline_imag_kq, &(matsubara_freqs_old[0]), &(sigma_imag_kq[0]), nfreq_min);
        gsl_spline_init (spline_real_qk, &(matsubara_freqs_old[0]), &(sigma_real_kq[0]), nfreq_min);
        gsl_spline_init (spline_imag_qk, &(matsubara_freqs_old[0]), &(sigma_imag_kq[0]), nfreq_min);
        
        double xi,yi;
        for (int i=0; i<nfreq_min; ++i){
          selfenergy_new(i,j,j+1,flavor) = std::complex<double>(gsl_spline_eval (spline_real_kq, matsubara_freqs_new[i], acc_real), gsl_spline_eval (spline_imag_kq, matsubara_freqs_new[i], acc_imag));
          selfenergy_new(i,j+1,j,flavor) = std::complex<double>(gsl_spline_eval (spline_real_qk, matsubara_freqs_new[i], acc_real), gsl_spline_eval (spline_imag_qk, matsubara_freqs_new[i], acc_imag));
        }
        for(int i=nfreq_min;i<nfreq_new;++i){
          selfenergy_new(i,j,j+1,flavor) = selfenergy_kspace_interp(nfreq_old-1,j,j+1,flavor); //TODO: append tail from high freq expansion.
          selfenergy_new(i,j+1,j,flavor) = selfenergy_kspace_interp(nfreq_old-1,j+1,j,flavor); //TODO: append tail from high freq expansion.
        }
      }
    }
  }
  gsl_spline_free (spline_real);
  gsl_spline_free (spline_imag);
  gsl_interp_accel_free (acc_real);
  gsl_interp_accel_free (acc_imag);
  
  gsl_spline_free (spline_real_kq);
  gsl_spline_free (spline_imag_kq);
  gsl_interp_accel_free (acc_real_kq);
  gsl_interp_accel_free (acc_imag_kq);
  gsl_spline_free (spline_real_qk);
  gsl_spline_free (spline_imag_qk);
  gsl_interp_accel_free (acc_real_qk);
  gsl_interp_accel_free (acc_imag_qk);

  //compute the self consistency: G0 from the interpolated sigma
  std::cout<<"...running self-consistency..."<<std::endl;
  {
    std::ofstream s_file("sigma_in.dat"); s_file<<std::make_pair(selfenergy_new, beta_new);
    std::ofstream g0_file("green0_in.dat"); g0_file<<std::make_pair(green0, beta_new);
    std::ofstream g_file("green_in.dat"); g_file<<std::make_pair(green, beta_new);
  }
  DCATransformer* clustertransformer_new=(DCATransformer*)ClusterTransformer::generate_transformer(parms_new);
  green0 = (*clustertransformer_new)(green, green0, selfenergy_new, mu_new, h_new, beta_new, true);
  {
    std::ofstream s_file("sigma_out.dat"); s_file<<std::make_pair(selfenergy_new, beta_new);
  }
  if(paramagnet){
    green0 = clustertransformer_new->symmetrize(green0, paramagnet);
    green0 = clustertransformer_new->symmetrize_k_space(green0);
  }
  
  std::string filename_new=parms_new["G0OMEGA_INPUT"];
  if(std::ifstream(filename_new.c_str(),std::ifstream::in)){ throw std::runtime_error("input file "+std::string(parms_new["G0OMEGA_INPUT"])+"already exists!"); }
  green0.write(filename_new.c_str());
  
  delete clustertransformer_new;
}

