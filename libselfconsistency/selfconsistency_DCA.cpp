/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2014 by Emanuel Gull <gullc@comp-phys.org>
 *
 * To be used with permission by the authors
 *
 *****************************************************************************/

#include <sys/types.h> 

#include <boost/tuple/tuple.hpp>

#include <alps/params.hpp>
#include <alps/gf/fourier.hpp>

#include"cluster_transform.h"
#include"impurity_solver.h"


void print_all_green_functions(int iteration_ctr,
    const cluster_matsubara_kspace_gf &G0_omega,
    const cluster_matsubara_kspace_gf &G_omega,
    const cluster_matsubara_kspace_gf &selfenergy,
    const cluster_itime_kspace_gf &G0_tau,
    const cluster_itime_kspace_gf &G_tau, double beta){

  std::ostringstream G0omega_name, G0tau_name, Gomega_name,
  Gtau_name, Gtau_name_2, Gomega_name_2, selfenergy_name;
  G0omega_name<<"G0_omega"<<"_"<<iteration_ctr;
  G0tau_name<<"G0_tau"<<"_"<<iteration_ctr;
  Gomega_name<<"G_omega"<<"_"<<iteration_ctr;
  Gtau_name<<"G_tau"<<"_"<<iteration_ctr;
  Gomega_name_2<<"G_omega";
  Gtau_name_2<<"G_tau";
  selfenergy_name<<"selfenergy"<<"_"<<iteration_ctr;

  std::ofstream G0_omega_file(G0omega_name.str().c_str());
  std::ofstream G0_tau_file(G0tau_name.str().c_str());
  std::ofstream G_omega_file(Gomega_name.str().c_str());
  std::ofstream G_tau_file(Gtau_name.str().c_str());
  std::ofstream G_omega_file_2(Gomega_name_2.str().c_str());
  std::ofstream G_tau_file_2(Gtau_name_2.str().c_str());
  std::ofstream selfenergy_file(selfenergy_name.str().c_str());

  G0_omega_file<<G0_omega;
  G0_tau_file<<G0_tau;
  G_omega_file<<G_omega;
  G_tau_file<<G_tau;
  G_omega_file_2<<G_omega;
  G_tau_file_2<<G_tau;
  selfenergy_file<<selfenergy;
}
bool converged(const cluster_matsubara_kspace_gf &G_omega, const cluster_matsubara_kspace_gf &G_omega_old, double tolerance){
  double max_diff=0;
  return ((G_omega-G_omega_old).norm() > max_diff)?false:true;
}

void selfconsistency_loop_DCA(const alps::params& parms, ImpuritySolver& solver,
    NormalStateClusterTransformer& clustertrans)
{
  int n_tau=parms["dca.N"];
  int n_orbital=parms["dca.FLAVORS"];
  double beta=parms["dca.BETA"];
  
  int max_it = parms["sc.MAX_IT"];
  double convergence_tolerance=parms["sc.SELFCONSISTENCY_CONVERGENCE"];

  //calculate initial G0
  cluster_matsubara_kspace_gf G0_omega=clustertrans.initial_G0();

  //define all needed matsubara and imaginary time Green's functions
  cluster_matsubara_kspace_gf G_omega=G0_omega;
  cluster_matsubara_kspace_gf G_omega_old=G0_omega;
  cluster_matsubara_kspace_gf G_omega_ksum=G0_omega;
  cluster_matsubara_kspace_gf G0_omega_new=G0_omega;
  cluster_matsubara_kspace_gf self_energy=G0_omega;

  cluster_itime_kspace_gf G_tau(alps::gf::itime_k_sigma_gf(alps::gf::itime_mesh(beta,n_tau),clustertrans.momentum_mesh(), alps::gf::index_mesh(n_orbital)));
  cluster_itime_kspace_gf G0_tau(alps::gf::itime_k_sigma_gf(alps::gf::itime_mesh(beta,n_tau),clustertrans.momentum_mesh(), alps::gf::index_mesh(n_orbital)));


  //This is the basic DCA self-consistency loop
  for(int iteration_ctr=0;iteration_ctr < max_it; ++iteration_ctr) {
    std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;
    alps::gf::fourier_frequency_to_time(G0_omega, G0_tau);
    solver.solve(G0_omega, G0_tau, G_omega, G_tau);

    clustertrans.spin_symmetrize(G_omega);
    clustertrans.symmetrize_k_space(G_omega);
    clustertrans.dyson_get_self_energy(self_energy, G_omega, G0_omega);
    clustertrans.matsubara_k_sum(G_omega_ksum, self_energy);
    clustertrans.dyson_get_green0(G0_omega_new, self_energy, G_omega_ksum);


    print_all_green_functions(iteration_ctr, G0_omega, G_omega, self_energy, G0_tau, G_tau, beta);
    if (converged(G_omega, G_omega_old,  convergence_tolerance)) break;

    G_omega_old=G_omega;
    G0_omega=G0_omega_new;
  }
}

void perform_cluster_loop(const alps::params &parms){
  std::string cluster_loop=parms["SC_LOOP"];
  std::string solver_path=parms["solver.SOLVER"];
  std::cout<<"cluster loop is: "<<cluster_loop<<std::endl;
  if (cluster_loop=="DCA"){
    std::cout<<"running a DCA calculation"<<std::endl;
    NormalStateClusterTransformer T(parms);
    ImpuritySolver solver(solver_path,parms);
    selfconsistency_loop_DCA(parms, solver, T);
  }else {
    boost::throw_exception(std::logic_error("there is no such cluster method: "+cluster_loop));
  }
}


/*void selfconsistency_loop_DCA_Nambu(const alps::params& parms, ImpuritySolver& solver,
                              DCATransformer& clustertrans)
{
  unsigned int n_tau=parms["N"];
  unsigned int n_matsubara=parms["NMATSUBARA"];
  unsigned int n_site=parms["SITES"]|1;

  double beta = parms["BETA"];
  double h    = parms["H"];
  double eta1 = parms["ETA1"];
  double eta2 = parms["ETA2"];
  double converged = parms["CONVERGED"];
  bool paramagnet = parms["PARAMAGNET"]|true;
  if(!paramagnet) throw std::logic_error("please implement mixed AFM & Nambu!");
  int max_it=parms["MAX_IT"]|1000;
  std::string basename=parms["BASENAME"];


  nambu_matsubara_green_function_t G0_omega=clustertrans.initial_G0_nambu(parms);
  nambu_matsubara_green_function_t G_omega=G0_omega;

  int iteration_ctr = 0;
  nambu_matsubara_green_function_t G_old_omega(G0_omega);
  nambu_matsubara_green_function_t G0_old_omega(G0_omega);
  nambu_matsubara_green_function_t selfenergy(G0_omega);
  nambu_itime_green_function_t G_tau(n_tau +1, n_site, 2, diagonal);
  nambu_itime_green_function_t G0_tau(n_tau+1, n_site, 2, diagonal);
  nambu_itime_green_function_t G_old_tau(n_tau+1, n_site, 2, diagonal);
  nambu_itime_green_function_t G0_old_tau(n_tau+1, n_site, 2, diagonal);
  //initialize space for multiple vectors
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  FourierTransformer::generate_transformer(parms, fourier_ptr);
  double max_diff;
  do {
    iteration_ctr++;
    std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;
    double mu = parms["MU"];
    G_old_omega = G_omega;
    G0_old_omega = G0_omega;
    G_old_tau = G_tau;
    G0_old_tau = G0_tau;
    fourier_ptr->frequency_to_time_ft(G0_tau, G0_omega);
    std::cout<<"running NAMBU solver in omega."<<std::endl<<std::flush;
    //running additional symmetrizations:
    G0_omega = clustertrans.symmetrize_k_space(G0_omega);
    boost::tie(G_omega, G_tau) = solver.solve_nambu(G0_omega, parms);
    G_tau = clustertrans.symmetrize(G_tau, paramagnet);
    G_tau = clustertrans.symmetrize_k_space(G_tau); //<- this would symmetrize away the d-wave behavior!
    G_omega = clustertrans.symmetrize(G_omega, paramagnet);
    G_omega = clustertrans.symmetrize_k_space(G_omega);

    std::cout<<"running cluster transform"<<std::endl<<std::flush;
    G0_omega = clustertrans(G_omega, G0_omega, selfenergy, mu, h, eta1, eta2, beta, false);


    G0_omega.write(std::string(parms["G0OMEGA_output"]|"G0omega_output").c_str());
    std::cout<<"comparing old and new result."<<std::endl<<std::flush;
    max_diff=0;
    if(iteration_ctr>1){
      for(int f1=0;f1<2;++f1){
        for(int f2=0;f2<2;++f2){
          //comparison for the dressed Green's function in Matsubara freq.
          for(unsigned int w=0; w<n_matsubara; w++) {
            for(unsigned int i1=0; i1<n_site; i1++) {
              for(unsigned int i2=0; i2<n_site; i2++) {
                if (std::abs(G_omega(w,i1, i2, f1,f2)-G_old_omega(w, i1, i2, f1,f2)) > max_diff)
                  max_diff = std::abs(G_omega(w,i1,i2,f1,f2)-G_old_omega(w,i1,i2,f1,f2));
              }
            }
          }
        }
      }
      std::cout<<"convergence loop: max diff in dressed Green (Matsubara freq): "
      <<max_diff<<"\tconverged: "<<converged<<std::endl;
    }
    nambu_matsubara_green_function_t real_G0_omega = clustertrans.transform_into_real_space(G0_omega);
    nambu_matsubara_green_function_t real_G_omega = clustertrans.transform_into_real_space(G_omega);
    nambu_itime_green_function_t real_G0_tau = clustertrans.transform_into_real_space(G0_tau);
    nambu_itime_green_function_t real_G_tau = clustertrans.transform_into_real_space(G_tau);
    print_all_nambu_green_functions(basename, iteration_ctr, G0_omega, G_omega, G0_tau, G_tau, beta, diagonal);
    if(!parms.defined("WRITE_SPARSE"))
      print_all_nambu_green_functions(basename, iteration_ctr, real_G0_omega, real_G_omega, real_G0_tau, real_G_tau, beta, nondiagonal, "_realspace");



  } while ((max_diff > converged || iteration_ctr <= 1) && iteration_ctr < max_it);
  // write G0
  G0_omega.write(std::string(parms["G0OMEGA_output"]|"G0omega_output").c_str());
}
 */

