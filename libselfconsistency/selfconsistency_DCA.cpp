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
#include <iomanip>

template<class GFW, class GFT> void print_all_green_functions(int iteration_ctr,
    const GFW &G0_omega,
    const GFW &G_omega,
    const GFW &selfenergy,  //the self-energy computed with the old G0
    const GFW &selfenergy2, //the self-energy computed with the new G0
    const GFT &G0_tau,
    const GFT &G_tau, double beta){

  std::ostringstream G0omega_name, G0tau_name, Gomega_name,
  Gtau_name, Gtau_name_2, Gomega_name_2, selfenergy_name, selfenergy_name_newG0, selfenergy_name2;
  G0omega_name<<"G0_omega"<<"_"<<iteration_ctr;
  G0tau_name<<"G0_tau"<<"_"<<iteration_ctr;
  Gomega_name<<"G_omega"<<"_"<<iteration_ctr;
  Gtau_name<<"G_tau"<<"_"<<iteration_ctr;
  Gomega_name_2<<"G_omega";
  Gtau_name_2<<"G_tau";
  selfenergy_name<<"selfenergy"<<"_"<<iteration_ctr;
  selfenergy_name_newG0<<"selfenergy_newG0"<<"_"<<iteration_ctr;
  selfenergy_name2<<"selfenergy"<<"_"<<iteration_ctr<<".h5";

  std::ofstream G0_omega_file(G0omega_name.str().c_str());
  std::ofstream G0_tau_file(G0tau_name.str().c_str());
  std::ofstream G_omega_file(Gomega_name.str().c_str());
  std::ofstream G_tau_file(Gtau_name.str().c_str());
  std::ofstream G_omega_file_2(Gomega_name_2.str().c_str());
  std::ofstream G_tau_file_2(Gtau_name_2.str().c_str());
  std::ofstream selfenergy_file(selfenergy_name.str().c_str());
  std::ofstream selfenergy_file_newG0(selfenergy_name_newG0.str().c_str());

  G0_omega_file<<std::setprecision(14)<<G0_omega;
  G0_tau_file<<std::setprecision(14)<<G0_tau;
  G_omega_file<<std::setprecision(14)<<G_omega;
  G_tau_file<<std::setprecision(14)<<G_tau;
  G_omega_file_2<<std::setprecision(14)<<G_omega;
  G_tau_file_2<<std::setprecision(14)<<G_tau;
  selfenergy_file<<std::setprecision(14)<<selfenergy;  
  selfenergy_file_newG0<<std::setprecision(14)<<selfenergy2;

  alps::hdf5::archive ars(selfenergy_name2.str(), alps::hdf5::archive::WRITE);
  selfenergy.save(ars, "/Sigma");
}
template<class GFW> bool converged(const GFW &G_omega, const GFW &G_omega_old, double tolerance){
  return ((G_omega-G_omega_old).norm() > tolerance)?false:true;
}
template<class GFW> void write_new_G0omega(const GFW &G0_omega, const std::string &path){
  alps::hdf5::archive ar(path, alps::hdf5::archive::WRITE);
  G0_omega.save(ar, "/G0");
}

template<class SOLV, class TRANS, class GFW, class GFT> void selfconsistency_loop_DCA(const alps::params& parms, SOLV& solver, TRANS& clustertrans, GFW &G0_omega, GFT &G0_tau)
{
  double beta=parms["BETA"];

  int max_it = parms["sc.MAX_IT"];
  double convergence_tolerance=parms["sc.SELFCONSISTENCY_CONVERGENCE"];
  std::string g0_output_path=parms["dca.G0_OMEGA_OUTPUT"].as<std::string>();

  //define all Matsubara and imaginary time Green's functions needed
  GFW G_omega=G0_omega;
  GFW G_omega_old=G0_omega;
  GFW G_omega_ksum=G0_omega;
  GFW G0_omega_new=G0_omega;
  GFW self_energy=G0_omega;
  GFW self_energy2=G0_omega;

  GFT G_tau(G0_tau);

  clustertrans.initial_G0(G0_omega, G0_tau);

  //This is the basic DCA self-consistency loop
  for(int iteration_ctr=0;iteration_ctr < max_it; ++iteration_ctr) {
    alps::gf::fourier_frequency_to_time(G0_omega, G0_tau);
    solver.solve(G0_omega, G0_tau, G_omega, G_tau);

    clustertrans.spin_symmetrize(G_omega);
    clustertrans.symmetrize_k_space(G_omega);
    clustertrans.dyson_get_self_energy(self_energy, G_omega, G0_omega);
    clustertrans.matsubara_k_sum(G_omega_ksum, self_energy);
    clustertrans.dyson_get_green0(G0_omega_new, self_energy, G_omega_ksum);
    clustertrans.dyson_get_self_energy(self_energy2, G_omega, G0_omega);

    print_all_green_functions(iteration_ctr, G0_omega, G_omega, self_energy, self_energy2, G0_tau, G_tau, beta);
    if (converged(G_omega, G_omega_old,  convergence_tolerance)) break;

    G_omega_old=G_omega;
    G0_omega=G0_omega_new;
    write_new_G0omega(G0_omega_new, g0_output_path);
  }
}
void perform_cluster_loop(const std::string &param_file_name, const alps::params &parms){
  std::string cluster_loop=parms["SC_LOOP"];
  bool nambu=parms["sc.NAMBU"];
  std::string solver_path=parms["solver.SOLVER"];
  int n_tau=parms["N"];
  int n_matsubara=parms["NMATSUBARA"];
  int n_orbital=parms["FLAVORS"];
  double beta=parms["BETA"];

  std::cout<<"cluster loop is: "<<cluster_loop<<std::endl;
  if (!nambu && cluster_loop=="DCA"){
    std::cout<<"running a DCA calculation in the normal state"<<std::endl;
    NormalStateClusterTransformer T(parms);
    ImpuritySolver solver(solver_path,param_file_name, parms);
    cluster_matsubara_kspace_gf G0_omega(alps::gf::omega_k_sigma_gf(alps::gf::matsubara_positive_mesh(beta,n_matsubara),T.momentum_mesh(), alps::gf::index_mesh(n_orbital)));
    cluster_itime_kspace_gf G0_tau(alps::gf::itime_k_sigma_gf(alps::gf::itime_mesh(beta,n_tau),T.momentum_mesh(), alps::gf::index_mesh(n_orbital)));
    selfconsistency_loop_DCA(parms, solver, T,G0_omega,G0_tau);
  }else if(cluster_loop=="DCA"){
    std::cout<<"running a DCA calculation allowing for superconducting order"<<std::endl;
    NambuClusterTransformer T(parms);
    ImpuritySolver solver(solver_path,param_file_name, parms);
    cluster_nambu_matsubara_kspace_gf G0_omega(alps::gf::omega_k_sigma1_sigma2_gf(alps::gf::matsubara_positive_mesh(beta,n_matsubara),T.momentum_mesh(), alps::gf::index_mesh(n_orbital),alps::gf::index_mesh(n_orbital)));
    cluster_nambu_itime_kspace_gf     G0_tau(alps::gf::itime_k_sigma1_sigma2_gf  (alps::gf::itime_mesh(beta,n_tau),T.momentum_mesh(), alps::gf::index_mesh(n_orbital),alps::gf::index_mesh(n_orbital)));
    selfconsistency_loop_DCA(parms, solver, T, G0_omega,G0_tau);
  }else{
    boost::throw_exception(std::logic_error("there is no such cluster method: "+cluster_loop));
  }
}




