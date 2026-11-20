/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005-2026 by Emanuel Gull
 *
 * ALPS Project: https://alps.comp-phys.org/
 * SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include <iomanip>

#include "alps/params.hpp"
#include "alps/gf/fourier.hpp"

#include "cluster_transform.h"
#include "impurity_solver.h"
#include "selfconsistency.h"

template<class GFW, class GFT> void print_all_green_functions(int iteration_ctr,
    const GFW &G0_omega,
    const GFW &G_omega,
    const GFW &selfenergy,
    const GFT &G0_tau,
    const GFT &G_tau, double beta){

  std::ostringstream G0omega_name, G0tau_name, Gomega_name,
  Gtau_name, Gtau_name_2, Gomega_name_2, selfenergy_name, Greens_h5_name;
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

  G0_omega_file<<std::setprecision(14)<<G0_omega;
  G0_tau_file<<std::setprecision(14)<<G0_tau;
  G_omega_file<<std::setprecision(14)<<G_omega;
  G_tau_file<<std::setprecision(14)<<G_tau;
  G_omega_file_2<<std::setprecision(14)<<G_omega;
  G_tau_file_2<<std::setprecision(14)<<G_tau;
  selfenergy_file<<std::setprecision(14)<<selfenergy;  

  Greens_h5_name << "Greens_" <<iteration_ctr<<".h5";
  alps::hdf5::archive ar(Greens_h5_name.str().c_str(), "w");
  G_omega.save(ar, "/G_komega");
  G0_omega.save(ar, "/G0_komega");
  selfenergy.save(ar, "/Sigma_komega");
  ar.close();
}
template<class GFW> void write_new_G0omega(const GFW &G0_omega, const std::string &path){
  alps::hdf5::archive ar(path, "w");
  G0_omega.save(ar, "/G0");
}

template<class SOLV, class TRANS, class GFW, class GFT> void selfconsistency_loop_DCA(const alps::params& parms, SOLV& solver, TRANS& clustertrans, GFW &G0_omega, GFT &G0_tau)
{
  double beta=parms["BETA"];
  bool spin_symmetrize = parms["sc.SPIN_SYMMETRIZATION"];
  std::cout<<"spin symmetrization is "<<spin_symmetrize<<std::endl;

  int max_it = parms["sc.MAX_IT"];
  bool adjust_mu = parms["sc.ADJUST_MU"];
  double target_density = adjust_mu?parms["sc.TARGET_DENSITY"]:0.;
  double convergence_tolerance=parms["sc.SELFCONSISTENCY_CONVERGENCE"];
  std::string g0_output_path=parms["dca.G0_OMEGA_OUTPUT"].as<std::string>();
  bool write_iterations = parms["sc.WRITE_ITERATIONS"];

  //define all Matsubara and imaginary time Green's functions needed
  GFW G_omega=G0_omega;
  GFW G_omega_old=G0_omega;
  GFW G_omega_ksum=G0_omega;
  GFW G0_omega_new=G0_omega;
  GFW self_energy=G0_omega;

  GFT G_tau(G0_tau);

  clustertrans.initial_G0(G0_omega, G0_tau);

  //This is the basic DCA self-consistency loop
  bool converged = false;
  double change = 0.;
  for(int iteration_ctr=0;iteration_ctr < max_it; ++iteration_ctr) {
    std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;

    std::cout<<"spin symmetrization is "<<spin_symmetrize<<std::endl;
    clustertrans.symmetrize_k_space(G0_omega);
    if(spin_symmetrize){ clustertrans.spin_symmetrize(G0_omega); }

    alps::gf::fourier_frequency_to_time(G0_omega, G0_tau);
    solver.solve(G0_omega, G0_tau, G_omega, G_tau);

    if(spin_symmetrize){ clustertrans.spin_symmetrize(G_omega); }
    clustertrans.symmetrize_k_space(G_omega);

    clustertrans.dyson_get_self_energy(self_energy, G_omega, G0_omega);
    if(!adjust_mu)
      clustertrans.matsubara_k_sum(G_omega_ksum, self_energy);
    else{
      double mu=clustertrans.matsubara_k_sum_adjust_mu(G_omega_ksum, self_energy, target_density);
      solver.adjust_mu(mu);
    }
    clustertrans.symmetrize_k_space(G_omega_ksum);
    if(spin_symmetrize){ clustertrans.spin_symmetrize(G_omega_ksum); }

    clustertrans.dyson_get_green0(G0_omega_new, self_energy, G_omega_ksum);

    if(write_iterations)
      print_all_green_functions(iteration_ctr, G0_omega, G_omega, self_energy, G0_tau, G_tau, beta);

    change = (G_omega-G_omega_old).norm();
    std::cout<<"iteration "<<iteration_ctr<<": change of G "<<change<<std::endl;
    if (change <= convergence_tolerance) { converged = true; break; }

    G_omega_old=G_omega;
    G0_omega=G0_omega_new;
    write_new_G0omega(G0_omega_new, g0_output_path);
  }
  // Running a fixed number of iterations is the normal way to use dca;
  // convergence is judged afterwards from the iterations written.
  if(converged)
    std::cout<<"dca stopped: the change of G, "<<change<<", is below sc.SELFCONSISTENCY_CONVERGENCE="
             <<convergence_tolerance<<"."<<std::endl;
  else if(max_it>0)
    std::cout<<"dca stopped after sc.MAX_IT="<<max_it<<" iterations; the last change of G was "<<change
             <<" (sc.SELFCONSISTENCY_CONVERGENCE="<<convergence_tolerance<<")."<<std::endl;
}
// Runs the loop in the framework of TRANSFORMER, which provides the Green's
// function types of its framework and the starting G0.
template<class TRANSFORMER> void run_dca(const std::string &param_file_name, const alps::params &parms){
  TRANSFORMER T(parms);
  ImpuritySolver solver(parms["solver.SOLVER"].as<std::string>(), param_file_name, parms);
  typename TRANSFORMER::matsubara_gf G0_omega = T.make_matsubara_gf();
  typename TRANSFORMER::itime_gf G0_tau = T.make_itime_gf(parms["N"]);
  selfconsistency_loop_DCA(parms, solver, T, G0_omega, G0_tau);
}

void perform_cluster_loop(const std::string &param_file_name, const alps::params &parms){
  const bool nambu=parms["sc.NAMBU"];
  const bool afm=parms["dca.AFM_SELFCONSISTENCY"];
  const bool charge_order=parms["dca.CO_SELFCONSISTENCY"];
  if(int(nambu)+int(afm)+int(charge_order)>1)
    throw std::invalid_argument("choose at most one of sc.NAMBU, dca.AFM_SELFCONSISTENCY and dca.CO_SELFCONSISTENCY");

  if(nambu){
    // In Nambu space, spin symmetrization imposes the particle-hole relations
    // of a G without Zeeman field (Re G_dn,dn = -Re G_up,up, a real F). With
    // a field it would remove the field's effect from G0 in every iteration,
    // while the k-sum applies it.
    const bool spin_symmetrize = parms["sc.SPIN_SYMMETRIZATION"];
    const double h = parms["dca.H"], h_init = parms["dca.H_INIT"];
    if(spin_symmetrize && (h!=0. || h_init!=0.))
      throw std::invalid_argument("sc.NAMBU with a Zeeman field (dca.H or dca.H_INIT) needs sc.SPIN_SYMMETRIZATION = false:"
                                  " the spin symmetrization assumes no field");
    // The pairing field is (dca.ETA1 + i dca.ETA2) times the d-wave form factor
    // cos kx - cos ky. Its phase is a global U(1) gauge choice: a complex field
    // describes the same physics as the real field |eta|. G0(tau) and its tails
    // are real, so a complex field would hand the solver a wrong G0(tau).
    const double eta2 = parms["dca.ETA2"], eta2_init = parms["dca.ETA2_INIT"];
    if(eta2!=0. || eta2_init!=0.)
      throw std::invalid_argument("sc.NAMBU: a complex pairing field dca.ETA1 + i dca.ETA2 with the d-wave form factor"
                                  " differs from the real field |eta| only by a global phase, and G0(tau) is real;"
                                  " set dca.ETA2 (and dca.ETA2_INIT) to 0 and dca.ETA1 to |eta|.");
    std::cout<<"running a DCA calculation allowing for superconducting order"<<std::endl;
    run_dca<NambuClusterTransformer>(param_file_name, parms);
  }else if(afm || charge_order){
    std::cout<<"running a DCA calculation in the double unit cell framework"<<std::endl;
    run_dca<DoubleCellClusterTransformer>(param_file_name, parms);
  }else{
    std::cout<<"running a DCA calculation in the normal state"<<std::endl;
    run_dca<NormalStateClusterTransformer>(param_file_name, parms);
  }
}
