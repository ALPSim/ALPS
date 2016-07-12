/*****************************************************************************
 *
 * Copyright (C) 2005-2015 Emanuel Gull <egull@umich.edu>
 *
 *****************************************************************************/

#include <fstream>
#include <iomanip>
#include "impurity_solver.h"
#include "selfconsistency.h"
#include "hilbert_transform.h"
#include "alps/gf/tail.hpp"
#include "alps/gf/fourier.hpp"

void define_selfconsistency_parameters(alps::params &p){

  FrequencySpaceHilbertTransformer::define_parameters(p);

  p.define<int>("N", "Number of imaginary time discretization points");
  p.define<int>("NMATSUBARA", "Number of matsubara frequency discretization points");
  p.define<int>("FLAVORS", 2, "Number of spins or (diagonal) orbitals");
  p.define<double>("BETA", "Inverse temperature");
  p.define<double>("H", "Staggered (antiferromagnetic) field");

  p.define<double>("sc.SELFCONSISTENCY_CONVERGENCE", 1.e-10, "Difference in Green's function at which the iteration is converged");
  p.define<int>("sc.MAX_IT", 1000, "Maximum number of DMFT iterations");
  p.define<double>("sc.RELAX_RATE", 1., "Convergence acceleration (or deceleration) parameter. Choose 1 if unsure.");
  p.define<bool>("sc.AFM_SELFCONSISTENCY", false, "Use the AFM self-consistency condition.");
  p.define<bool>("sc.SPIN_SYMMETRIZATION", true, "symmetrize 'up' and 'down' spins after each solver call.");
  p.define<std::string>("G0_OMEGA_OUTPUT", "Output filename for G0");
}
void run_convergence_acceleration(single_site_matsubara_gf &G_new, const single_site_matsubara_gf &G_old, double relax_rate){
  if(relax_rate !=1){
    std::cout<<"using over/underrelaxation with rate: "<<relax_rate<<std::endl;
    for(alps::gf::matsubara_index i(0); i<G_new.mesh1().extent(); ++i) {
      for(alps::gf::index flavor(0);flavor<G_new.mesh2().extent(); ++flavor) {
        G_new(i,flavor)=relax_rate*G_new(i,flavor)+(1-relax_rate)*G_old(i,flavor);
      }
    }
  }
}
bool converged(const single_site_matsubara_gf &G_omega, const single_site_matsubara_gf &G_omega_old, double tolerance){
  double max_diff=0;
  return ((G_omega-G_omega_old).norm() > max_diff)?false:true;
}
template<class GFW> void write_new_G0omega(const GFW &G0_omega, const std::string &path){
  alps::hdf5::archive ar(path, alps::hdf5::archive::WRITE);
  G0_omega.save(ar, "/G0_omega");
}
void print_all_green_functions(int iteration_ctr,
    const single_site_matsubara_gf &G0_omega,
    const single_site_matsubara_gf &G_omega,
    const single_site_matsubara_gf &selfenergy,
    const single_site_itime_gf &G0_tau,
    const single_site_itime_gf &G_tau, double beta){

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

  G0_omega_file<<std::setprecision(14)<<G0_omega;
  G0_tau_file<<std::setprecision(14)<<G0_tau;
  G_omega_file<<std::setprecision(14)<<G_omega;
  G_tau_file<<std::setprecision(14)<<G_tau;
  G_omega_file_2<<std::setprecision(14)<<G_omega;
  G_tau_file_2<<std::setprecision(14)<<G_tau;
  selfenergy_file<<std::setprecision(14)<<selfenergy;
}
void compute_selfenergy(const single_site_matsubara_gf &G0_omega,
    const single_site_matsubara_gf &G_omega,
    single_site_matsubara_gf &selfenergy){


  for(alps::gf::matsubara_index i(0);i<G0_omega.mesh1().extent();++i){
    for(alps::gf::index s(0);s<G0_omega.mesh2().extent();++s){
      selfenergy(i,s)=1./G0_omega(i,s)-1./G_omega(i,s);
    }
  }
  std::cout<<"computing selfenergy from G0: "<<G0_omega(alps::gf::matsubara_index(0), alps::gf::index(0))<<" and G: "<<G_omega(alps::gf::matsubara_index(0), alps::gf::index(0))<<" to obtain "<<selfenergy(alps::gf::matsubara_index(0), alps::gf::index(0))<<std::endl;
}

void selfconsistency_loop_single_site(const alps::params& parms, ImpuritySolver& solver, FrequencySpaceHilbertTransformer& hilbert)
{
  int n_tau=parms["N"];
  int n_matsubara=parms["NMATSUBARA"];
  int n_orbital=parms["FLAVORS"];

  double beta = parms["BETA"];
  double h = parms["H"];
  bool spin_symmetrize = parms["sc.SPIN_SYMMETRIZATION"];
  double relax_rate=parms["sc.RELAX_RATE"];
  int max_it=parms["sc.MAX_IT"];
  double convergence_tolerance=parms["sc.SELFCONSISTENCY_CONVERGENCE"];
  std::string g0_output_path = parms["G0_OMEGA_OUTPUT"];


  single_site_matsubara_gf G0_omega = hilbert.initial_G0();
  single_site_matsubara_gf G_omega = G0_omega;
  single_site_matsubara_gf selfenergy = G0_omega;

  int iteration_ctr = 0;
  //define the Green's functions
  single_site_matsubara_gf G_old_omega(G0_omega);
  single_site_matsubara_gf G0_old_omega(G0_omega);
  single_site_itime_gf G_tau(alps::gf::itime_sigma_gf(alps::gf::itime_mesh(beta,n_tau),alps::gf::index_mesh(n_orbital)));
  single_site_itime_gf G0_tau(alps::gf::itime_sigma_gf(alps::gf::itime_mesh(beta,n_tau),alps::gf::index_mesh(n_orbital)));
  single_site_itime_gf G_old_tau(alps::gf::itime_sigma_gf(alps::gf::itime_mesh(beta,n_tau),alps::gf::index_mesh(n_orbital)));
  single_site_itime_gf G0_old_tau(alps::gf::itime_sigma_gf(alps::gf::itime_mesh(beta,n_tau),alps::gf::index_mesh(n_orbital)));

  do {
    iteration_ctr++;
    std::cout<<"starting iteration nr. "<<iteration_ctr<<std::endl;
    G_old_omega = G_omega;
    G0_old_omega = G0_omega;
    G_old_tau = G_tau;
    G0_old_tau = G0_tau;
    alps::gf::fourier_frequency_to_time(G0_omega,G0_tau);
    std::cout<<"running solver."<<std::endl<<std::flush;
    solver.solve(G0_omega,G0_tau, G_omega, G_tau);
    if(spin_symmetrize){
      hilbert.spin_symmetrize(G_tau);
      hilbert.spin_symmetrize(G_omega);
    }

    compute_selfenergy(G0_omega, G_omega, selfenergy);

    G0_omega = hilbert(G_omega, G0_omega);

    //relaxation to speed up /slow down convergence
    run_convergence_acceleration(G0_omega, G0_old_omega, relax_rate);
    write_new_G0omega(G0_omega, g0_output_path);
    print_all_green_functions(iteration_ctr, G0_omega, G_omega, selfenergy, G0_tau, G_tau, beta);
  }while (iteration_ctr < max_it && !converged(G_omega,G_old_omega,convergence_tolerance));
}



//construct the Hilbert transformer and run a single site loop
void perform_single_site_loop(const std::string &param_file_name, const alps::params &parms){

  FrequencySpaceHilbertTransformer *transform_ptr; //super ugly
  if(parms.exists("sc.TWODBS"))
    transform_ptr=setup_2D_hilbert_transform(parms);
  else if (parms.exists("sc.DOSFILE"))
    transform_ptr=setup_DOS_hilbert_transform(parms);
  else
    transform_ptr=setup_Bethe_hilbert_transform(parms);

  std::string solver_path=parms["solver.SOLVER"];
  ImpuritySolver solver(solver_path,param_file_name, parms);
  selfconsistency_loop_single_site(parms, solver, *transform_ptr);
  delete transform_ptr;
}

