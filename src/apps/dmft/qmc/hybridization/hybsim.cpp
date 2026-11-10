/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <egull@umich.edu>,
 *                       Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/
#ifndef HYB_SIM_MAIN
#define HYB_SIM_MAIN
#endif

#include <iomanip>
#include"hyb.hpp"

void hybridization::define_parameters(alps::params &parms){
  alps::mcbase::define_parameters(parms); //SEED
  alps::input::define_hybridization_parameters(parms);
}

hybridization::hybridization(const alps::params &parms, int crank_)
: hybridization(parms, crank_, alps::input::read_hybridization_parameters(parms))
{
}

hybridization::hybridization(const alps::params &parms, int crank_, const alps::input::hybridization_parameters &input)
: alps::mcbase(parms, crank_),
crank(crank_),
local_config(input,crank),
hyb_config(input)
{
  VERBOSE = input.verbose;
  show_info(input,crank);
  
  //initializing general simulation constants
  nacc.assign(7,0.);
  nprop.assign(7,0.);
  sweep_count=0;
  output_period=input.output_period;
  
  update_type.clear();
  update_type.push_back("change zero state   ");
  update_type.push_back("insert segment      ");
  update_type.push_back("remove segment      ");
  update_type.push_back("insert anti-segment ");
  update_type.push_back("remove anti-segment ");
  update_type.push_back("swap segment        ");
  update_type.push_back("global flip         ");
  
  sweeps=0;                                                                        //Sweeps currently done
  thermalization_sweeps = input.thermalization;                                 //Sweeps to be done for thermalization
  total_sweeps = input.sweeps;                                                  //Sweeps to be done in total
  n_orbitals = input.n_orbitals;                                                //number of orbitals
  sign = 1.;                                                                       //fermionic sign. plus or minus one.
  
  //initializing physics parameters
  beta = input.beta;                                                            //inverse temperature
  
  //initializing updates parameters
  N_meas = input.n_meas;                                                        //number of updates per measurement
  N_hist_orders = input.n_histogram_orders;                                  //number of orders that are measured for the order histogram
  MEASURE_timeseries = static_cast<bool>(input.timeseries);
  NUM_BINS = input.num_bins;
  //std::cerr << "NUM_BINS = " << NUM_BINS << std::endl;
  //initializing measurement parameters
  spin_flip = input.spinflip;                                                //whether to perform local spin-flip updates
  global_flip = input.globalflip;                                                //whether to perform global spin-flip updates
  MEASURE_nnt = input.measure_nnt;                                           //measure density-density correlation function in imaginary time
  MEASURE_nnw = input.measure_nnw;                                           //measure density-density correlation function in frequency
  MEASURE_nn = input.measure_nn;                                              //measure density-density correlation function at equal times
  MEASURE_g2w = input.measure_g2w;                                            //measure two-particle Green function
  MEASURE_h2w = input.measure_h2w;                                            //measure higher-order two-particle correlator
  MEASURE_time = input.measure_time;                                          //measure in imaginary time (ON by default)
  MEASURE_freq = input.measure_freq;                                          //measure in frequency space
  MEASURE_legendre = input.measure_legendre;                                  //measure in legendre polynomials
  MEASURE_sector_statistics = input.measure_sector_statistics;                //measure sector statistics
  N_w = input.n_matsubara;                                                    //number of Matsubara frequencies for gw
  N_l = input.n_legendre;                                                     //number of Legendre polynomial coefficients
  N_t = input.n_tau;                                                            //number of tau slices for gt
  N_nn = input.n_nn;                                                          //number of tau-points on which density density correlator is measured
  N_w2 = input.n_w2;                                                          //number of Matsubara frequency points for two-particle measurements
  N_W = input.n_W;                                                            //number of bosonic Matsubara frequency points for two-particle measurements
  N_w_aux = (N_w2+N_W>1 ? N_w2+N_W-1 : 0);                                         //number of Matsubara frequency points for the measurment of M(w1,w2)
  
  //create measurement objects
  create_measurements();
  
  if(crank==0){
    std::cout<<"Hybridization Expansion Simulation CT-HYB"<<std::endl;
    std::cout<<"Part of the ALPS DMFT Project"<<std::endl;
  }
  
  start_time=clock();
  end_time=start_time+ CLOCKS_PER_SEC*((long)input.max_time);

  
  //std::cout<<"process " << crank << " starting simulation"<<std::endl;
  csize=1;
 //we don't have a nice way of getting the MPI size from ALPS, because we don't know about the communicator at this point.
 //here is a safe way of getting the pool size into csize.
#ifdef ALPS_HAVE_MPI
  int mpi_init;
  MPI_Initialized(&mpi_init);
  if(mpi_init){
     MPI_Comm_size(MPI_COMM_WORLD, &csize);
  }
#endif
  std::cout<<"process " << crank << " of total: "<<csize<<" starting simulation"<<std::endl;
}

void hybridization::show_info(const alps::input::hybridization_parameters &input, int crank){
  if(!static_cast<bool>(input.verbose)) return;

  //provide info on what is measured and how long the simulation will run
  if(!crank){
    if(static_cast<bool>(input.measure_time)) std::cout << "measuring gt" << std::endl;
    if(static_cast<bool>(input.measure_freq)) std::cout << "measuring gw" << std::endl << "measuring fw" << std::endl;
    if(static_cast<bool>(input.measure_legendre)) std::cout << "measuring gl" << std::endl << "measuring fl" << std::endl;
    if(static_cast<bool>(input.measure_g2w)) std::cout << "measuring g2w" << std::endl;
    if(static_cast<bool>(input.measure_h2w)) std::cout << "measuring h2w" << std::endl;
    if(static_cast<bool>(input.measure_nn)) std::cout << "measuring nn" << std::endl;
    if(static_cast<bool>(input.measure_nnt)) std::cout << "measuring nnt" << std::endl;
    if(static_cast<bool>(input.measure_nnw)) std::cout << "measuring nnw" << std::endl;
    if(static_cast<bool>(input.measure_sector_statistics)) std::cout << "measuring sector statistics" << std::endl;
    if(static_cast<bool>(input.compute_vertex)) std::cout << "vertex will be computed" << std::endl;
    if(!input.ret_int_k.empty()) std::cout << "using retarded interaction" << std::endl;
    if(!input.u_matrix.empty()) std::cout << "reading U matrix from file " << input.u_matrix << std::endl;
    if(!input.mu_vector.empty()) std::cout << "reading MU vector from file " << input.mu_vector << std::endl;
    std::cout << "Simulation scheduled to run " << input.max_time << " seconds" << std::endl << std::endl;
  }
  return;
}

//this is a debug function that recomputes the full weight of the local and the hybridization configuration
double hybridization::full_weight() const{
  return local_config.full_weight()*hyb_config.full_weight();
}

std::ostream &operator<<(std::ostream &os, const hybridization &hyb){
  os<<cred<<"-----------------------------------------------------------------------------------"<<cblack<<std::endl;
  os<<hyb.local_config<<std::endl;
  os<<hyb.hyb_config<<std::endl;
  os<<cred<<"-----------------------------------------------------------------------------------"<<cblack<<std::endl;
  return os;
}
std::ostream &operator<<(std::ostream &os, const segment &s){
  os<<"( "<<s.t_start_<<" , "<<s.t_end_<<" ) ";
  return os;
}
double hybridization::fraction_completed()const{
  if(!is_thermalized()) return 0.;
  double work_fraction= (sweeps-thermalization_sweeps)/(double)total_sweeps;
  double time_fraction= (clock()-start_time)/(double)(end_time-start_time);
  //return max of sweeps done and time used. Divide time used by the number of processes in pool (all work done will be added up)
  return std::max(work_fraction, time_fraction/csize);
  //return work_fraction;
}
