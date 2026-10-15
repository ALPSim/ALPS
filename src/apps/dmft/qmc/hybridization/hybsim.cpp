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

namespace {
//gives a parameter its default value if the input does not set it; the
//description documents the parameter
template<typename T> void define(alps::params &parms, std::string const &name, T const &value, std::string const &){
  if(!parms.defined(name)) parms[name]=value;
}
}

void hybridization::define_parameters(alps::params &parms){
  //simulation
  define<int>(parms, "OUTPUT_PERIOD", 100000, "number of sweeps between progress reports");
  define<int>(parms, "N_HISTOGRAM_ORDERS", 50, "number of orders measured for the order histogram");
  define<bool>(parms, "TIMESERIES", false, "use fixed binning instead of detailed binning");
  define<int>(parms, "NUM_BINS", 0, "number of bins for the measurements");
  define<bool>(parms, "VERBOSE", false, "print extra information");
  //updates
  define<bool>(parms, "SPINFLIP", false, "perform local spin-flip updates");
  define<bool>(parms, "GLOBALFLIP", false, "perform global spin-flip updates");
  //measurements
  define<bool>(parms, "MEASURE_time", true, "measure the Green's function in imaginary time");
  define<bool>(parms, "MEASURE_freq", false, "measure the Green's function on Matsubara frequencies");
  define<bool>(parms, "MEASURE_legendre", false, "measure the Green's function in Legendre polynomials");
  define<bool>(parms, "MEASURE_nn", false, "measure equal-time density-density correlations");
  define<bool>(parms, "MEASURE_nnt", false, "measure the density-density correlation function in imaginary time");
  define<bool>(parms, "MEASURE_nnw", false, "measure the density-density correlation function on Matsubara frequencies");
  define<bool>(parms, "MEASURE_g2w", false, "measure the two-particle Green's function");
  define<bool>(parms, "MEASURE_h2w", false, "measure the higher-order two-particle correlator");
  define<bool>(parms, "MEASURE_sector_statistics", false, "measure sector statistics");
  define<bool>(parms, "COMPUTE_VERTEX", false, "compute the vertex function");
  define<bool>(parms, "ACCURATE_COVARIANCE", true, "compute the covariance of the imaginary-time Green's function accurately");
  define<int>(parms, "N_MATSUBARA", 0, "number of Matsubara frequencies; required with MEASURE_freq or MEASURE_legendre");
  define<int>(parms, "N_LEGENDRE", 0, "number of Legendre coefficients; required with MEASURE_legendre");
  define<int>(parms, "N_nn", 0, "number of imaginary-time points for the density-density correlator; required with MEASURE_nnt");
  define<int>(parms, "N_W", 0, "number of bosonic Matsubara frequencies; required with MEASURE_nnw, MEASURE_g2w or MEASURE_h2w");
  define<int>(parms, "N_w2", 0, "number of fermionic Matsubara frequencies for two-particle functions; required with MEASURE_g2w or MEASURE_h2w");
  //input and output
  define<std::string>(parms, "BASEPATH", "", "path in the output file under which results are stored");
  define<bool>(parms, "TEXT_OUTPUT", false, "also write results as text files");
  define<bool>(parms, "DELTA_IN_HDF5", false, "read the hybridization function DELTA from an HDF5 file");
  define<bool>(parms, "DMFT_FRAMEWORK", false, "read HDF5 input in the layout of the DMFT framework");
  define<std::string>(parms, "U_MATRIX", "", "file with the interaction matrix; replaces U, U' and J");
  define<bool>(parms, "UMATRIX_IN_HDF5", false, "read U_MATRIX from an HDF5 file");
  define<std::string>(parms, "MU_VECTOR", "", "file with the orbital chemical potentials; replaces MU");
  define<bool>(parms, "MU_IN_HDF5", false, "read MU_VECTOR from an HDF5 file");
  define<std::string>(parms, "RET_INT_K", "", "file with the retarded interaction K(tau); enables the retarded interaction");
  define<bool>(parms, "K_IN_HDF5", false, "read RET_INT_K from an HDF5 file");
  //interaction
  define<double>(parms, "J", 0., "Hund's coupling");
}

hybridization::hybridization(const alps::params &parms, int crank_)
: alps::mcbase(parms, crank_),
crank(crank_),
local_config(parms,crank),
hyb_config(parms)
{
  sanity_check(parms); //before doing anything, check whether the input parameters make sense
  show_info(parms,crank);
  
  //initializing general simulation constants
  nacc.assign(7,0.);
  nprop.assign(7,0.);
  sweep_count=0;
  output_period=parms["OUTPUT_PERIOD"];
  //lasttime = boost::chrono::steady_clock::now();
  //delay = boost::chrono::seconds(parms["OUTPUT_PERIOD"]|600);
  
  update_type.clear();
  update_type.push_back("change zero state   ");
  update_type.push_back("insert segment      ");
  update_type.push_back("remove segment      ");
  update_type.push_back("insert anti-segment ");
  update_type.push_back("remove anti-segment ");
  update_type.push_back("swap segment        ");
  update_type.push_back("global flip         ");
  
  sweeps=0;                                                                        //Sweeps currently done
  thermalization_sweeps = parms["THERMALIZATION"];                                 //Sweeps to be done for thermalization
  total_sweeps = parms["SWEEPS"];                                                  //Sweeps to be done in total
  n_orbitals = parms["N_ORBITALS"];                                                //number of orbitals
  sign = 1.;                                                                       //fermionic sign. plus or minus one.
  
  //initializing physics parameters
  beta = parms["BETA"];                                                            //inverse temperature
  
  //initializing updates parameters
  N_meas = parms["N_MEAS"];                                                        //number of updates per measurement
  N_hist_orders = parms["N_HISTOGRAM_ORDERS"];                                  //number of orders that are measured for the order histogram
  MEASURE_timeseries = static_cast<bool>(parms["TIMESERIES"]);
  NUM_BINS = parms["NUM_BINS"];
  //std::cerr << "NUM_BINS = " << NUM_BINS << std::endl;
  //initializing measurement parameters
  spin_flip = parms["SPINFLIP"];                                                //whether to perform local spin-flip updates
  global_flip = parms["GLOBALFLIP"];                                                //whether to perform global spin-flip updates
  MEASURE_nnt = parms["MEASURE_nnt"];                                           //measure density-density correlation function in imaginary time
  MEASURE_nnw = parms["MEASURE_nnw"];                                           //measure density-density correlation function in frequency
  MEASURE_nn = parms["MEASURE_nn"];                                              //measure density-density correlation function at equal times
  MEASURE_g2w = parms["MEASURE_g2w"];                                            //measure two-particle Green function
  MEASURE_h2w = parms["MEASURE_h2w"];                                            //measure higher-order two-particle correlator
  MEASURE_time = parms["MEASURE_time"];                                          //measure in imaginary time (ON by default)
  MEASURE_freq = parms["MEASURE_freq"];                                          //measure in frequency space
  MEASURE_legendre = parms["MEASURE_legendre"];                                  //measure in legendre polynomials
  MEASURE_sector_statistics = parms["MEASURE_sector_statistics"];                //measure sector statistics
  N_w = parms["N_MATSUBARA"];                                                    //number of Matsubara frequencies for gw
  N_l = parms["N_LEGENDRE"];                                                     //number of Legendre polynomial coefficients
  N_t = parms["N_TAU"];                                                            //number of tau slices for gt
  N_nn = parms["N_nn"];                                                          //number of tau-points on which density density correlator is measured
  N_w2 = parms["N_w2"];                                                          //number of Matsubara frequency points for two-particle measurements
  N_W = parms["N_W"];                                                            //number of bosonic Matsubara frequency points for two-particle measurements
  N_w_aux = (N_w2+N_W>1 ? N_w2+N_W-1 : 0);                                         //number of Matsubara frequency points for the measurment of M(w1,w2)
  
  //create measurement objects
  create_measurements();
  
  if(crank==0){
    std::cout<<"Hybridization Expansion Simulation CT-HYB"<<std::endl;
    std::cout<<"Part of the ALPS DMFT Project"<<std::endl;
  }
  
  start_time=clock();
  end_time=start_time+ CLOCKS_PER_SEC*((long)parms["MAX_TIME"]);

  
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

void hybridization::sanity_check(const alps::params &parms){
  //check whether the input parameters make sense before computing
  //NOTE: these checks are likely not to be complete, passing all checks does not guarantee all parameters to be meaningful!
  
  //first check that all mandatory parameters are defined
  if(!parms.defined("N_TAU")) throw std::invalid_argument("please specify the parameter N_TAU");
  if(!parms.defined("BETA")) throw std::invalid_argument("please specify parameter BETA for inverse temperature");
  if(!parms.defined("N_MEAS")) throw std::invalid_argument("please specify parameter N_MEAS for measurement interval");
  if(!parms.defined("THERMALIZATION") ||
     !parms.defined("SWEEPS") ||
     !parms.defined("N_ORBITALS") ) throw std::invalid_argument("please specify parameters THERMALIZATION, SWEEPS, and N_ORBITALS");
  
  //check paramater that are conditionally required
  if(static_cast<bool>(parms["MEASURE_freq"]) && static_cast<int>(parms["N_MATSUBARA"])==0) throw std::invalid_argument("please specify parameter N_MATSUBARA for # of Matsubara frequencies to be measured");
  
  if(static_cast<bool>(parms["MEASURE_legendre"]) && static_cast<int>(parms["N_LEGENDRE"])==0) throw std::invalid_argument("please specify parameter N_LEGENDRE for # of Legendre coefficients to be measured");
  if(static_cast<bool>(parms["MEASURE_legendre"]) && static_cast<int>(parms["N_MATSUBARA"])==0) throw std::invalid_argument("please specify parameter N_MATSUBARA for # of Matsubara frequencies");
  if(static_cast<bool>(parms["MEASURE_nnt"]) && static_cast<int>(parms["N_nn"])==0) throw std::invalid_argument("please specify the parameter N_nn for # of imaginary time points for the density-density correlator");
  if(static_cast<bool>(parms["MEASURE_nnw"]) && static_cast<int>(parms["N_W"])==0) throw std::invalid_argument("please specify the parameter N_W for # of bosonic frequencies for the density-density correlator");
  if(static_cast<bool>(parms["MEASURE_g2w"]) || static_cast<bool>(parms["MEASURE_h2w"])){
    if(static_cast<int>(parms["N_w2"])==0) throw std::invalid_argument("please specify the parameter N_w2 for # of fermionic Matsubara frequencies for two-particle functions");
    if(static_cast<int>(parms["N_W"])==0) throw std::invalid_argument("please specify the parameter N_W for # of bosonic Matsubara frequencies for two-particle functions");
    if((int)parms["N_w2"]%2!=0) throw std::invalid_argument("parameter N_w2 must be even");
  }
  if(static_cast<bool>(parms["COMPUTE_VERTEX"])){
    if(!static_cast<bool>(parms["MEASURE_freq"])) throw std::invalid_argument("frequency measurement is required for computing the vertex, please set MEASURE_freq=1");
    
    if(!(static_cast<bool>(parms["MEASURE_g2w"]) || static_cast<bool>(parms["MEASURE_h2w"]))) throw std::invalid_argument("at least one two-particle quantity is required for computing the vertex, set MEASURE_g2w=1 or MEASURE_h2w=1");
    if((int) parms["N_MATSUBARA"] < ((int)parms["N_w2"]/2 + (int)parms["N_W"] - 1) ) throw std::invalid_argument("for computing the vertex, N_MATSUBARA must be at least N_w2/2+N_W-1");
  }
  VERBOSE = parms["VERBOSE"];
  
  return;
}


void hybridization::show_info(const alps::params &parms, int crank){
  if(!static_cast<bool>(parms["VERBOSE"])) return;

  //provide info on what is measured and how long the simulation will run
  if(!crank){
    if(static_cast<bool>(parms["MEASURE_time"])) std::cout << "measuring gt" << std::endl;
    if(static_cast<bool>(parms["MEASURE_freq"])) std::cout << "measuring gw" << std::endl << "measuring fw" << std::endl;
    if(static_cast<bool>(parms["MEASURE_legendre"])) std::cout << "measuring gl" << std::endl << "measuring fl" << std::endl;
    if(static_cast<bool>(parms["MEASURE_g2w"])) std::cout << "measuring g2w" << std::endl;
    if(static_cast<bool>(parms["MEASURE_h2w"])) std::cout << "measuring h2w" << std::endl;
    if(static_cast<bool>(parms["MEASURE_nn"])) std::cout << "measuring nn" << std::endl;
    if(static_cast<bool>(parms["MEASURE_nnt"])) std::cout << "measuring nnt" << std::endl;
    if(static_cast<bool>(parms["MEASURE_nnw"])) std::cout << "measuring nnw" << std::endl;
    if(static_cast<bool>(parms["MEASURE_sector_statistics"])) std::cout << "measuring sector statistics" << std::endl;
    if(static_cast<bool>(parms["COMPUTE_VERTEX"])) std::cout << "vertex will be computed" << std::endl;
    if(!parms["RET_INT_K"].cast<std::string>().empty()) std::cout << "using retarded interaction" << std::endl;
    if(!parms["U_MATRIX"].cast<std::string>().empty()) std::cout << "reading U matrix from file " << parms["U_MATRIX"] << std::endl;
    if(!parms["MU_VECTOR"].cast<std::string>().empty()) std::cout << "reading MU vector from file " << parms["MU_VECTOR"] << std::endl;
    std::cout << "Simulation scheduled to run " << parms["MAX_TIME"] << " seconds" << std::endl << std::endl;
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
