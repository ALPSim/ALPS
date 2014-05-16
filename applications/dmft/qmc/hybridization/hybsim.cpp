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
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/
#ifndef HYB_SIM_MAIN
#define HYB_SIM_MAIN
#endif

#include <iomanip>
#include"hyb.hpp"

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
  output_period=parms["OUTPUT_PERIOD"]|100000;
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
//  U_ = parms["U"];
//  MU_ = parms["MU"];
  
//  if ((parms["USE_FRACTION"]| 0)) {
//    fraction = beta/((U_>MU_)?MU_:U_);
////    std::cerr << "Using only " << fraction << " of segments" << std::endl;
//  }
//  else fraction = beta;
  
  //initializing updates parameters
  N_meas = parms["N_MEAS"];                                                        //number of updates per measurement
  N_hist_orders = parms["N_HISTOGRAM_ORDERS"]|50;                                  //number of orders that are measured for the order histogram
  MEASURE_timeseries = parms["TIMESERIES"]|0;
  NUM_BINS = parms["NUM_BINS"]|0;
  //std::cerr << "NUM_BINS = " << NUM_BINS << std::endl;
  //initializing measurement parameters
  spin_flip = parms["SPINFLIP"]| 0;                                                //whether to perform local spin-flip updates
  global_flip = parms["GLOBALFLIP"]| 0;                                                //whether to perform global spin-flip updates
  MEASURE_nnt = parms["MEASURE_nnt"]| 0;                                           //measure density-density correlation function in imaginary time
  MEASURE_nnw = parms["MEASURE_nnw"]| 0;                                           //measure density-density correlation function in frequency
  MEASURE_nn = parms["MEASURE_nn"]|0;                                              //measure density-density correlation function at equal times
  MEASURE_g2w = parms["MEASURE_g2w"]|0;                                            //measure two-particle Green function
  MEASURE_h2w = parms["MEASURE_h2w"]|0;                                            //measure higher-order two-particle correlator
  MEASURE_time = parms["MEASURE_time"]|1;                                          //measure in imaginary time (ON by default)
  MEASURE_freq = parms["MEASURE_freq"]|0;                                          //measure in frequency space
  MEASURE_legendre = parms["MEASURE_legendre"]|0;                                  //measure in legendre polynomials
  MEASURE_sector_statistics = parms["MEASURE_sector_statistics"]|0;                //measure sector statistics
  N_w = parms["N_MATSUBARA"]|0;                                                    //number of Matsubara frequencies for gw
  N_l = parms["N_LEGENDRE"]|0;                                                     //number of Legendre polynomial coefficients
  N_t = parms["N_TAU"];                                                            //number of tau slices for gt
  N_nn = parms["N_nn"]|0;                                                          //number of tau-points on which density density correlator is measured
  N_w2 = parms["N_w2"]|0;                                                          //number of Matsubara frequency points for two-particle measurements
  N_W = parms["N_W"]|0;                                                            //number of bosonic Matsubara frequency points for two-particle measurements
  N_w_aux = (N_w2+N_W>1 ? N_w2+N_W-1 : 0);                                         //number of Matsubara frequency points for the measurment of M(w1,w2)
  
  //create measurement objects
  create_measurements();
  
  if(crank==0){
    std::cout<<"Hybridization Expansion Simulation CT-HYB"<<std::endl;
    std::cout<<"Part of the ALPS DMFT Project"<<std::endl;
    std::cout<<"Usage requires citation of the ALPS CT-HYB paper and the ALPS paper"<<std::endl;
    std::cout<<"Refer to the documentation for more information."<<std::endl;
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
  if(parms["MEASURE_freq"]|false && !parms.defined("N_MATSUBARA")) throw std::invalid_argument("please specify parameter N_MATSUBARA for # of Matsubara frequencies to be measured");
  
  if(parms["MEASURE_legendre"]|false && !parms.defined("N_LEGENDRE")) throw std::invalid_argument("please specify parameter N_LEGENDRE for # of Legendre coefficients to be measured");
  if(parms["MEASURE_legendre"]|false && !parms.defined("N_MATSUBARA")) throw std::invalid_argument("please specify parameter N_MATSUBARA for # of Matsubara frequencies");
  if(parms["MEASURE_nnt"]|false && !parms.defined("N_nn")) throw std::invalid_argument("please specify the parameter N_nn for # of imaginary time points for the density-density correlator");
  if(parms["MEASURE_nnw"]|false && !parms.defined("N_W")) throw std::invalid_argument("please specify the parameter N_W for # of bosonic frequencies for the density-density correlator");
  if(parms["MEASURE_g2w"]|false || parms["MEASURE_h2w"]|false ){
    if(!parms.defined("N_w2") ) throw std::invalid_argument("please specify the parameter N_w2 for # of fermionic Matsubara frequencies for two-particle functions");
    if(!parms.defined("N_W") ) throw std::invalid_argument("please specify the parameter N_W for # of bosonic Matsubara frequencies for two-particle functions");
    if((int)parms["N_w2"]%2!=0) throw std::invalid_argument("parameter N_w2 must be even");
  }
  if(parms["COMPUTE_VERTEX"]|false){
    if( !(parms["MEASURE_freq"]|false) ) throw std::invalid_argument("frequency measurement is required for computing the vertex, please set MEASURE_freq=1");
    
    if(! (parms["MEASURE_g2w"]|false || parms["MEASURE_h2w"]|false ) ) throw std::invalid_argument("at least one two-particle quantity is required for computing the vertex, set MEASURE_g2w=1 or MEASURE_h2w=1");
    if((int) parms["N_MATSUBARA"] < ((int)parms["N_w2"]/2 + (int)parms["N_W"] - 1) ) throw std::invalid_argument("for computing the vertex, N_MATSUBARA must be at least N_w2/2+N_W-1");
  }
  VERBOSE = (parms["VERBOSE"]|false);
  
  return;
}


void hybridization::show_info(const alps::params &parms, int crank){
  if(!(parms["VERBOSE"]|false)) return;

  //provide info on what is measured and how long the simulation will run
  if(!crank){
    if(parms["MEASURE_time"]|true) std::cout << "measuring gt" << std::endl;
    if(parms["MEASURE_freq"]|false) std::cout << "measuring gw" << std::endl << "measuring fw" << std::endl;
    if(parms["MEASURE_legendre"]|false) std::cout << "measuring gl" << std::endl << "measuring fl" << std::endl;
    if(parms["MEASURE_g2w"]|false) std::cout << "measuring g2w" << std::endl;
    if(parms["MEASURE_h2w"]|false) std::cout << "measuring h2w" << std::endl;
    if(parms["MEASURE_nn"]|false) std::cout << "measuring nn" << std::endl;
    if(parms["MEASURE_nnt"]|false) std::cout << "measuring nnt" << std::endl;
    if(parms["MEASURE_nnw"]|false) std::cout << "measuring nnw" << std::endl;
    if(parms["MEASURE_sector_statistics"]|false) std::cout << "measuring sector statistics" << std::endl;
    if(parms["COMPUTE_VERTEX"]|false) std::cout << "vertex will be computed" << std::endl;
    if(parms.defined("RET_INT_K")) std::cout << "using retarded interaction" << std::endl;
    if(parms.defined("U_MATRIX")) std::cout << "reading U matrix from file " << parms["U_MATRIX"] << std::endl;
    if(parms.defined("MU_VECTOR")) std::cout << "reading MU vector from file " << parms["MU_VECTOR"] << std::endl;
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
