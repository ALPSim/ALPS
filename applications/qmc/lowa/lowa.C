/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Si_mulations
*
* ALPS Applications
*
* Copyright (C) 2006-2010 by Lode Pollet <lpollet@physics.harvard.edu>,
*                            Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id: Lowa.C 3520 2010-03-03 16:55:00Z tamama $ */


/*
 *
 * 1) Code modification      -- partly done
 * 2) Replacing raw pointers -- not done yet
 *
 */


#include "./lowa.h"
#include "./include/lowa.insertion_deletion.h"
#include "./include/lowa.propagation.h"
#include "./include/lowa.measurement.h"
#include "./include/lowa.dostep.h"
#include "./include/lowa.io.h"
#include <alps/alea/detailedbinning.h>


#include <iostream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <fstream>


namespace alps {
namespace applications {


void Lowa::print_copyright(std::ostream& out)
{
  out << "**************************************************************************************************\n"
      << "Locally-Optimized Worm Algorithm (LOWA)  -- version 2.0\n"
      << "  Continuous Time Worldline Path Integral Quantum Monte Carlo Si_mulation\n"
      << "\n"
      << "  copyright(c) 2006-2010 by Lode Pollet     <lpollet@physics.harvard.edu>     (original LOWA code),\n"
      << "                            Ping Nang Ma    <pingnang@itp.phys.ethz.ch>       (extension into ALPS),\n"
      << "                            Matthias Troyer <troyer@itp.phys.ethz.ch>         (ALPS library)       \n"
      << "**************************************************************************************************\n\n";
}


Lowa::Lowa(const alps::ProcessList& where,const alps::Parameters& p,int node)
  : alps::scheduler::MCRun(where,p,node)

  , label(static_cast<uint32_t>(p["LABEL"]))
  , Nmeasure(static_cast<uint64_t>(p["N_MEASURE"]))
  , Nmeasure_green(static_cast<uint64_t>(p["N_MEASURE_GREEN"]))
  , Ntest(static_cast<uint64_t>(p["N_TEST"]))
  , Nsave(static_cast<uint64_t>(p["N_SAVE"]))
  , sweeps(0)
  , thermalization_sweeps(static_cast<uint64_t>(p["THERMALIZATION"]))
  , total_sweeps(static_cast<uint64_t>(p["SWEEPS"]))

  , dim(static_cast<int>(p["DIM"]))
  , beta(static_cast<parm_type>(p["BETA"]))
  , _mu(static_cast<parm_type>(p["CP"]))
  , _as(static_cast<parm_type>(p["SCATTERING_LENGTH"]))
  , mass(static_cast<parm_type>(p["MASS"]))
  , _tof(static_cast<parm_type>(p["TIME_OF_FLIGHT"]))
  , Ncan(static_cast<site_type>(p["CANONICAL_VALUE"]))
  , _statemax(static_cast<int>(p["MAXIMUM_NUMBER_PER_SITE"]))
  , _Eoffset(static_cast<parm_type>(p["E_OFF"]))
  , is_doublon_treated_as_hole(static_cast<bool>(p["IS_DOUBLON_OR_HIGHER_TREATED_AS_HOLE"]))
  , _binsize(static_cast<time_type>(p["BIN_SIZE"]))
{
  // Further initialization
  zcoord = 2*dim;

  if (dim == 1) {
    Ls     = new site_type;
    v0     = new parm_type;
    vc     = new parm_type;
    lambda = new parm_type;
    waist  = new parm_type;
    phase  = new parm_type;
  }
  else {
    Ls     = new site_type [dim];
    v0     = new parm_type [dim];
    vc     = new parm_type [dim];
    lambda = new parm_type [dim];
    waist  = new parm_type [dim];
    phase  = new parm_type [dim];
  }

  for (int k = 0; k < dim; ++k)    {  Ls[k] = static_cast<site_type>(p["LENGTH"]); }
  _N = 1;    for (int k=0; k < dim; ++k)      {  _N   *= Ls[k];  }
  _Nxy = 1;  for (int k=0; k < (dim-1); ++k)  {  _Nxy *= Ls[k];  }

  if (dim == 3) {
    v0[0] = static_cast<parm_type>(p["V0_x"]);   v0[1] = static_cast<parm_type>(p["V0_y"]);   v0[2] = static_cast<parm_type>(p["V0_z"]);
    vc[0] = static_cast<parm_type>(p["VC_x"]);   vc[1] = static_cast<parm_type>(p["VC_y"]);   vc[2] = static_cast<parm_type>(p["VC_z"]);
    lambda[0] = static_cast<parm_type>(p["WAVELENGTH_x"]);   lambda[1] = static_cast<parm_type>(p["WAVELENGTH_y"]);   lambda[2] = static_cast<parm_type>(p["WAVELENGTH_z"]);
    waist[0]  = static_cast<parm_type>(p["WAIST_x"]);        waist[1]  = static_cast<parm_type>(p["WAIST_y"]);        waist[2] = static_cast<parm_type>(p["WAIST_z"]);
    phase[0] = static_cast<parm_type>(p["MASS"]) * amu * lambda[0] * lambda[0] / (8*_tof*hbar) * 1e-8;
    phase[1] = static_cast<parm_type>(p["MASS"]) * amu * lambda[1] * lambda[1] / (8*_tof*hbar) * 1e-8;
    phase[2] = static_cast<parm_type>(p["MASS"]) * amu * lambda[2] * lambda[2] / (8*_tof*hbar) * 1e-8;
  }
  else if (dim == 2) {  
    // Tama : to be added...
  }
  else if (dim == 1) {
    // Tama : to be added...
  }

  // lattice initialization
  _lattice.init(Ls[0],Ls[1],Ls[2]);

  // band structure calculations (BHM)
  std::cout << "\nCalculating band parameters for boson Hubbard model...\t";
  _band_structure.init(label,Ls[0],v0[0],v0[1],v0[2],_as,mass,lambda[0],lambda[1],lambda[2],waist[0],waist[1],waist[2],vc[0],vc[1],vc[2]);
  _site0 = _band_structure.site0();
  std::cout << "\t ...done.";

  //initialise parameters 
  init();
  print_params(std::cout);
}


void Lowa::init()
{
// ### INTRINSIC SIMULATION VARIABLES

  counter_MEASURE = 0;
  counter_MEASURE_GREEN = 0;
  counter_TEST    = 0;
  MCstep_total = 0.;
  MCold = MCstep + MCstep_total;

  time (&times1);


// ### PRINT and I/O (ASCII/hdf5)
  // definitions of filenames
  filename1       = obtain_filename("qmc_saveconf");
  filename1_trial = obtain_filename("tqmc_saveconf");

  filename_N           = obtain_filename("N");
  filename_mN          = obtain_filename("mN");
  filename_dns0        = obtain_filename("dnscenter");
  filename_mdns0       = obtain_filename("mdnscenter");
  filename_proj_cymdns = obtain_filename("proj_cymdns");
  filename_cs_cymdns   = obtain_filename("cs_cymdns");

  filename_mdns        = obtain_filename_hdf5("mdns");
  filename_mdnsmat     = obtain_filename("mdnsmat");
  filename_mdnsmatinf  = obtain_filename("mdnsmatinf");


// ### WOLDLINE DESCRIPTION

  _kinks_on_lattice = new kinks_type [_N];
  _mu_eff            = new parm_type [_N];

  if (dim==3)
  {
    // setup chemical potential
#ifdef TRAPPEDSYSTEM
    for (site_type p = 0; p < _N; ++p)  {  _mu_eff[p] = (V(p) - _mu);  }
#else
    for (site_type p = 0; p < _N; ++p)  {  _mu_eff[p] = -_mu;          }
#endif
  }
  else if (dim==2)
  {
    // setup chemical potential
#ifdef TRAPPEDSYSTEM
    // *** Tama : this has to be changed...
    double my = Ls[1]/2.0 - 0.5;
    double mx = Ls[0]/2.0 - 0.5 ;
    for (site_type j =0; j < Ls[1]; j++) {
      site_type y = j*Ls[0];
      for (site_type i =0; i < Ls[0]; i++) {
        site_type p = i + y;
        double d = vc[1]*(j-my)*(j-my) + vc[0]*(i-mx)*(i-mx);
        if ( (std::abs(vc[0]) > tol)) {
          if ((j==0) || (j==Ls[1]-1) || (i==0) || (i==Ls[0]-1) )  { d = 1000000.*std::abs(d); }
        }
        d -= _mu;
        _mu_eff[p] = d;
        //std::cout << "\n chempot : " << i << "\t" << j <<  "\t" << d;
      }
    }
#else
    for (site_type p = 0; p < _N; ++p)  {  _mu_eff[p] = -_mu;          }
#endif
  }
  else if (dim==1)
  {
    // setup chemical potential
#ifdef TRAPPEDSYSTEM
    // *** Tama : this has to be changed
    double mx = Ls[0]/2.0 - 0.5;
    for (site_type i =0; i < Ls[0]; i++) {
      double d = vc[0]*(i-mx)*(i-mx);
      d -= _mu ;
      if ( (std::abs(vc[0]) > tol) && (i==0)) {
        d = 1000000. * std::abs(d);
      }
      else if ( (std::abs(vc[0]) > tol) && (i==Ls[0]-1)) {
        d = 1000000. * std::abs(d);
      }
      _mu_eff[i] = d;
      //std::cout << "\n chempot : " << i << "\t" << d;
    }
#else
    for (site_type i =0; i < Ls[0]; ++i)  {  _mu_eff[i] = -_mu;  }
#endif
  }


// ### MEASUREMENT

  new_measurement = true;
  _nrvertex       = 0;

  _state        = new fock_basis_type [_N];

  av_dnsmat     = new obs_type [_N];
  av_dnsmat_inf = new obs_type [_N];
  reset_av_dnsmat();

  _nrbin        = std::ceil(std::sqrt((Ls[0]-1-mid(0))*(Ls[0]-1-mid(0)) + (Ls[1]-1-mid(1))*(Ls[1]-1-mid(1)))/_binsize);

  _proj_binnr   = new site_type [_Nxy];
  _proj_binfreq = new site_type [_nrbin];

  for (site_type index=0; index < _nrbin; ++index)  {  _proj_binfreq[index] = 0;  }
  for (site_type p=0, j=0; j < Ls[1]; ++j) {
    for (site_type i=0; i < Ls[0]; ++i, ++p) {
      site_type cur_binnr = std::floor(std::sqrt((i-mid(0))*(i-mid(0)) + (j-mid(1))*(j-mid(1)))/_binsize);
      _proj_binnr[p]            = cur_binnr;
      _proj_binfreq[cur_binnr] += 1;
    }
  }

  _proj_binstate.resize(_nrbin);
  _cs_binstate.resize(_nrbin);

  //measurements << alps::RealObservable("Total Particle Number (Actual)");
  //measurements << alps::RealObservable("Total Particle Number (Measured)");
  //measurements << alps::RealObservable("Density at center (Actual)");
  //measurements << alps::RealObservable("Density at center (Measured)");
  //measurements << alps::RealObservable("Kinetic Energy");
  //measurements << alps::RealObservable("Potential Energy");
  //measurements << alps::RealObservable("Energy");
  //measurements << alps::RealVectorObservable("Columnn integrated Density (cylindrically binned)");
  //measurements << alps::RealVectorObservable("Cross sectional Density (cylindrically binned)");


// ### ABOUT WORM DESCRIPTION

  site_it       = new kinks_iterator_type [_N];
  dummy_kink_it = new kinks_iterator_type [_N];
  for (site_type i=0 ; i < _N; ++i)  {  site_it[i] = _kinks_on_lattice[i].begin(); }

  is_worm_diagonal = true;

  worm_head.init(zcoord);
  worm_tail.init(zcoord);


// ### ABOUT WORM INSERTION, PROPAGATION, DELETION PROCESSING

  nbs                                  = new fock_basis_type [zcoord];
  local_transition_weight              = new obs_type [zcoord + 1];
  cummulative_local_transition_weight  = new obs_type [zcoord + 1];

  no_of_accepted_worm_insertions = 0.;
  no_of_proposed_worm_insertions = 0.;

  MCstep = 0.;


// final stage of initialization

  load_lowa();  // this will check automatically if we can load from file...

  std::cout << "\n# Finished intializing.\n";
}


bool Lowa::test_conf()
{ 
  if ((!is_worm_diagonal) && (worm_it != site_it[worm_head.to()])) {
    std::cout << "\n Error with site_it and worm_it";
    return(0);
  }
  for (site_type i = 0; i < _N; i++) {
    int n = _kinks_on_lattice[i].size();
    int ii = 0;
    for (kinks_iterator_type it = _kinks_on_lattice[i].begin(); it != _kinks_on_lattice[i].end(); ++it, ++ii) {
      if (it->before() < 0 || it->before() > _statemax) {
	std::cout << "\n Error with density  << " << i << "\t "<< &*it << "\t" << it->time();
	return (0);
      }
      if (it->after() < 0 || it->after() > _statemax) {
	std::cout << "\n Error with density  << " << i << "\t "<< &*it << "\t" << it->time();
	return (0);
      }
    }
  }
  return (1);
}


bool Lowa::change_parameter(const std::string& name,const alps::StringValue& value)
{
  if(name=="SWEEPS")
    total_sweeps=static_cast<uint32_t>(value);
  else if (name=="THERMALIZATION" && !is_thermalized())
    thermalization_sweeps=static_cast<uint32_t>(value);
  else
    return false; // cannot do it
  return true; // can do it
}



} // ending namespace applications
} // ending namespace alps

