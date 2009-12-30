#include "lowa.hpp"
#include "BH_parm.hpp"
#include <alps/alea.h>
//#include <alps/alea/detailedbinning.h>
#include <cmath>
#include <iostream.h>
#include <fstream.h>


void lowa::print_copyright(std::ostream& out)
{
  out << "************************************************************************************\n"
      << "Directed-Loop Locally Optimized Worm Algorithm in Continuous Time World-Line Path\n"
      << "  Integral Quantum Monte Carlo Method using ALPS library\n"
      << "     --Version 1.5\n"
      << "  copyright(c) 2006-2008 by Lode Pollet     <pollet@itp.phys.ethz.ch>         (lowa)\n"
      << "                            Matthias Troyer <troyer@comp-phys.org>            (ALPS)\n"
      << "                  edited by Tama Ma         <tamama@hkusua.hku.hk>                  \n"
      << "************************************************************************************\n\n";
}


lowa::lowa(const alps::ProcessList& where,const alps::Parameters& p,int node)
  : alps::scheduler::MCRun(where,p,node),

    //Setup of general parameters
    label(static_cast<uint32_t>(p["LABEL"])),
    read_configuration(static_cast<int>(p["READ_CONFIGURATION"])),
    read_statistics(static_cast<int>(p["READ_STATISTICS"])),
    a_s(static_cast<double>(p["SCATTERING_LENGTH"])),
    mass(static_cast<double>(p["MASS"])),
    mu(static_cast<double>(p["CP"])),
    beta(static_cast<double>(p["BETA"])),
    time_of_flight(static_cast<double>(p["TIME_OF_FLIGHT"])),
    nmax(static_cast<SiteType>(p["MAXIMUM_NUMBER_PER_SITE"])),
    E_off(static_cast<double>(p["E_OFF"])),
    Ntest(static_cast<int32_t>(p["N_TEST"])),
    Ncan(static_cast<int32_t>(p["CANONICAL_VALUE"])),
    Nmeasure(static_cast<int32_t>(p["N_MEASURE"])),
    Nmeasure_green(static_cast<int32_t>(p["N_MEASURE_GREEN"])),

    //Configuration setup
    sweeps(0), thermalization_sweeps(static_cast<uint64_t>(p["THERMALIZATION"])),
    total_sweeps(static_cast<uint64_t>(p["SWEEPS"]))

{
  // Further initialization
  if (dim == 1) {
    vc     = new double;
    v0     = new double;
    waist  = new double;
    lambda = new double;
    phase  = new double;
    mid    = new double;
    Ls     = new SiteType;
  }
  else {
    vc = new double [dim];
    v0 = new double [dim];
    waist  = new double [dim];
    lambda = new double [dim];
    phase  = new double [dim];
    mid    = new double [dim];
    Ls = new SiteType [dim];
  }

// Please add accordingly to v0, waist, lambda here for dim other than 3 ...

  if (dim == 3) {
     vc[0] = static_cast<double>(p["VC_x"]);
     vc[1] = static_cast<double>(p["VC_y"]);
     vc[2] = static_cast<double>(p["VC_z"]);

     v0[0] = static_cast<double>(p["V0_x"]);
     v0[1] = static_cast<double>(p["V0_y"]);
     v0[2] = static_cast<double>(p["V0_z"]);

     waist[0] = static_cast<double>(p["WAIST_x"]);
     waist[1] = static_cast<double>(p["WAIST_y"]);
     waist[2] = static_cast<double>(p["WAIST_z"]);

     lambda[0] = static_cast<double>(p["WAVELENGTH_x"]);
     lambda[1] = static_cast<double>(p["WAVELENGTH_y"]);
     lambda[2] = static_cast<double>(p["WAVELENGTH_z"]);

     phase[0] = 87 * amu * lambda[0] * lambda[0] / (8*time_of_flight*hbar) * 1e-8;
     phase[1] = 87 * amu * lambda[1] * lambda[1] / (8*time_of_flight*hbar) * 1e-8;
     phase[2] = 87 * amu * lambda[2] * lambda[2] / (8*time_of_flight*hbar) * 1e-8;
  }

  for (int k = 0; k < dim; k++)    Ls[k] = static_cast<SiteType>(p["LENGTH"]);

  Nsites = 1;
  for (int k =0; k < dim; k++)  Nsites *= Ls[k];

  projected_Nsites = 1;
  for (int k=0; k < (dim-1); k++)  projected_Nsites *= Ls[k];

  mid[0] = ((Ls[0] % 2 == 0) ? Ls[0]/2 - 0.5 : (Ls[0]-1)/2. );
  mid[1] = ((Ls[1] % 2 == 0) ? Ls[1]/2 - 0.5 : (Ls[1]-1)/2. );
  mid[2] = ((Ls[2] % 2 == 0) ? Ls[2]/2 - 0.5 : (Ls[2]-1)/2. );


  // BH parm initialization
  t_x_plus  = new double [Nsites];
  t_x_minus = new double [Nsites];
  t_y_plus  = new double [Nsites];
  t_y_minus = new double [Nsites];
  t_z_plus  = new double [Nsites];
  t_z_minus = new double [Nsites];
  U         = new double [Nsites];
  epsilon   = new double [Nsites];

  cout << "\n Initializing Bose-Hubbard parameters...\t";
  BH_PARM BH_parm(label,Ls[0],v0[0],v0[1],v0[2],a_s,mass,lambda[0],lambda[1],lambda[2],waist[0],waist[1],waist[2],vc[0],vc[1],vc[2]);
  BH_parm.t_x_plus(t_x_plus);
  BH_parm.t_x_minus(t_x_minus);
  BH_parm.t_y_plus(t_y_plus);
  BH_parm.t_y_minus(t_y_minus);
  BH_parm.t_z_plus(t_z_plus);
  BH_parm.t_z_minus(t_z_minus);
  BH_parm.U(U);
  BH_parm.epsilon(epsilon);

  //t_hop = BH_parm.t_0;
  //U_on  = BH_parm.U_0;
  centre_site = BH_parm.centre_site;
  cout << "\t ...done.";

  // definition of bin parameters
  bin_size         = static_cast<TimeType>(p["BIN_SIZE"]); 
  total_no_of_bins = std::ceil(sqrt((Ls[0]-1-mid[0])*(Ls[0]-1-mid[0]) + (Ls[1]-1-mid[1])*(Ls[1]-1-mid[1]))/bin_size);
 
  projected_bin_no = new SiteType [projected_Nsites];
  bin_freq.resize(total_no_of_bins);

  for (SiteType j=0; j < Ls[1]; ++j) {
    SiteType p = j*Ls[0];
    for (SiteType i=0; i < Ls[0]; ++i) {
      SiteType cur_bin_no = std::floor(sqrt((i-mid[0])*(i-mid[0]) + (j-mid[1])*(j-mid[1]))/bin_size);
      projected_bin_no[p+i] = cur_bin_no;
      bin_freq[cur_bin_no] += 1;
    }
  }

  projected_binned_state.resize(total_no_of_bins);  
  projected_binned_state_times_N.resize(total_no_of_bins);

  // definitions of filenames
  filename0 = obtain_filename("qmc_profile_chem");
  filename1 = obtain_filename("qmc_saveconf");
  filename2 = obtain_filename("qmc_n_xyz_profile");
  filename3 = obtain_filename("qmc_nint_xy_profile");

  filename1_trial = obtain_filename("tqmc_saveconf");
  filename2_trial = obtain_filename("tqmc_profile");
  filename3_trial = obtain_filename("tqmc_nint_xy_profile");

  filenameN  = obtain_filename("N");

  //initialise parameters 
  time (&times1);
  timesm1 = times1;
  times_tot = 0;
  times_overall = 0;
  times_run = 0;
  init_internal();
  print_params();
  MCold = MCstep+MCstep_run;

  // alps measurement declaration
/*
  measurements << alps::RealObservable("Kinetic Energy");
  measurements << alps::RealObservable("Potential Energy");
  measurements << alps::RealObservable("Energy");
  measurements << alps::RealObservable("Particle Number");
  measurements << alps::RealObservable("Condensate Fraction");
  measurements << alps::RealObservable("Density at Centre");
*/
  measurements << alps::RealObservable("N");          // particle density
  measurements << alps::RealVectorObservable("n");    // density in cylindrical radial bins
  measurements << alps::RealVectorObservable("n^2");  // squared density "
  measurements << alps::RealVectorObservable("nN");   // density times N "


#ifndef TRAPPEDSYSTEM
  measurements << alps::RealObservable("Superfluid Density");
  measurements << alps::RealVectorObservable("Winding^2");

  if (dim == 1) {
    this_winding = new double;
    this_winding_squared = new double;
  }
  else {
    this_winding = new double[dim];
    this_winding_squared = new double[dim];
  }
#endif
}


void lowa::print_params() const
{
  std::cout << "\n Parameters : ";
  std::cout << "\n -------------";
  std::cout << "\n label           :  " << label;
  std::cout << "\n dim             :  " << dim;
  std::cout << "\n v0       [Er]   :  "; for(int k =0; k < dim; k++) cout << v0[k] << "\t";
  std::cout << "\n mu       [nK]   :  " << mu;
  std::cout << "\n vc       [nK]   :  "; for(int k =0; k < dim; k++) cout << vc[k] << "\t";
  std::cout << "\n waist    [um]   :  "; for(int k =0; k < dim; k++) cout << waist[k] << "\t";
  std::cout << "\n lambda   [nm]   :  "; for(int k =0; k < dim; k++) cout << lambda[k] << "\t";
  std::cout << "\n beta     [nK-1] :  " << beta;
  std::cout << "\n Ls              :  "; for(int k =0; k < dim; k++) cout << Ls[k] << "\t";
  std::cout << "\n Nsites          :  " << Nsites;
  std::cout << "\n nmax            :  " << nmax;
  std::cout << "\n E_off           :  " << E_off;
  //std::cout << "\n Ntherm          :  " << Ntherm;
  //std::cout << "\n Nloop           :  " << Nloop;
  std::cout << "\n Ntest           :  " << Ntest;
  if (Ncan >= 0) {
    std::cout << "\n Ncan            :  " << Ncan;
  }
  else {
    std::cout << "\n Grand can. ";
  }
  std::cout << "\n Nmeasure : " << Nmeasure;
  std::cout << "\n Nmeas_gr : " << Nmeasure_green;
  std::cout << "\n Nsave    : " << Nsave;
  std::cout << std::endl << std::endl << std::endl;
}


void lowa::init() { }


void lowa::init_internal()
{
  counter.resize(TEST+1);
  for (int i = 0; i < TEST+1; i++) counter[i] = 0;
  MCstep = 0;
  MCstep_run = 0.;
  MCmeas = 0;
  MCstep_total = 0.;
  acc_insert = 0.;
  tot_insert = 0.;
  worm_diag = true;
  new_measurement = true;
  if (dim == 1) {
    mWinding = new TimeType;
  }
  else {
    mWinding = new TimeType[dim];
  }
  for (int k = 0; k < dim; k++) mWinding[k] = 0.;

  trans.resize(zcoord + 1);
  nrvertex = 0;

  bool* border_site;
  border_site = new bool [Nsites];
  for (SiteType i = 0; i < Nsites; i++) border_site[i] = 0;


  //operator_string.resize(Nsites);
  operator_string = new list<element_type> [Nsites];
  //site_it.resize(Nsites);
  //dummy_it.resize(Nsites);
  site_it = new list<element_type>::iterator [Nsites];
  dummy_it = new list<element_type>::iterator [Nsites];
  for (SiteType i = 0 ; i < Nsites; i++) {
    site_it[i] = operator_string[i].begin();
  }
  mu_eff    = new double [Nsites];


  //state.resize(Nsites);
  //dns.resize(Nsites);
  //av_dns_inf.resize(Nsites);
  //av_state.resize(Nsites);
  //av_state_sq.resize(Nsites);
  state       = new SiteType [Nsites];
  //dns         = new TimeType [Nsites];

  av_dns        = new TimeType [Nsites];
  av_dns_inf    = new TimeType [Nsites];
/*
  av_state         = new TimeType [Nsites];
  av_state_sq      = new TimeType [Nsites];

  av_state_projected    = new TimeType [projected_Nsites];
  av_state_sq_projected = new TimeType [projected_Nsites];  
*/
  cdist.init(Nsites, dim);

  reset_av();

  nb.init(Nsites, zcoord);

  winding_element.init(zcoord, dim);

  //bond.init(Nsites, Nsites);
  //for (SiteType i = 0; i < Nsites; i++) {
  //  for (SiteType j = 0; j < Nsites; j++) {
  //    bond(i,j) = -1;
  //  }
  //}

  // specify lattice
  if (dim==1)
  {
    // version is not updated for dim =1 ***

    nb(0,0) = 1;
    nb(0,1) = Ls[0]-1;
    for(SiteType i = 1; i < Ls[0]-1; i++)
    {
      nb(i,0) = i+1;
      nb(i,1) = i-1;
    }
    nb(Ls[0]-1,0) = 0;
    nb(Ls[0]-1,1) = Ls[0] - 2;
    // read chemical potential from file
    double mx = Ls[0]/2.0 - 0.5;
    for (SiteType i =0; i < Ls[0]; i++) {
      double d = vc[0]*(i-mx)*(i-mx);
      d -= mu ;
      if ( (abs(vc[0]) > tol) && (i==0)) {
        d = 1000000. * abs(d);
      }
      else if ( (abs(vc[0]) > tol) && (i==Ls[0]-1)) {
        d = 1000000. * abs(d);
      }
      mu_eff[i] = d;
      //cout << "\n chempot : " << i << "\t" << d;
    }
    //for (SiteType i = 0; i < Ls[0]; i++) {
    // for (SiteType j = 0; j < Ls[0]; j++) {
    // cout << "\n dist : " << i << "\t" << j << "\t" <<  dist(i,j);
    //}
    //}
    winding_element(0,0) = 1.;
    winding_element(1,0) = -1.;
    //for (SiteType i = 0; i < Nsites; i++) {
    //  for (SiteType j = 0; j < zcoord; j++) {
        // bond(i, nb(i, j)) = j;
      // }
    // }
    border_site[0] = 1;
    border_site[Ls[0]-1] = 1;
  }
  else if (dim==2)
  {
    // version is not updated or dim =2 ***

    for (SiteType j =0; j < Ls[1]; j++) {
      SiteType y = j*Ls[0];
      for(SiteType i = 0; i < Ls[0]; i++) {
        if (i == 0) border_site[i+y] = 1;
        if (i == Ls[0]-1) border_site[i+y] = 1;
        if (j == 0) border_site[i+y] = 1;
        if (j == Ls[1]-1) border_site[i+y] = 1;
        i == Ls[0]-1 ? nb(i+y, 0) = y : nb(i+y, 0) = i+y+1;
        i == 0 ? nb(i+y, 2) = y + Ls[0]-1 : nb(i+y, 2) = i+y-1;
        j == Ls[1] - 1 ? nb(i+y, 1) = i : nb(i+y,1) = (j+1)*Ls[0] + i;
        j == 0 ? nb(i+y, 3) = (Ls[1]-1)*Ls[0] + i : nb(i+y,3) = (j-1)*Ls[0] + i;
      }
    }
    // read chemical potential from file
    double my = Ls[1]/2.0 - 0.5;
    double mx = Ls[0]/2.0 - 0.5 ;
    for (SiteType j =0; j < Ls[1]; j++) {
      SiteType y = j*Ls[0];
      for (SiteType i =0; i < Ls[0]; i++) {
        SiteType p = i + y;
        double d = vc[1]*(j-my)*(j-my) + vc[0]*(i-mx)*(i-mx);
        if ( (abs(vc[0]) > tol)) {
          if ((j==0) || (j==Ls[1]-1) || (i==0) || (i==Ls[0]-1) ) d = 1000000.*abs(d);
        }
        d -= mu;
        mu_eff[p] = d;
        //cout << "\n chempot : " << i << "\t" << j <<  "\t" << d;
      }
    }
    winding_element(0,0) = 1.;
    winding_element(0,1) = 0.;
    winding_element(1,0) = 0.;
    winding_element(1,1) = 1.;
    winding_element(2,0) = -1.;
    winding_element(2,1) = 0.;
    winding_element(3,0) = 0.;
    winding_element(3,1) = -1.;
    // for (SiteType i = 0; i < Nsites; i++) {
      // for (SiteType j = 0; j < zcoord; j++) {
        // bond(i, nb(i, j)) = j;
      // }
    // }
  }
  else // dim==3
  {
    for (SiteType k = 0; k < Ls[2]; k++) {
      SiteType z = k*Ls[0]*Ls[1];
      for (SiteType j =0; j < Ls[1]; j++) {
        SiteType y = j*Ls[0];
        for (SiteType i =0; i < Ls[0]; i++) {
          SiteType p = i + y + z;
          if (i == 0) border_site[p] = 1;
          if (i == Ls[0]-1) border_site[p] = 1;
          if (j == 0) border_site[p] = 1;
          if (j == Ls[1]-1) border_site[p] = 1;
          if (k == 0) border_site[p] = 1;
          if (k == Ls[2]-1) border_site[p] = 1;
          i == Ls[0]-1 ? nb(p,0) = y+z : nb(p,0) = p+1;
          i == 0 ? nb(p,3) = y + z + Ls[0]-1 : nb(p,3) = p- 1;
          j == Ls[1]-1 ? nb(p,1) = z+ i : nb(p,1) = z + i + (j+1)*Ls[0] ;
          j == 0 ? nb(p,4) = z + i + (Ls[1]-1)*Ls[0] :  nb(p,4) = z + i + (j-1)*Ls[0];
          k == Ls[2]-1 ? nb(p,2) = y+i : nb(p,2) = y+i+ (k+1)*Ls[0]*Ls[1];
          k == 0 ? nb(p,5) = y + i + (Ls[2]-1)*Ls[0]*Ls[1] : nb(p,5) = y + i + (k-1)*Ls[0]*Ls[1];
        }
      }
    }
    // read chemical potential from file
/*
    double mz = Ls[2]/2.0 - 0.5;
    double my = Ls[1]/2.0 - 0.5;
    double mx = Ls[0]/2.0 - 0.5;
    for (SiteType k = 0; k < Ls[2]; k++) {
      SiteType z = k*Ls[0]*Ls[1];
      for (SiteType j =0; j < Ls[1]; j++) {
        SiteType y = j*Ls[0];
        for (SiteType i =0; i < Ls[0]; i++) {
          SiteType p = i + y + z;
          double d = vc[2]*(k-mz)*(k-mz) + vc[1]*(j-my)*(j-my) + vc[0]*(i-mx)*(i-mx);
          if ( (abs(vc[0]) > tol)) {
            if ((j==0) || (j==Ls[1]-1) || (i==0) || (i==Ls[0]-1) || (k==0) || (k==Ls[2]-1) ) d = 1000000.;
          }
          d -= mu ;
          mu_eff[p] = d;
          //cout << "\n chempot : " << i << "\t" << j << "\t" << k << "\t" <<d;
        }
      }
    }
*/
    for (SiteType k = 0; k < Ls[2]; k++) {
       SiteType z = k*Ls[0]*Ls[1];
       for (SiteType j = 0; j < Ls[1]; j++) {
          SiteType y = j*Ls[0];
          for (SiteType i = 0; i < Ls[0]; i++) {
             SiteType p = i + y + z;
              
             cdist(p,0) = (i-mid[0]) * (i-mid[0]);
             cdist(p,1) = (j-mid[1]) * (j-mid[1]);
             cdist(p,2) = (k-mid[2]) * (k-mid[2]);

          }
       }
    }


    for (SiteType p = 0; p < Nsites; p++) {
        mu_eff[p] = (epsilon[p] - mu);
    }

    winding_element(0,0) = 1.;
    winding_element(0,1) = 0.;
    winding_element(0,2) = 0.;
    winding_element(1,0) = 0.;
    winding_element(1,1) = 1.;
    winding_element(1,2) = 0.;
    winding_element(2,0) = 0.;
    winding_element(2,1) = 0.;
    winding_element(2,2) = 1.;
    winding_element(3,0) = -1.;
    winding_element(3,1) = 0.;
    winding_element(3,2) = 0.;
    winding_element(4,0) = 0.;
    winding_element(4,1) = -1.;
    winding_element(4,2) = 0.;
    winding_element(5,0) = 0.;
    winding_element(5,1) = 0.;
    winding_element(5,2) = -1.;
    // for (SiteType i = 0; i < Nsites; i++) {
      // for (SiteType j = 0; j < zcoord; j++) {
        // bond(i, nb(i, j)) = j;
      // }
    // }
  }


  itherm = 0;
  iloop = 0;

  std::ifstream inFile(filename1.c_str(), ios::binary);

  if (inFile.good()) {
    std::cout << "\n# Retrieving old configuration....";
    load_internal(inFile, counter);
    if (!read_statistics) reset_av();
    inFile.close();
    cout << "\n# Time elapsed in previous jobs : " << times_run << "\t" << times_overall;
    cout << "\n# Itherm : " << itherm << "\t" << " iloop : " << iloop;
    std::cout << "\n# Retrieving old configuration done. Testing...";
    if (!test_conf()) {char ch; cin >> ch;};
    cout << "\n# Testing OK...\n";
    update_en();
    cout << "# Number of particles : " << number_of_particles << endl;
    cout << "# Potential Energy    : " << Epot << endl;
    cout << "# Kinetic Energy      : " << Ekin << endl << endl << endl;
    cout << "# Worm head           : " << worm_head << endl;
    cout << "# Worm tail           : " << worm_tail << endl;
    cout << "# wormdiag - wormdir  : " << worm_diag << "\t" << worm_dir << endl;
    cout << "# wormrising          : " << worm_rising << endl;
  }

  else {
    SiteType i = 0;
#ifdef HARDCORE
    for (SiteType s = 0; s < Nsites; s++) {
      state[s] = ( random_real() < 0.5 ? 0 : 1);
    }
#else
    //for (vector<SiteType>::iterator it = state.begin(); it !=state.end(); ++it, ++i)   *it = random_int();
    //for (vector<SiteType>::iterator it = state.begin(); it !=state.end(); ++it, ++i)   *it = 0;

    for (SiteType s = 0; s < Nsites; s++) {
       state[s] = 0;
    }
#endif


#ifdef TRAPPEDSYSTEM
    for (SiteType i = 0; i < Nsites; i++) {
      if (border_site[i] == 1) {
        state[i] = 0;
        mu_eff[i] = abs(mu_eff[i])*10000;
      }
    }
#endif
    // insert dummy elements on all sites, keeping the code ler and making a line of where to measure diag properties
    for (SiteType i = 0; i < Nsites; i++) {
      Element new_elem(state[i],state[i], i, i, beta);
      site_it[i] = operator_string[i].insert(site_it[i], new_elem);
      dummy_it[i] = site_it[i];
#ifdef DEBUGMODE
      cout << "\nInserting dummy : " << i << "\t " << site_it[i]->time() << "\t" << site_it[i]->before();
#endif
    }
    for (SiteType i = 0; i < Nsites; i++) {
      for (SiteType j = 0; j < zcoord; j++) {
        site_it[i]->set_assoc(j, site_it[nb(i,j)]);
      }
    }
#ifdef DEBUGMODE
    //print_conf(cout);
#endif
  }

  std::ifstream inFile2(filename2.c_str(), ios::in);
  std::ifstream inFile3(filename3.c_str(), ios::in);
  if (!inFile2.good() || !inFile3.good()) {
     reset_av();
  }
  

#ifdef TRAPPEDSYSTEM
  delete [] border_site;
#endif
  cout << "\n# Finished intializing.\n";
}


int lowa::dostep_internal()
{
  if (worm_diag) {
  char ch;
  int start_occ;

  //for (int32_t i = 0; i < Nsites; i++) dns[i] = 0.;
  mZ_dns += 1.;
  random_real() < 0.5 ? worm_dir = -1 : worm_dir = +1;
  SiteType isite=int32_t(random_real()*Nsites);
  TimeType start_time=random_real()*beta;

  // find suitable place to make a worm pair
  if (operator_string[isite].size() != 1) // dummy element always there
  {
    list<element_type>::iterator prevworm_it;
    prevworm_it = ( site_it[isite] == operator_string[isite].begin() ? operator_string[isite].end() : site_it[isite] );
    --prevworm_it;

    // you cannot create a worm at the same place where the previous one was removed...
    // you really need to go around in space-time. maybe for insulating phases
    // the next couple of lines are time consuming
    while (!t_between(start_time, prevworm_it->time(), site_it[isite]->time()))
    {
      prevworm_it = site_it[isite];
      site_right(isite);
    }
  }
  start_occ = site_it[isite]->before();
  worm_it = site_it[isite];

  bool l_t, mkworm = 1;
  tot_insert += 1.;
  // detailed balance slightly differs for hard-core bosons
  // see also density matrix for that
#ifdef HARDCORE
  l_t = 1;
  worm_rising = (start_occ == 0 ? 0 : 1);
#else

  if ((start_occ > 0) && (start_occ < nmax)) {
    random_real() < 0.5 ? worm_rising = 1 : worm_rising = 0;
  }
  else if (start_occ == 0) {
    if (random_real() < 0.5) {
      worm_rising = 0;
    }
    else {
      mkworm = 0;
      l_t = 1;
    }
  }
  else { // startocc equals nmax, worm should only be created with prob 1/2
    if (random_real() < 0.5) {
      worm_rising = 1;
    }
    else {
      mkworm = 0;
      l_t = 0;
    }
  }

  MCstep++;

  if (!mkworm) // density matrix is a bit tricky
  {
//  dns[0] = (l_t ? 0. : nmax);
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
    av_dns[0] += (l_t ? 0. : nmax);
    av_dns_inf[0] += (l_t ? 0. : nmax);
#else
    av_dns[0] += (l_t ? 0. : nmax);
    av_dns_inf[0] += (l_t ? 0. : nmax);
#endif
#else
#ifdef TRAPPEDSYSTEM
    av_dns[0] += (l_t ? 0. : nmax);
    av_dns_inf[0] += (l_t ? 0. : nmax); //dns[0];
#else
    av_dns[0] += (l_t ? 0. : nmax);
    av_dns_inf[0] += (l_t ? 0. : nmax); //dns[0];
#endif
#endif
    //cout << "\n No worm inserted, mkworm = 0";

    return(0);
  }

#endif

  acc_insert += 1.;
  new_measurement = true;

  // we can now really create the worms
  worm_tail.from(isite);
  worm_tail.to(isite);
  worm_tail.time(start_time);
  worm_head = worm_tail;

  worm_diag = 0;

  // the worm tail is added to the interaction list.
  // the worm head is not, but we always need to know
  // where it is and what the density just before/after
  // the worm head is
  if (worm_dir == +1) {
    if ( worm_rising ) {
      start_bos = start_occ;
      worm_tail.before(start_occ);
      worm_tail.after(--start_occ);
      worm_head.before(worm_tail.after());
      worm_head.after(worm_tail.before());
    }
    else {
      worm_tail.before(start_occ);
      worm_tail.after(++start_occ);
      worm_head.before(worm_tail.after());
      worm_head.after(worm_tail.before());
      start_bos = start_occ;
    }
    site_it[isite] = operator_string[isite].insert(site_it[isite], worm_tail);
    find_assoc_insert(isite, site_it[isite]);
    site_right(isite);
    worm_it = site_it[isite];
  }
  else {
    if ( worm_rising ) {
      start_bos = start_occ;
      worm_tail.after(start_occ);
      worm_tail.before(--start_occ);
      worm_head.before(worm_tail.after());
      worm_head.after(worm_tail.before());
    }
    else {
      worm_tail.after(start_occ);
      worm_tail.before(++start_occ);
      worm_head.before(worm_tail.after());
      worm_head.after(worm_tail.before());
      start_bos = start_occ;
    }
    site_it[isite] = operator_string[isite].insert(site_it[isite], worm_tail);
    find_assoc_insert(isite, site_it[isite]);
    site_left(isite);
    worm_it = site_it[isite];
  }
#ifdef DEBUGMODE
  cout << "\nInitial head : " << worm_head;
  cout << "\nInitial tail : " << worm_tail;
  cout << "\nInitialized worm..." << worm_tail << endl;
  cout << "\nworm_rising " << worm_rising << "\tdir " << worm_dir;
  //cin>>ch;
#endif


#ifdef DEBUGMODE
  //print_conf(cout);
#endif
  }


  do {
    TimeType pexp = -log(random_real());
#ifdef DEBUGMODE
    //cout << "\nExponential jump time : " << pexp  << " dir " << worm_dir<< "\t" << worm_rising << endl;
#endif
    while (pexp!=0.) {  //exponential time jump

      MCstep ++;
      worm_cycle(pexp);

    }
    if (MCstep >= Nsave) {
      return (1);
    }

  } while (!worm_diag);


#ifdef DEBUGMODE
  cout << "\n*** DIAGONAL *** ";
  print_conf(cout);
  cout << "\n*****************";
  if (!test_conf()) {char ch; cin >> ch;};
#endif
  return(0);
}



void lowa::worm_cycle(TimeType& p)
{

  TimeType El, Er, dE;
  TimeType segm_tau;

  // we distinguish between moving to the right and moving to the left
  if (worm_dir > 0) {

    // find out the energies before and after the worm
    if (worm_rising) {
      dE = U[worm_head.to()] * (site_it[worm_head.to()]->before() - 1) + mu_eff[worm_head.to()];
      dE = ( dE > 0 ? E_off : -dE + E_off);
    }
    else {
      dE =  -U[worm_head.to()] * (site_it[worm_head.to()]->before())- mu_eff[worm_head.to()];
      dE = (dE > 0 ? E_off : -dE + E_off);
    }
#ifdef DEBUGMODE
    cout << "\n moving right p = " << p << "\tdE " << dE << "\tdens " << site_it[worm_head.to()]->before() << "\tmu " << mu_eff[worm_head.to()] << "\tworm_it " << *worm_it << "\tsite_it " << *site_it[worm_head.to()];
    cout << "\n Worm : " << worm_head.to() << "\t" << worm_head.time();
    if (worm_it != site_it[worm_head.to()]) cout << "\nWRONG assignment of worm_it";
#endif

    if (shift_right_kink(p, dE, segm_tau) ) {
#ifdef DEBUGMODE
      cout << "\n moving right p = " << p << "\tdE " << dE << "\tsegm_tau " << segm_tau;
#endif
      if ((t_between(worm_tail.time(), worm_head.time(), worm_it->time())) && (worm_tail.time()!=(worm_it->time() ))) {
#ifdef DEBUGMODE
	cout << "\nMoving upwards, measuring density matrix " << worm_tail.time() << "\t" << worm_head.time() << "\t" << worm_it->time();
#endif

	// measure density matrix
	worm_head.time(worm_tail.time());
	//dns[dist(worm_head.to(), worm_tail.to())] += start_bos;
	if (worm_head.to() != worm_tail.to()) {
          double cosphi=1.; for (int k=0; k < dim; k++) cosphi *= cos(phase[k]* (cdist(worm_head.to(),k) - cdist(worm_tail.to(),k)));
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dns[dist(worm_head.to(), worm_tail.to())] += 0.25*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.25*start_bos;
#else
          av_dns[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos;
#endif
#else
#ifdef TRAPPEDSYSTEM       
          av_dns[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos * cosphi;
          av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos;
#else
          av_dns[dist(worm_head.to(), worm_tail.to())] += 1.*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 1.*start_bos;
#endif 
#endif
        }
	p = 0.;
	return;
      }

      worm_head.time(worm_it->time());  //worm jumps to next element

      if (worm_it->from() == worm_it->to()) // the time of the worm head is the same as the time of the worm tail
      {
#ifdef DEBUGMODE
        cout << "\nMoving upwards, time head = time tail" << "\t" << worm_head << "\tworm_it " << *worm_it;
#endif
        if (worm_head.time() == worm_tail.time()) { // worm bites in tail (same sites)
#ifdef DEBUGMODE
	  cout << "\n ...removing worm";
#endif
	  //dns[0] = worm_it->after();  
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dns[0] += worm_it->after();
	  av_dns_inf[0] += worm_it->after();
#else
          av_dns[0] += worm_it->after();
	  av_dns_inf[0] += worm_it->after();
#endif
#else
#ifdef TRAPPEDSYSTEM           
          av_dns[0] += worm_it->after();
	  av_dns_inf[0] += worm_it->after();
#else
          av_dns[0] += worm_it->after();
	  av_dns_inf[0] += worm_it->after();
#endif 
#endif


	  find_assoc_delete(worm_head.to(), site_it[worm_head.to()]);
          site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
          if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();
	  worm_it = site_it[worm_head.to()];
          worm_diag = 1;

          p = 0.;
          return;
        }
        else {  // worm head passes dummy
#ifdef DEBUGMODE
	  cout << "\n ...passing dummy" << worm_head.to() << "\t" << worm_tail.to() << "\trising : "<< worm_rising;
#endif

	  if (worm_rising) {
	    site_it[worm_head.to()]->before(site_it[worm_head.to()]->before()-1);
	    site_it[worm_head.to()]->after( site_it[worm_head.to()]->after() -1);
	  }
	  else {
	    site_it[worm_head.to()]->before(site_it[worm_head.to()]->before()+1);
	    site_it[worm_head.to()]->after( site_it[worm_head.to()]->after() +1);
	  }
#ifdef DEBUGMODE
	  cout << "\n ...passed dummy" << *worm_it << "\t" << *site_it[worm_head.to()];
#endif
	  site_right(worm_head.to());
	  worm_it = site_it[worm_head.to()];
#ifdef DEBUGMODE
	  cout << "\n ...passed dummy" << *worm_it;
#endif
          p -= segm_tau * dE;
          return;
        } // else ... worm passes tail
      } // time of worm head equals time of worm tail  
      else if (worm_it->from() == worm_head.to())
      {
#ifdef DEBUGMODE
        cout << "\nMoving upwards, next.from = head.to";
#endif
        if (!worm_rising) {  //pass the interaction with probability 1
	  // no iterators change, worm remains on same site
	  site_it[worm_head.to()]->before(site_it[worm_head.to()]->before()+1);
	  site_it[worm_head.to()]->after( site_it[worm_head.to()]->after() +1);
	  site_right(worm_head.to());
	  worm_it = site_it[worm_head.to()];
          p -= segm_tau*dE;
          return;
        }
        else // try to delete/relink interaction
        {
	  for (SiteType j = 0; j < zcoord; j++) {
	    list<element_type>::iterator it = site_it[worm_it->from()]->get_assoc(j);
	    if (worm_it->time() == it->time()) {
	      site_it[worm_it->to()] = it;
	      break;
	    }
	  }
	  //cout << "\n calc_del_weight iterator for to" << site_it[worm_it->to()]->time() << "\t" << &*site_it[worm_it->to()];
          calc_del_weight(worm_it->to()); // this should also set the iterators right
          SiteType i3 = heatbath_search(trans);
          if (i3 == zcoord) { //annihilate interaction
#ifdef DEBUGMODE
            cout << "\n  annihilate " << i3;
#endif
	    SiteType ib = worm_it->to();

	    // delete interaction on the from site
	    find_assoc_delete(worm_it->from(), site_it[worm_it->from()]);
	    find_assoc_delete(worm_it->to(), site_it[worm_it->to()]);
	    site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();

	    // delete interaction on the to site
            worm_head.to(ib);
            site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();
	    worm_it = site_it[worm_head.to()]; // worm jumps
            p -= segm_tau * dE;
	    nrvertex--;
	    //print_conf(cout);
            return;
          }
          else //relink
          {
#ifdef DEBUGMODE
            cout << "\n relink or bounce " << i3 << "\t" << nb(worm_it->to(),i3);
#endif
	    // define the new element
	    SiteType ito = worm_it->to();
            SiteType ib = nb(ito,i3);


	    site_it[ito]->from(ib);

	    // delete the old element on the from site
	    find_assoc_delete(worm_it->from(), site_it[worm_it->from()]);
	    site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();

	    // insert new element
	    int n0 = site_it[ib]->before();
	    element_type new_elem(n0+1,n0, ib, ito, worm_head.time() );
	    site_it[ib] = operator_string[ib].insert(site_it[ib], new_elem);
	    //cout << "\nafter insert, before find_assoc!!! " << site_it[ib]->time() << "\t" << site_it[ito]->time();
	    find_assoc_insert(ib, site_it[ib]);

	    worm_head.to(ib); // worm jumps
	    site_left(ib);
	    worm_it = site_it[ib];
	    worm_dir = -worm_dir;
	    worm_rising = !worm_rising;
            p -= segm_tau * dE;
            return;
          }
        } // ... else delete/relink interaction
      }
      else if (worm_it->to() == worm_head.to()) {
#ifdef DEBUGMODE
	cout << "\nMoving upwards next.to = head.to";
#endif
        if (worm_rising) { // pass interaction with probability 1
	  // no iterators change, worm remains on same site
	  site_it[worm_head.to()]->before(site_it[worm_head.to()]->before()-1);
	  site_it[worm_head.to()]->after( site_it[worm_head.to()]->after() -1);
	  site_right(worm_head.to());
	  worm_it = site_it[worm_head.to()];
          p -= segm_tau * dE;
          return;
        }
        else { // try to delete/relink interaction
	  for (SiteType j = 0; j < zcoord; j++) {
	    list<element_type>::iterator it = site_it[worm_it->to()]->get_assoc(j);
	    if (worm_it->time() == it->time()) {
	      site_it[worm_it->from()] = it;
	      break;
	    }
	  }
	  calc_del_weight(worm_it->from());
          SiteType i3 = heatbath_search(trans);

          if (i3 == zcoord) { //annihilate interaction

            SiteType ib = worm_it->from();

#ifdef DEBUGMODE
	    cout << "\n Annihilate " << i3 << endl;
#endif

	    // delete element on the to site
	    find_assoc_delete(worm_it->to(), site_it[worm_it->to()]);
	    find_assoc_delete(worm_it->from(), site_it[worm_it->from()]);
	    site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();

	    // delete element on the from site
            worm_head.to(ib);
            site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();
	    worm_it = site_it[worm_head.to()]; // worm jumps

            p -= segm_tau * dE;
	    nrvertex--;
            return;
          }
          else {// relink the interaction
#ifdef DEBUGMODE
	    cout << "\n Relink : " << i3 << endl;
#endif
	    // define the new element
	    SiteType ifrom = worm_it->from();
            SiteType ib = nb(ifrom,i3);

	    site_it[ifrom]->to(ib);

	    // delete the old element on the to site
	    find_assoc_delete(worm_it->to(), site_it[worm_it->to()]);
	    site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();

	    // insert new element
	    int n0 = site_it[ib]->before();
	    element_type new_elem(n0-1,n0, ifrom, ib, worm_head.time() );
	    site_it[ib] = operator_string[ib].insert(site_it[ib], new_elem);
	    find_assoc_insert(ib, site_it[ib]);

	    worm_head.to(ib); // worm jumps
	    site_left(ib);
	    worm_it = site_it[ib];
	    worm_dir = -worm_dir;
	    worm_rising = !worm_rising;
            p -= segm_tau*dE;
	    return;
          }
        }  //....try to delete/relink interaction
      }
      else // worm passes the interaction (no sites in common), but diagonal energy might change
      {
#ifdef DEBUGMODE
        cout << "\nMoving upwards, worm passes remote interaction. This should not occur at all!!";
	//print_conf(cout);
	char ch; cin >> ch;
#endif
      }
    } // shift_right_kink

    else { // no interaction reached, try to insert new interaction

#ifdef DEBUGMODE
      cout << "\n moving right p = " << p << "\t dE " << dE << "\t segm_tau" << segm_tau;
      cout << "\nMoving to the right, inserting interaction " << worm_head.to() << "\t" << worm_rising ;
#endif
      TimeType newtime = worm_head.time() + p/dE;
      while (newtime > beta) {
         newtime -= beta;
      }

      if (  t_between(worm_tail.time(), worm_head.time(), newtime) ) { 
#ifdef DEBUGMODE
	cout << "\nMoving to the right, not inserting interaction, but measuring the density matrix ";
#endif

	// measure density matrix

	worm_head.time(worm_tail.time());
	//dns[dist(worm_head.to(), worm_tail.to())] += start_bos;
	if (worm_head.to() != worm_tail.to()) {
          double cosphi=1.; for (int k=0; k < dim; k++) cosphi *= cos(phase[k]* (cdist(worm_head.to(),k) - cdist(worm_tail.to(),k)));
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dns[dist(worm_head.to(), worm_tail.to())] += 0.25*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.25*start_bos;
#else
          av_dns[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos;
#endif
#else
#ifdef TRAPPEDSYSTEM           
          av_dns[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos;
#else
          av_dns[dist(worm_head.to(), worm_tail.to())] += 1.*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 1.*start_bos;
#endif 
#endif
	}
	p = 0.;
	//print_conf(cout);
	return;
      }

      calc_ins_weight(dE, worm_head.to(), newtime); // this should also set the iterators right...
      worm_head.time(newtime);

      SiteType i3 = heatbath_search(trans);
#ifdef DEBUGMODE
      cout << "\n ..." << i3 << "\t" << nb(worm_head.to(),i3);
#endif
      if ( i3 < zcoord)
      {
        SiteType ib = nb(worm_head.to(),i3);

        element_type new_elem;
        new_elem.time(worm_head.time());
        if (worm_rising) {
	  new_elem.to(worm_head.to());
	  new_elem.from(ib);
	  new_elem.after(site_it[worm_head.to()]->before());
	  new_elem.before(site_it[worm_head.to()]->before()-1);
	  site_it[worm_head.to()] = operator_string[worm_head.to()].insert(site_it[worm_head.to()], new_elem);
	  new_elem.after(site_it[ib]->before()-1);
	  new_elem.before(site_it[ib]->before());
	  site_it[ib] = operator_string[ib].insert(site_it[ib], new_elem);
        }
	else {
          new_elem.to(ib);
          new_elem.from(worm_head.to());
	  new_elem.after(site_it[worm_head.to()]->before());
	  new_elem.before(site_it[worm_head.to()]->before()+1);
	  site_it[worm_head.to()] = operator_string[worm_head.to()].insert(site_it[worm_head.to()], new_elem);
	  new_elem.after(site_it[ib]->before()+1);
	  new_elem.before(site_it[ib]->before());
	  site_it[ib] = operator_string[ib].insert(site_it[ib], new_elem);
	}
	find_assoc_insert(worm_head.to(), site_it[worm_head.to()]);
	find_assoc_insert(ib, site_it[ib]);
	site_right(worm_head.to());
	site_right(ib);

	// worm jumps
        worm_head.to(ib);
	worm_it = site_it[ib];
        p = 0.;
	nrvertex++;
        return;
      }
      else // bounce back
      {
#ifdef DEBUGMODE
        cout << "\nInserting interaction failed, bouncing back, moving upwards";
#endif

	site_left(worm_head.to());
	worm_it = site_it[worm_head.to()];
	worm_dir = -worm_dir;
	worm_rising = !worm_rising;
        p = 0.;
        return;
      }
    }
   } // worm_dir > 0
   
//*************************************************************************************
  else { // worm_dir < 0; this is mirror symmetry of the case worm_dir > 0. it's all about getting the +1 and -1 right
    if (worm_rising) {
      dE = U[worm_head.to()] * (site_it[worm_head.to()]->after() - 1) + mu_eff[worm_head.to()];
      dE = (dE > 0 ? E_off : -dE + E_off);
    }
    else {
      dE = -U[worm_head.to()] * (site_it[worm_head.to()]->after()) - mu_eff[worm_head.to()];
      dE = (dE > 0 ? E_off : -dE + E_off);
    }
#ifdef DEBUGMODE
    cout << "\n moving left p = " << p << "\tdE " << dE << "\tdens " << site_it[worm_head.to()]->after() << "\tmu " << mu_eff[worm_head.to()] << "\tworm_it " << *worm_it << "\tsite_it " << *site_it[worm_head.to()];
    cout << "\n Worm : " << worm_head.to() << "\t" << worm_head.time();
    //char chr; cin>> chr;
    if (worm_it != site_it[worm_head.to()]) cout << "\nWRONG assignment of worm_it";
#endif

    if (shift_left_kink(p, dE, segm_tau)) {

      if ( ( t_between(worm_tail.time(), worm_it->time(), worm_head.time() )) && ( worm_tail.time()  !=  worm_it->time() ) && (worm_tail.time() != worm_head.time()) ) { 
	// measure density matrix
#ifdef DEBUGMODE
	cout << "\n measuring density matrix";
#endif

	worm_head.time(worm_tail.time());
	//dns[dist(worm_head.to(), worm_tail.to())] += start_bos;
	if (worm_head.to() != worm_tail.to()) {
          double cosphi=1.; for (int k=0; k < dim; k++) cosphi *= cos(phase[k]* (cdist(worm_head.to(),k) - cdist(worm_tail.to(),k)));
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dns[dist(worm_head.to(), worm_tail.to())] += 0.25*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.25*start_bos;
#else
          av_dns[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos;
#endif
#else
#ifdef TRAPPEDSYSTEM           
	  av_dns[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos;
#else
          av_dns[dist(worm_head.to(), worm_tail.to())] += 1.*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 1.*start_bos;
#endif 
#endif
	}
	p = 0.;
	return;
      }

      worm_head.time(worm_it->time());  //worm jumps to next element
      if (worm_it->from() == worm_it->to()) { // the time of the worm head is the same as the time of the worm tail
#ifdef DEBUGMODE
        cout << "\nMoving downwards, time head equals time tail";
#endif
        if (worm_head.time() == worm_tail.time()) // worm bites in tail (same sites)
        {
#ifdef DEBUGMODE
          cout << "\n  remove tail in moving downwards";
#endif
          //dns[0] = worm_it->before(); 
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dns[0] += worm_it->before();
	  av_dns_inf[0] += worm_it->before();
#else
          av_dns[0] += worm_it->before();
	  av_dns_inf[0] += worm_it->before();
#endif
#else
#ifdef TRAPPEDSYSTEM           
          av_dns[0] += worm_it->before();
	  av_dns_inf[0] += worm_it->before();
#else
          av_dns[0] += worm_it->before();
	  av_dns_inf[0] += worm_it->before();
#endif 
#endif

	  find_assoc_delete(worm_head.to(), site_it[worm_head.to()]);

	  site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
	  if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();
	  worm_it = site_it[worm_head.to()];
          worm_dir = 1;
          worm_diag = 1;
          p = 0.;
          return;
        }
	else {  // worm head passes tail on other site, update green function (different sites)
#ifdef DEBUGMODE
	  cout << "\n  passing dummy in moving downwards";
#endif

	  if (worm_rising) {
	    site_it[worm_head.to()]->before(site_it[worm_head.to()]->before()-1);
	    site_it[worm_head.to()]->after(site_it[worm_head.to()]->after() -1);
	  }
	  else {
	    site_it[worm_head.to()]->before(site_it[worm_head.to()]->before()+1);
	    site_it[worm_head.to()]->after(site_it[worm_head.to()]->after() +1);
	  }
	  site_left(worm_head.to());
	  worm_it = site_it[worm_head.to()];
	  p -= segm_tau * dE;
          //print_conf(cout);
	  return;
	} // else ... worm passes tail
      } // time of worm head equals time of worm tail
      else if (worm_it->from() == worm_head.to()) {
#ifdef DEBUGMODE
	cout << "\nMoving downwards, prev.from = head.to " << worm_rising;
#endif
        if (worm_rising) { // pass interaction with probability 1
	  site_it[worm_head.to()]->after( site_it[worm_head.to()]->after()  - 1);
	  site_it[worm_head.to()]->before(site_it[worm_head.to()]->before() - 1);
	  site_left(worm_head.to());
	  worm_it = site_it[worm_head.to()];
          p -= segm_tau * dE;
          return;
        }
        else { // try to delete/relink interaction
	  for (SiteType j = 0; j < zcoord; j++) {  
	    list<element_type>::iterator it = site_it[worm_it->from()]->get_assoc(j);
	    if (worm_it->time() == it->time()) {
	      site_it[worm_it->to()] = it;
	      break;
	    }
	  }
	  calc_del_weight(worm_it->to());        // This should also put the iterators right
          SiteType i3 = heatbath_search(trans);

          if (i3 == zcoord) {  
	    //annihilate interaction
#ifdef DEBUGMODE
            cout << "\n  annihilate " << i3;
#endif

	    SiteType ib = worm_it->to();


	    // delete interaction on the from site
	    find_assoc_delete(worm_it->from(), site_it[worm_it->from()]);
	    find_assoc_delete(worm_it->to(), site_it[worm_it->to()]);
	    site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();
	    site_left(worm_head.to());

	    // delete interaction on the to site
            worm_head.to(ib);
            site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();
	    site_left(worm_head.to());

	    worm_it = site_it[worm_head.to()]; // worm jumps
            p -= segm_tau * dE;
	    nrvertex--;
	    //print_conf(cout);
            return;
          }
          else { //relink
#ifdef DEBUGMODE
	    cout << "\n relink or bounce " << i3 << "\t" << nb(worm_it->to(),i3);
#endif
	    // define the new element
	    SiteType ito = worm_it->to();
            SiteType ib  = nb(ito,i3);

	    site_it[ito]->from(ib);

	    // delete the old element on the from site
	    find_assoc_delete(worm_it->from(), site_it[worm_it->from()]);
	    site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();

	    // insert new element
	    int n0 = site_it[ib]->before(); 
	    element_type new_elem(n0,n0-1, ib, ito, worm_head.time() );
	    site_it[ib] = operator_string[ib].insert(site_it[ib], new_elem);
	    find_assoc_insert(ib, site_it[ib]);

	    worm_head.to(ib); // worm jumps
	    site_right(ib);
	    worm_it = site_it[ib];
	    worm_dir = -worm_dir;
	    worm_rising = !worm_rising;
            p -= segm_tau * dE;
            return;
          }
        }  //....try to delete/relink interaction
      }
      else if (worm_it->to() == worm_head.to()) {
#ifdef DEBUGMODE
        cout << "\nMoving downwards, prev.to = head.to";
#endif
        if (!worm_rising) {  //pass the interaction with probability 1
	  site_it[worm_head.to()]->after( site_it[worm_head.to()]->after()  + 1);
	  site_it[worm_head.to()]->before(site_it[worm_head.to()]->before() + 1);
	  site_left(worm_head.to());
	  worm_it = site_it[worm_head.to()];
          p -= segm_tau * dE;
          return;
        }
        else { // try to delete/relink interaction
	  for (SiteType j = 0; j < zcoord; j++) {  
	    list<element_type>::iterator it = site_it[worm_it->to()]->get_assoc(j);
	    if (worm_it->time() == it->time()) {
	      site_it[worm_it->from()] = it;
	      break;
	    }
	  }
	  //cout << "\n calc_del_weight iterator for from" << site_it[worm_it->from()]->time() << "\t" << &*site_it[worm_it->from()];
	  calc_del_weight(worm_it->from());
          SiteType i3 = heatbath_search(trans);
          if (i3 == zcoord) { //annihilate interaction
#ifdef DEBUGMODE
            cout << "\n  annihilate " << i3;
#endif
	    SiteType ib = worm_it->from();

	    // delete interaction on the from site
	    find_assoc_delete(worm_it->to(),   site_it[worm_it->to()]  );
	    find_assoc_delete(worm_it->from(), site_it[worm_it->from()]);
	    site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();
	    site_left(worm_head.to());

	    // delete interaction on the to site
            worm_head.to(ib);
            site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();
	    site_left(worm_head.to());
	    worm_it = site_it[worm_head.to()];
            p -= segm_tau * dE;
	    nrvertex--;
            return;
          }
          else { // relink the interaction
#ifdef DEBUGMODE
	    cout << "\n relink or bounce " << i3 << "\t" << nb(worm_it->from(),i3);
#endif
	    // define the new element
	    SiteType ifrom = worm_it->from();
            SiteType ib  = nb(ifrom,i3);

	    site_it[ifrom]->to(ib);

	    // delete the old element on the to site
	    find_assoc_delete(worm_it->to(), site_it[worm_it->to()]);
	    site_it[worm_head.to()] = operator_string[worm_head.to()].erase(site_it[worm_head.to()]);
            if (site_it[worm_head.to()] == operator_string[worm_head.to()].end()) site_it[worm_head.to()] = operator_string[worm_head.to()].begin();

	    // insert new element
	    int n0 = site_it[ib]->before();
	    element_type new_elem(n0,n0+1, ifrom, ib, worm_head.time() );
	    site_it[ib] = operator_string[ib].insert(site_it[ib], new_elem);
	    find_assoc_insert(ib, site_it[ib]);

	    worm_head.to(ib); // worm jumps
	    site_right(ib);
	    worm_it = site_it[ib];
	    worm_dir = -worm_dir;
	    worm_rising = !worm_rising;
            p -= segm_tau * dE;
            return;
          }
        } // ... else delete/relink interaction
      }
      else {// worm passes the interaction (no sites in common), but diagonal energy might change
#ifdef DEBUGMODE
        cout << "\nMoving downwards, worm passes remote interaction. This should not occur at all!!";
	char ch; cin >> ch;
#endif
      }
    } // shift_left_kink
    else { // no interaction reached, try to insert new interaction
      TimeType newtime = worm_head.time() - p/dE;
      while (newtime <= 0.) {
         newtime += beta;
      }
#ifdef DEBUGMODE
      cout << "\nMoving to the left, inserting interaction " << worm_head.to() << "\t" << worm_rising << "\t " << newtime;
#endif
      if (  (worm_tail.time() != worm_head.time() ) && (t_between(worm_tail.time(), newtime, worm_head.time() )) ) { 
	// measure density matrix
#ifdef DEBUGMODE
	cout << "\nMoving to the left, measuring density matrix ";
#endif

	worm_head.time(worm_tail.time());
	//dns[dist(worm_head.to(), worm_tail.to())] += start_bos;
	if (worm_head.to() != worm_tail.to()) {
          double cosphi=1.; for (int k=0; k < dim; k++) cosphi *= cos(phase[k]* (cdist(worm_head.to(),k) - cdist(worm_tail.to(),k)));
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
          av_dns[dist(worm_head.to(), worm_tail.to())] += 0.25*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.25*start_bos;
#else
          av_dns[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos;
#endif
#else
#ifdef TRAPPEDSYSTEM           
          av_dns[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 0.5*start_bos;
#else
          av_dns[dist(worm_head.to(), worm_tail.to())] += 1.*start_bos * cosphi;
	  av_dns_inf[dist(worm_head.to(), worm_tail.to())] += 1.*start_bos;
#endif 
#endif
	}

	p = 0.;
	return;
      }

      site_right(worm_head.to());   // then the iterator is the same as for right_insert
      worm_head.time(newtime);
      calc_ins_weight(dE, worm_head.to(), newtime); // this should also set the iterators right...

      SiteType i3 = heatbath_search(trans);

      if ( i3 < zcoord)
      {
        SiteType ib = nb(worm_head.to(),i3);
#ifdef DEBUGMODE
	cout << "\n ..." << i3 << "\t" << ib << "\t" << worm_head.to();
#endif
        element_type new_elem;
        new_elem.time(worm_head.time());
        if (worm_rising) {
	  new_elem.to(ib);
          new_elem.from(worm_head.to());
	  new_elem.after(site_it[worm_head.to()]->before() );
	  new_elem.before(site_it[worm_head.to()]->before() + 1  );
	  site_it[worm_head.to()] = operator_string[worm_head.to()].insert(site_it[worm_head.to()], new_elem);

	  new_elem.after(site_it[ib]->before());
	  new_elem.before(site_it[ib]->before()-1);
	  site_it[ib] = operator_string[ib].insert(site_it[ib], new_elem);
        }
        else {
          new_elem.to(worm_head.to());
	  new_elem.from(ib);
	  new_elem.after(  site_it[worm_head.to()]->before()  );
	  new_elem.before( site_it[worm_head.to()]->before() - 1    );
	  site_it[worm_head.to()] = operator_string[worm_head.to()].insert(site_it[worm_head.to()], new_elem);
	  
	  new_elem.after( site_it[ib]->before()  );
	  new_elem.before(site_it[ib]->before()+1);
	  site_it[ib] = operator_string[ib].insert(site_it[ib], new_elem);
	}
        find_assoc_insert(worm_head.to(), site_it[worm_head.to()]);
	find_assoc_insert(ib, site_it[ib]);
	site_left(worm_head.to());  //restore
	site_left(ib);
	
	// worm jumps
        worm_head.to(ib);
	worm_it = site_it[ib];
        p = 0.;
	nrvertex++;
	//print_conf(cout);
        return;
      }
      else // bounce back
      {
#ifdef DEBUGMODE
        cout << "\nInserting interaction failed, bouncing back, moving upwards";
#endif
	// site_right(worm_head.to()); is already done
	worm_it = site_it[worm_head.to()];
	worm_dir = -worm_dir;
	worm_rising = !worm_rising;
        p = 0.;
        return;
      }
    }
  } // worm_dir < 0  
}

/*
long lowa::calc_number_of_particles(double& d)
{
  long nr=0;
  d = 0.;
  for (SiteType j = 0; j < Nsites; j++) {
    int n = dummy_it[j]->before();
    nr += n; 
    state[j] = n;
    if (n==2) d+= 2.;
  }
  if (nr != 0) d /= (nr*1.);
  return (nr);
}
*/


TimeType lowa::calc_potential_energy()
{
  int32_t En1 = 0;
  TimeType En2 = 0.;

  number_of_particles=0;
  for (SiteType j = 0; j < Nsites; j++) {
    int n = dummy_it[j]->before();
    //En += 0.5*U_on*n*(n-1) + mu_eff[j]*n;
    En1 += U[j] *n*(n-1);
    En2 += mu_eff[j]*n;
    number_of_particles += n; 
    state[j] = n;
  }
  return (0.5*En1 + En2);
}



#ifdef BOUND_WORM
TimeType lowa::calc_M()               // impt for debug...
{
  TimeType dummy_sum = 0;
  for (SiteType j = 0; j < Nsites; j++) {
     int n = dummy_it[j]->before();
     dummy_sum += n;
  }
  dummy_sum *= beta;

  if (worm_diag == 1) 
     return dummy_sum;

  else {
     if (worm_dir == 1) 
        dummy_sum += (worm_rising ? (worm_tail.time() - worm_head.time()) : (worm_head.time() - worm_tail.time()));
     else 
        dummy_sum -= (worm_rising ? (worm_tail.time() - worm_head.time()) : (worm_head.time() - worm_tail.time()));

     if ((worm_dir == 1) && (worm_head.time() == beta))
        dummy_sum += (worm_rising ? beta : -beta);        
  
     return dummy_sum;
  }
}    
#endif

bool lowa::test_conf()
{ 
  if ((!worm_diag) && (worm_it != site_it[worm_head.to()])) {
    cout << "\n Error with site_it and worm_it";
    return(0);
  }
  for (SiteType i = 0; i < Nsites; i++) {
    int n = operator_string[i].size();
    int ii = 0;
    for (list<element_type>::iterator it = operator_string[i].begin(); it != operator_string[i].end(); ++it, ++ii) {
      if (it->before() < 0 || it->before() > nmax) {
	cout << "\n Error with density  << " << i << "\t "<< &*it << "\t" << it->time();
	return (0);
      }
      if (it->after() < 0 || it->after() > nmax) {
	cout << "\n Error with density  << " << i << "\t "<< &*it << "\t" << it->time();
	return (0);
      }
    }
  }
  return (1);
}

void lowa::calc_winding() {
/*
  for (int k = 0; k < dim; k++) mWinding[k] = 0.;
  for (SiteType s = 0; s < Nsites; s++) {
    list<element_type>::const_iterator it=operator_string[s].begin();
    if (operator_string[s].size() > 1) {
      int i;
      for (; it !=operator_string[s].end(); ++it) {
	if (it->from() != it->to()) {
	  for (int k = 0; k < dim; k++) {
	    mWinding[k] += winding_element(bond(it->from(), it->to()), k);
	  }
	} //if
      } // ...iterator over operator_string[site]
    } // if
  } // ...over all sites
  for (int k = 0; k < dim; k++) mWinding[k] /= (2.*Ls[k]);
*/
}


void lowa::save_internal(ostream& os, vector<int32_t>& counter)
{
  os << itherm << "\t" << iloop << endl;
  os << counter.size() << "\n";
  for (int i = 0; i < counter.size(); i++) os << counter[i] << "\t";
  os << "\n";
  os << MCstep << "\t" << MCstep_total << "\t" << MCstep_run << "\t" << MCmeas << "\n";
  os << worm_head << endl;
  os << worm_tail << endl;
  os << start_bos << endl;
  os << *worm_it << endl;
  for (SiteType s = 0; s < Nsites; s++) {
    int n = 0;
    os << operator_string[s].size() << endl;
    list<element_type>::iterator it = dummy_it[s];
    os << *it << endl;
    it->name(n); n++;
    ++it; if (it == operator_string[s].end()) it = operator_string[s].begin();
    while (it != dummy_it[s]) {
      os << *it << endl;
      it->name(n); n++;
      ++it; if (it == operator_string[s].end()) it = operator_string[s].begin();
    }
  }
  // save associations
  for (SiteType s = 0; s < Nsites; s++) {
    list<element_type>::iterator it = dummy_it[s];
    for (int j = 0; j < zcoord; j++) os << (it->get_assoc(j))->name() << "\t";
    os << endl;
    ++it; if (it == operator_string[s].end()) it = operator_string[s].begin();
    while (it != dummy_it[s]) {
      for (int j = 0; j < zcoord; j++) os << (it->get_assoc(j))->name() << "\t";
      os << endl;
      ++it; if (it == operator_string[s].end()) it = operator_string[s].begin();
    }
  }

/*
  for (SiteType s = 0; s < Nsites; s++) {
    os << av_state[s] << "\t" << av_state_sq[s] << "\t" << av_dns[s] << "\t" << av_dns_inf[s] << endl;
  }

  for (SiteType s = 0; s < projected_Nsites; s++) {
    os << av_state_projected[s] << "\t" << av_state_sq_projected[s] << endl;
  }

  os << mZ_state << "\t" << mZ_dns << endl;
*/

  for (SiteType s = 0; s < Nsites; s++) {
    os << av_dns[s] << "\t" << av_dns_inf[s] << endl;
  }
  os << mZ_dns << endl;


  os << worm_diag << "\t" << worm_dir << "\t" << worm_rising << endl;
  os << nrvertex << endl;
#ifdef BOUND_WORM
  os << M__ << endl;
  os << correct_lattice_structure << endl;
#endif
  os << Ekin << "\t" << Epot << "\t" << new_measurement << endl;
  os << acc_insert << "\t" << tot_insert << endl;
  os << times_run << "\t" << times_overall << endl;
  os << sweeps << endl;
  return;
}


void lowa::load_internal(istream& is, vector<int32_t >& counter)
{
  list<element_type>::iterator** v;
  v = new list<element_type>::iterator* [Nsites];
  element_type a_worm_elem;
  is >> itherm >> iloop;
  int tmp_dim;
  is >> tmp_dim;
  counter.resize(tmp_dim);
  for (int i = 0; i < counter.size(); i++)  is >> counter[i];
  is >> MCstep >> MCstep_total >> MCstep_run >> MCmeas;

  MCstep = 0.;  MCstep_total = 0. ; MCstep_run = 0.; MCmeas = 0.;

  is >> worm_head;
  is >> worm_tail;
  is >> start_bos;
  is >> a_worm_elem;
  for (SiteType s = 0; s < Nsites; s++) {
    int64_t opsize_n;
    is >> opsize_n;
    v[s] = new list<element_type>::iterator [opsize_n];
    for (int32_t i = 0; i < opsize_n; i++) {
      element_type elem;
      is >> elem;
      operator_string[s].push_back(elem);
      if (elem == a_worm_elem) { worm_it = operator_string[s].end(); --worm_it;}
    }
    int32_t i = 0;
    for (list<element_type>::iterator it=operator_string[s].begin(); it != operator_string[s].end(); ++it, i++) {
      v[s][i]=it;
    }
    dummy_it[s] = operator_string[s].begin();
    site_it[s] = dummy_it[s];
  }
  // load associations
  for (SiteType i = 0; i < Nsites; i++) {
    int ass;
    list<element_type>::iterator it = dummy_it[i];
    for (int j = 0; j < zcoord; j++) {
      SiteType s = nb(i,j);
      is >> ass;
      it->set_assoc(j, v[s][ass] );
    }
    ++it; if (it == operator_string[i].end()) it = operator_string[i].begin();
    while (it != dummy_it[i]) {
      for (int j = 0; j < zcoord; j++) {
	SiteType s = nb(i,j);
	is >> ass;
	it->set_assoc(j, v[s][ass] );
      }
      ++it; if (it == operator_string[i].end()) it = operator_string[i].begin();
    }
  } // ... associations

/*
  for (SiteType s = 0; s < Nsites; s++) {
    is >> av_state[s] >> av_state_sq[s] >> av_dns[s] >> av_dns_inf[s];
  }
  for (SiteType s = 0; s < projected_Nsites; s++) {
    is >> av_state_projected[s] >> av_state_sq_projected[s];
  }
  is >> mZ_state >> mZ_dns;
*/

  for (SiteType s = 0; s < Nsites; s++) {
    is >> av_dns[s] >> av_dns_inf[s];
  }
  is >> mZ_dns;

  is >> worm_diag >>  worm_dir >> worm_rising;
  is >> nrvertex;
#ifdef BOUND_WORM
  is >> M__;
  is >> correct_lattice_structure;
#endif
  is >> Ekin >> Epot >> new_measurement;
  is >> acc_insert >> tot_insert;
  is >> times_run >> times_overall;

  uint32_t dummy_sweeps;
  is >> dummy_sweeps;

  calc_potential_energy();
  calc_winding();
  for (int i = 0; i < Nsites; i++) {
    delete [] v[i];
  }
  delete [] v;
  site_it[worm_head.to()] = worm_it;
#ifdef DEBUGMODE
  cout << "\n dummy :\n ";
  for (SiteType i = 0; i < Nsites; i++) {
    cout << dummy_it[i]->time() << "\t";
  }
  cout << "\nWorm iterator : " << worm_it->time() << "\t" << &*worm_it << "\n";
  print_conf(cout);
#endif
}


void lowa::find_assoc_insert(SiteType cursite, list<element_type>::iterator it)
{
  // go one up on the operator string for the current element and look to what is pointing
  list<element_type>::iterator ito = it;
  ++ito;
  if (ito == operator_string[cursite].end()) ito = operator_string[cursite].begin();

  for (SiteType j = 0; j < zcoord; j++) {
    SiteType s = nb(cursite, j);
    list<element_type>::iterator itl;
    
    itl = ito->get_assoc(j);
    list<element_type>::iterator itp = itl;
    if (itp == operator_string[s].begin()) itp = operator_string[s].end();
    --itp;
    it->set_assoc(j, itl); 
    while (!t_between(it->time(), itp->time(), itl->time())) {
      itl = itp;
      if (itp == operator_string[s].begin()) itp = operator_string[s].end();
      --itp;
      it->set_assoc(j, itl); 
    }
  }

  // interactions on the neighbors below the time of it might also change...
  for (SiteType j = 0; j < zcoord; j++) {
    SiteType oppdir = (j + dim) % zcoord;
    SiteType s = nb(cursite, j);
    list<element_type>::iterator itl;
    itl = it->get_assoc(j);
    list<element_type>::iterator itp = itl;
    if (itp == operator_string[s].begin()) itp = operator_string[s].end();
    --itp;
    list<element_type>::iterator itw = itp->get_assoc(oppdir);
    if (itl->time() == it->time()) itl->set_assoc(oppdir, it);
    while ((itp->time() != itw->time())  && ( t_between(it->time(), itp->time(), itw->time()) )) {
      itp->set_assoc(oppdir, it);
      if (itp == operator_string[s].begin()) itp = operator_string[s].end();
      --itp;
      itw = itp->get_assoc(oppdir);
    }
  }
}

void lowa::find_assoc_delete(SiteType cursite, list<element_type>::iterator it)
{
  if (operator_string[cursite].size() == 1) {
     for (SiteType j = 0; j < zcoord; j++) {
      SiteType oppdir = (j + dim) % zcoord;
      SiteType s = nb(cursite, j);
	  for (list<element_type>::iterator itt=operator_string[s].begin(); itt != operator_string[s].end(); ++itt) 
	    itt->set_assoc(oppdir, operator_string[cursite].begin());
	  }
  }
  else {
    list<element_type>::iterator newpoint=it;
    ++newpoint;
    if (newpoint == operator_string[cursite].end()) newpoint = operator_string[cursite].begin();
    list<element_type>::iterator itp = it;
    if (itp == operator_string[cursite].begin()) itp=operator_string[cursite].end();
    --itp;
    // on nb sites, go down until you find an element that points at the current element
    for (SiteType j = 0; j < zcoord; j++) {
      SiteType oppdir = (j + dim) % zcoord;
      SiteType s = nb(cursite, j);
      list<element_type>::iterator it_end=it->get_assoc(j);
      list<element_type>::iterator it_begin=itp->get_assoc(j);
      list<element_type>::iterator itt=it_begin;
      if ((itt == it_end) && (itt == dummy_it[s])) {
	for (list<element_type>::iterator itt=operator_string[s].begin(); itt != operator_string[s].end(); ++itt) {
	  if (itt->get_assoc(oppdir) == it) itt->set_assoc(oppdir, newpoint);
	}
      }
      else {
	while (itt != it_end) {
	  if (itt->get_assoc(oppdir) == it) {
	    itt->set_assoc(oppdir, newpoint);
	  }
	  ++itt;
	  if (itt==operator_string[s].end()) itt=operator_string[s].begin();
	}
	if (it_end->get_assoc(oppdir) == it) {
	  it_end->set_assoc(oppdir, newpoint);
	}
      } //else
    }
  }
}

void lowa::print_conf(ostream& out)
{
  out << "\n\n Printing operator string";
  for (SiteType i = 0; i < Nsites; i++) {
    out << "\nSite : " << i;
    int n = operator_string[i].size();
    int ii = 0;
    for (list<element_type>::iterator it = operator_string[i].begin(); it != operator_string[i].end(); ++it, ++ii) {
      it->print();
      if (ii > n) {
	cout << "\n Error with list!\n";
	char ch; cin >> ch;
      }
    }
    out << "\n--------------------\n\n\n";
  }
  cout << "\n Worm site : " << worm_head.to() << "\t time " << worm_head.time() << "\t dir "<< worm_dir << "\t rising " << worm_rising << "\n";
  //char ch; cin >> ch;
}


void lowa::dostep() {
   // dostep is the core unit. It generates a new diagonal configuration
   int is_save = dostep_internal();
 
   if (is_save) {
      cout << "\n Saving...\t";
      MCstep_run   += MCstep*1.;
      MCstep_total += MCstep*1.;
      MCstep       =  0;
      time (&times2);
      dtimes = times2 - times1;
      times1 = times2;
      times_tot += dtimes;
      times_overall += dtimes;
      times_run += dtimes;

      update_av();

      ofstream outFile(filename1_trial.c_str(), ios::binary);
      outFile << setprecision(20);
      save_internal(outFile, counter);
      outFile.close();

      int bk1;
      bk1 = std::rename(filename1_trial.c_str(),filename1.c_str());

#ifdef PRINT_INFOMATION_WHILE_SAVING
/*
      ofstream outDens(filename2_trial.c_str());
      outDens << setprecision(20);
      print_av(outDens);
      outDens.close();

      int bk2;
      bk2 = std::rename(filename2_trial.c_str(),filename2.c_str());


      ofstream outDens_projected(filename3_trial.c_str());
      outDens_projected << setprecision(20);
      print_av_projected(outDens_projected);
      outDens_projected.close();

      int bk3;
      bk3 = std::rename(filename3_trial.c_str(),filename3.c_str());
*/
#endif

      cout << "\t...done\n";


/*
      if (times_tot + 1.5*dtimes > times_runtime ) {
         cout << "\nApproaching run time limit. Saving and Quitting....\n\n";
	 update_av();      
	 ofstream outDens(filename2.c_str());
	 outDens << setprecision(20);
	 print_av(outDens);
	 outDens.close();
	 ofstream outDensChem(filename0.c_str());
	 outDensChem << setprecision(20);
	 print_av_chem(outDensChem);
	 outDens.close();
	 cout << "# Worm insert ratio (run)     : " << get_worm_insert_ratio() << endl;
	 cout << "# MCsteps (overall)           : " << MCstep_total << endl;
	 cout << "# MCsteps (run)               : " << MCstep_run << endl;
	 cout << "# Nr MCstep per s (run)       : " << MCstep_run / times_run << endl;
	 double d = MCmeas / MCstep_run;
	 cout << "# ratio meas/ non-diag        : " << d << endl;
	 if ( ( d < 0.005) || (d > 0.1) ) {
	    cout << "# CHANGE FREQUENCY OF MEASURING to Nmeasure " << 0.01/d/Nmeasure << "\t" << endl;
	 }
         cout << "# CHEMICAL POTENTIAL = " << mu << " ; <N> = " << dummy_av_Npart << endl;
	 return(0);
      }
*/
   }
   else {
      counter[MEASURE]++;
      counter[MEASURE_GREEN]++;
      counter[TEST]++;

      if (counter[MEASURE] >= Nmeasure) { 

         if ((Ncan < 0) || (get_Npart() == Ncan)) {
            ++sweeps;

            MCmeas++;
            time (&timesm2);
            cout << "\nMeasuring..." << get_Npart() << "\t MCstep : " << MCstep_run+MCstep << "\t worm insert ratio : " << get_worm_insert_ratio() << " Average nr MCstep per second  : " << (MCstep_run+MCstep - MCold) / ( timesm2 - timesm1) << " Sweeps : " << sweeps;
            MCold = MCstep_run + MCstep;
            timesm1 = timesm2;
            update_av();      


            // particle number measurements
            measurements["N"]     << get_Npart();

            N_out.open(filenameN.c_str(),ios::app);
            N_out << get_Npart() << endl;
            N_out.close();

            // projected binned density measurements
            measurements["n"]    << projected_binned_state;
            measurements["nN"]   << projected_binned_state_times_N;


#ifndef TRAPPEDSYSTEM
            get_winding(this_winding);       // O(N) operation, expensive to measure this way
            for (int k = 0; k < dim; k++) {this_winding_squared[k] = this_winding[k] * this_winding[k];};
            measurements["Winding^2"] << this_winding_squared;
            measurements["Superfluid Density"] << get_rho_sf();
#endif
         }
         counter[MEASURE] = 0;
      }

      if (counter[TEST] >= Ntest) {
         cout << "\n Testing...\t";
         if (!test_conf()) {char ch; cin >> ch;}; 
         cout << "\t...OK\n";
         counter[TEST] = 0;
      }
   }
}


inline void lowa::calc_del_weight(const SiteType cur_site)
// interaction elements to delete an interaction
{
  // cur_site is here NOT the worm site, it is the other bond site
  std::vector<TimeType> weight;
  std::vector<int> nbs;
  std::vector<SiteType>::iterator it;
  std::vector<TimeType>::iterator iter;
  SiteType bond;

  weight.resize(zcoord+1);
  nbs.resize(zcoord);

  TimeType Ebis;

  SiteType i = 0;
  int cur_dens = (worm_dir == 1 ? site_it[cur_site]->before() : site_it[cur_site]->after());


  // find the densities on the neigboring sites
  // site_it[cur_site] should point to the element which we want to delete
  // all neighboring iterators should be put right as well
  for ( it = nbs.begin(); it != nbs.end(); ++it, ++i) {
    SiteType s = nb(cur_site,i);
    site_it[s] = site_it[cur_site]->get_assoc(i);
    *it = site_it[s]->before();
    if (s == worm_head.to()) {
      *it = (worm_dir == 1 ? site_it[worm_head.to()]->after() : site_it[worm_head.to()]->before() );
      // when moving upwards, we need the density just below the worm; which is the same as after the interaction.

    }
    //cout << "\n calc_del_w : " << cur_site << "\t "<< &*site_it[cur_site] << "\t" << s << "\t" << site_it[s]->time();
  }

  i = 0;
  if (worm_rising)
  { // only work for 3 dimensions...
    weight[0] = (nbs[0] == nmax ? 0 : (t_x_plus[cur_site] * (nbs[0]+1)));
    weight[1] = (nbs[1] == nmax ? 0 : (t_y_plus[cur_site] * (nbs[1]+1)));
    weight[2] = (nbs[2] == nmax ? 0 : (t_z_plus[cur_site] * (nbs[2]+1)));
    weight[3] = (nbs[3] == nmax ? 0 : (t_x_minus[cur_site] * (nbs[3]+1)));
    weight[4] = (nbs[4] == nmax ? 0 : (t_y_minus[cur_site] * (nbs[4]+1)));
    weight[5] = (nbs[5] == nmax ? 0 : (t_z_minus[cur_site] * (nbs[5]+1)));

    Ebis = U[cur_site] * cur_dens + mu_eff[cur_site];
    weight[6] = (Ebis > 0 ? (E_off+Ebis) : E_off);
  }
  else
  {
    weight[0] = t_x_plus[cur_site] * nbs[0];
    weight[1] = t_y_plus[cur_site] * nbs[1];
    weight[2] = t_z_plus[cur_site] * nbs[2];
    weight[3] = t_x_minus[cur_site] * nbs[3];
    weight[4] = t_y_minus[cur_site] * nbs[4];
    weight[5] = t_z_minus[cur_site] * nbs[5];

    Ebis = -U[cur_site] * (cur_dens-1) - mu_eff[cur_site];
    weight[zcoord] = (Ebis > 0 ? (E_off + Ebis) : E_off);
  }

  // heatbath weights
  // this seems enough... please change to locally optimal ones if you are not happy with autocorr
  TimeType sw = accumulate(weight.begin(), weight.end(), 0.);
  trans[0] = weight[0] / sw;
  i=1;
  iter = weight.begin();
  for (++iter; iter != weight.end(); ++iter, ++i)
  {
    trans[i] = weight[i] / sw + trans[i-1];
  }
}


inline void lowa::calc_ins_weight(const TimeType E, const SiteType cur_site, const TimeType instime)
// matrix elements to insert an interaction
{
  std::vector<TimeType> weight(zcoord+1);
  std::vector<int> nbs(zcoord);
  std::vector<SiteType> :: iterator it;

  // density on the current state is irrelevant since we know E?
  weight[6] = E; //  when bouncing
  SiteType i = 0;
  for ( it = nbs.begin(); it != nbs.end(); ++it, ++i) {
    SiteType s = nb(cur_site,i);
    site_it[s] = site_it[cur_site]->get_assoc(i);
    list<element_type>::iterator previt = site_it[s];
    //cout << "\n calc_ins_weight : " << cur_site << "\t" << &*site_it[cur_site] << "\t" << s << *site_it[s];

    if (previt == operator_string[s].begin()) previt = operator_string[s].end();
    --previt;
    while (!t_between(instime, previt->time(), site_it[s]->time() ) ) {
      site_it[s] = previt;
      if (previt == operator_string[s].begin()) previt = operator_string[s].end();
      --previt;
    }
    *it = site_it[s]->before();
  }

  if (worm_rising)
  {
    weight[0] = t_x_plus[cur_site] * nbs[0];
    weight[1] = t_y_plus[cur_site] * nbs[1];
    weight[2] = t_z_plus[cur_site] * nbs[2];
    weight[3] = t_x_minus[cur_site] * nbs[3];
    weight[4] = t_y_minus[cur_site] * nbs[4];
    weight[5] = t_z_minus[cur_site] * nbs[5];
  }
  else
  {
    weight[0] = (nbs[0] == nmax ? 0 : (t_x_plus[cur_site] * (nbs[0]+1)));
    weight[1] = (nbs[1] == nmax ? 0 : (t_y_plus[cur_site] * (nbs[1]+1)));
    weight[2] = (nbs[2] == nmax ? 0 : (t_z_plus[cur_site] * (nbs[2]+1)));
    weight[3] = (nbs[3] == nmax ? 0 : (t_x_minus[cur_site] * (nbs[3]+1)));
    weight[4] = (nbs[4] == nmax ? 0 : (t_y_minus[cur_site] * (nbs[4]+1)));
    weight[5] = (nbs[5] == nmax ? 0 : (t_z_minus[cur_site] * (nbs[5]+1)));
  }
  TimeType sw = accumulate(weight.begin(), weight.end(), 0.);


  //heatbath weights
  i=1;
  trans[0] = weight[0] / sw;
  std::vector<TimeType>::const_iterator iter = weight.begin();
  for (++iter; iter != weight.end(); ++iter, ++i)
  {
    trans[i] = weight[i] / sw + trans[i-1];
  }
}


void lowa::load(alps::IDump& dump)   {}
void lowa::save(alps::ODump& dump) const   {}


bool lowa::change_parameter(const std::string& name,const alps::StringValue& value)
{
  if(name=="SWEEPS")
    total_sweeps=static_cast<uint32_t>(value);
  else if (name=="THERMALIZATION" && !is_thermalized())
    thermalization_sweeps=static_cast<uint32_t>(value);
  else
    return false; // cannot do it
  return true; // can do it
}

bool lowa::is_thermalized() const
{
  return (sweeps >= thermalization_sweeps);
}

double lowa::work_done() const 
{
  return (is_thermalized() ? (sweeps-thermalization_sweeps)/double(total_sweeps) : 0.);
}
