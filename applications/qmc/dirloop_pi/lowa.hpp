/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2010 by Lode Pollet <pollet@itp.phys.ethz.ch>
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


#ifndef lowa_HPP
#define lowa_HPP

#include <alps/scheduler/montecarlo.h>
#include "lowa.element.hpp"
#include <valarray>
#include <numeric>
#include <time.h>
#include <iomanip>
#include <sstream>
#include <vector>
#include <list>

// dim and zcoord are defined in lowa.Element.hpp
//#define DEBUGMODE
//#define HARDCORE
#define TRAPPEDSYSTEM      // important for distance
#define PRINT_INFOMATION_WHILE_SAVING


const double tol = 1e-10;
const long Nsave = 100000000;
const double hbar = 1.05457148;
const double amu  = 1.66053886;



using namespace std;


class lowa : public alps::scheduler::MCRun
{
 public:
  enum {MEASURE, MEASURE_GREEN, TEST};
  vector<int32_t> counter;
  uint32_t label;
  double MCstep_total, MCstep_run;
  long MCstep;
  long MCold;
  int32_t MCmeas;
  double acc_insert, tot_insert;
  typedef Element element_type;

  lowa(const alps::ProcessList&,const alps::Parameters&,int);
/*
  ~lowa() {
    if (dim == 1) {
      delete Ls;
      delete mWinding;
      delete vc;
      delete v0;
      delete mid;
      delete phase;
      delete waist;
      delete lambda;
    }
    else {
      delete [] Ls;
      delete [] mWinding;
      delete [] vc;
      delete [] v0;
      delete [] mid;
      delete [] phase;
      delete [] waist;
      delete [] lambda;
    }
    delete [] operator_string;
    delete [] dummy_it;
    delete [] site_it;
    delete [] t_x_plus;
    delete [] t_x_minus;
    delete [] t_y_plus;
    delete [] t_y_minus;
    delete [] t_z_plus;
    delete [] t_z_minus;
    delete [] U;
    delete [] mu_eff;
    //delete [] dns, av_dns_inf, av_state, av_state_sq;
    delete [] av_dns, av_dns_inf;
    //delete [] av_state, av_state_sq;
    //delete [] av_state_projected, av_state_sq_projected;
  }
*/
  static void print_copyright(std::ostream&);

  void init();
  void init_internal();
  void dostep();
  int dostep_internal();
  void read_params();
  void print_params() const;

  void change_CP(double);

  SiteType dist(const SiteType s1, const SiteType s2);
  //double dist_realsq(const SiteType s1);
  TimeType get_kinetic_energy() {if (new_measurement) update_en(); return (Ekin);} // -nrvertex/beta);}
  TimeType get_potential_energy() {if (new_measurement) update_en(); return (Epot);} //{return (calc_potential_energy());}
  TimeType get_energy() {return (Epot + Ekin);} //{return (calc_potential_energy() - nrvertex/beta);}
  double get_Npart() { if (new_measurement) update_en();  return (number_of_particles); }// {double d; return(calc_number_of_particles(d));}
  long get_Npart(double& d) { if (new_measurement) update_en(); return (number_of_particles); } //{return(calc_number_of_particles(d));}
  double get_worm_insert_ratio() {return (acc_insert / tot_insert);}
  TimeType calc_potential_energy();
  void calc_winding();
  //long calc_number_of_particles(double& );
  TimeType calc_local_energy(const SiteType) const;
  void geometry();
  void update_system(const SiteType, int&);
  void worm_cycle(TimeType&);
  bool worm_diag;

  void find_assoc_insert(SiteType, list<element_type>::iterator);
  void find_assoc_delete(SiteType, list<element_type>::iterator);
  bool shift_right_kink(TimeType&, double, TimeType&);
  bool shift_left_kink(TimeType&, double, TimeType&);
  void site_right(const SiteType s);
  void site_left(const SiteType s);
  bool t_between(const TimeType, const TimeType, const TimeType) const;
  void worm_left();
  void worm_right();
  void worm_pass_next();
  void worm_pass_prev();
  void worm_changedir_toleft();
  void worm_changedir_toright();
  void calc_del_weight(const SiteType);
  void calc_ins_weight(const TimeType , const SiteType, const TimeType);
  SiteType calc_dist(const SiteType, const SiteType);

  void print_conf(ostream&);
  void save_internal(ostream&, vector <int32_t >& );
  void load_internal(istream&, vector <int32_t >& );

  void save(alps::ODump&) const;
  void load(alps::IDump&);
  bool is_thermalized() const;
  double work_done() const;
  bool change_parameter(const std::string& name,const alps::StringValue& value);

  bool test_conf();
  SiteType heatbath_search(std::vector<TimeType>&);

  uint32_t itherm, iloop;


  void update_av() {
/*
    mZ_dns += 1.;
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
    for (SiteType i = 1; i < Nsites; i++) av_dns_inf[i] += 0.25*dns[i];
    av_dns_inf[0] += dns[0];
    //for (SiteType i = 0; i < Nsites; i++) m[0] = dns[0];
#else
    for (SiteType i = 1; i < Nsites; i++) av_dns_inf[i] += 0.5*dns[i];
    av_dns_inf[0] += dns[0];
    //for (SiteType i = 0; i < Nsites; i++) m[0] = dns[0];
#endif
#else
#ifdef TRAPPEDSYSTEM
    for (SiteType i = 1; i < Nsites; i++) av_dns_inf[i] += 0.5*dns[i];
    av_dns_inf[0] += dns[0];
#else
    for (SiteType i = 0; i < Nsites; i++) av_dns_inf[i] += dns[i];
#endif
#endif
*/

/*
    mZ_state += 1.;
    for (SiteType i = 0; i < Nsites; i++) av_state[i] += state[i];
    for (SiteType i = 0; i < Nsites; i++) av_state_sq[i] += state[i]*state[i];
    if (new_measurement) update_en();

    for (SiteType i = 0; i < projected_Nsites; i++) {
       TimeType dummy_projected_state = 0.;
       for (SiteType dummy = 0; dummy < Ls[2]; dummy++) {
          SiteType current_dummy_site =  i + (dummy * projected_Nsites); 
          dummy_projected_state       += state[current_dummy_site]; 
       }
       dummy_projected_state /= Ls[2];
       av_state_projected[i]    += dummy_projected_state;
       av_state_sq_projected[i] += dummy_projected_state * dummy_projected_state;
    }
*/

    if (new_measurement) update_en();
    to_projected_binned_state_from_state();
  }


  void update_en() {
    Ekin = -nrvertex / beta;
    Epot = calc_potential_energy();
    new_measurement = false;
  }


  void reset_av() {
    mZ_dns = 0.;
    //mZ_state = 0.;
    for (int32_t i = 0; i < Nsites; i++) {
      av_dns[i]           = 0.;
      av_dns_inf[i]       = 0.;
      //av_state[i]         = 0.;
      //av_state_sq[i]      = 0.;
    }
    /*
    for (int32_t i=0; i < projected_Nsites; i++) {
      av_state_projected[i]    = 0.;
      av_state_sq_projected[i] = 0.;
    }
    */
    /*
    tot_insert = 0.;
    acc_insert = 0.;
    MCstep = 0;
    MCstep_run = 0.;
    MCmeas = 0;
    times_run = 0;
    */
  }


/*
  void print_av(std::ostream& os)   {
    if (dim == 3) {
      for (int k = 0; k < Ls[2]; k++) {
        for (int j = 0; j < Ls[1]; j++) {
          for (int i = 0; i < Ls[0]; i++) {
            int32_t s = i + j*Ls[0] + k * Ls[0]*Ls[1];
            os << i << "\t" << j << "\t" << k << "\t" << av_state[s]/mZ_state << "\t" << av_state_sq[s]/mZ_state << "\t" << av_dns[s]/mZ_dns << "\t" << av_dns_inf[s]/mZ_dns << "\n" ;
          }
        }
      }
    }
  }
*/

/*
  void print_av_projected(std::ostream& os)   {
    if (dim == 3) {
      for (int j = 0; j < Ls[1]; j++) {
        for (int i = 0; i < Ls[0]; i++) {
          int32_t s = i + j*Ls[0];
          os << i << "\t" << j << "\t" << av_state_projected[s]/mZ_state << "\t" << av_state_sq_projected[s]/mZ_state << "\n";
        }
      }
    }
  }
*/

/*
  void print_av_chem(std::ostream& os) {
    if (dim == 3) {
      //int32_t Npoints = min( min (Ls[0], Ls[1]), Ls[2]);
      //Npoints = (Npoints % 2 == 0 ? Npoints / 2 : (Npoints + 1) / 2);
      int32_t vol = 1;
      for (int k = 0; k < dim; k++) vol += Ls[k]*Ls[k]/4;
      int32_t Npoints = int(sqrt(vol*1.))+1;
      vector<int32_t> Nentry(Npoints);
      vector<double> prof(Npoints);
      for (int32_t i = 0; i < Npoints; i++) {
        Nentry[i] = 0;
        prof[i] = 0.;
      }

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
            d = sqrt(d/vc[1]);
            int32_t pos = int(d);
            if (pos >= Npoints) {
              os << "Increase Npoints " << pos << "\t" << Npoints << "\n";
              pos = Npoints-1;
            }
            Nentry[pos]++;
            prof[pos] += av_state[p];
          }
        }
      }

      for (int32_t i = 0; i < Npoints; i++) {
        if (Nentry[i] > 0) os << i <<  "\t" << prof[i]/mZ_state/Nentry[i] << "\t" << Nentry[i] << "\n";
      }

    }
    else {
      os << "print_av_chem only implemented for dim == 3.\n";
    }
  }
*/

  void get_green(TimeType* m)   {
/*
#ifdef HARDCORE
#ifdef TRAPPEDSYSTEM
    for (SiteType i = 0; i < Nsites; i++) m[i] = 0.25*dns[i];
    for (SiteType i = 0; i < Nsites; i++) m[0] = dns[0];
#else
    for (SiteType i = 0; i < Nsites; i++) m[i] = 0.5*dns[i];
    for (SiteType i = 0; i < Nsites; i++) m[0] = dns[0];
#endif
#else
#ifdef TRAPPEDSYSTEM           // corrected after comparison with Nikolai was done for the 16x16x16 sample; I did not check it for HARDCORE yet
    for (SiteType i = 1; i < Nsites; i++) m[i] = 0.5*dns[i];
    m[0] = dns[0];
#else
    for (SiteType i = 0; i < Nsites; i++) m[i] = dns[i];
#endif
#endif
*/
  }


  void get_density(TimeType* m) {for (SiteType i = 0; i < Nsites; i++) m[i] = state[i];}
  //void get_winding(TimeType* m)  {calc_winding(); for (int i =0; i < dim; i++) m[i] = mWinding[i]/(Ls[i]*1.);}
  void get_winding(TimeType* m)  {calc_winding(); for (int i =0; i < dim; i++) m[i] = mWinding[i];}
  TimeType get_rho_sf() {
    // only call this function after having called get_winding or calc_winding;
    TimeType r = 0.;
    for (int i =0; i < dim; i++) {
      //r += (mWinding[i]/(Ls[i]*1.)) * (mWinding[i]/(Ls[i]*1.));
      r+= mWinding[i] * mWinding[i];
    }
    r /= (2. * beta);
    return (r);
  }

  std::string obtain_filename(std::string original_filename) {
     std::ostringstream ss;
     ss << label;
     std::string new_filename = original_filename + ss.str();
     return new_filename;
  }

  void to_projected_binned_state_from_state(void)
  {
    for (SiteType i=0; i < total_no_of_bins; ++i)   
      projected_binned_state[i] = 0.;
    
    for (SiteType j=0; j < Ls[1]; ++j) {
      SiteType p = j * Ls[0];
      for (SiteType i=0; i < Ls[0]; ++i) {
        SiteType cur_projected_site = p + i;
        SiteType cur_projected_bin  = projected_bin_no[cur_projected_site];

        for (SiteType k=0; k < Ls[2]; ++k) {
          projected_binned_state[cur_projected_bin] += state[cur_projected_site + k*projected_Nsites];
        }
      }
    }

    for (SiteType i=0; i < total_no_of_bins; ++i)
      projected_binned_state[i] /= bin_freq[i];

    for (SiteType i=0; i < total_no_of_bins; ++i)
      projected_binned_state_times_N[i] = projected_binned_state[i] * get_Npart();
  } 

  time_t times_tot, times1, times2, dtimes, times_overall, times_run, times_log1, times_log2, timesm1, timesm2;
  int read_configuration, read_statistics;
  int32_t Ntest, Nmeasure, Nmeasure_green;

private :
  // system parameters
  uint32_t centre_site;

  double t_hop;
  double U_on;
  double* vc;
  double mu;
  double beta;
  double time_of_flight;
  double* v0;
  double* waist;
  double a_s;
  double mass;
  double* lambda, *phase, *mid;
  SiteType* Ls;    // lattice dimensions
  SiteType Nsites; // total number of sites
  SiteType projected_Nsites;
  SiteType nmax;   // maximum occupation number
  double E_off;    // energy offset, internal parameter of the code
  long nrvertex;
  int32_t Ncan;       // canonical number of particles
  
  double* t_x_plus;
  double* t_x_minus;
  double* t_y_plus;
  double* t_y_minus;
  double* t_z_plus;
  double* t_z_minus;
  double* U;
  double* epsilon;

  double* mu_eff;

  long number_of_particles;

  // more intrinsic variables
  std::string filename0, filename1, filename2, filename3;
  std::string filename1_trial, filename2_trial, filename3_trial;
  std::string filenameN;
  std::ofstream N_out;

#ifndef TRAPPEDSYSTEM
  double *this_winding, *this_winding_squared;
#endif

  double Ekin, Epot;
  bool new_measurement;

  // system variables
  SiteType* state; // occupation per site
  list<element_type>* operator_string;
  TimeType En; // energy
  std::vector<TimeType> trans; // transition probabilies
  Matrix<SiteType, SiteType> nb;  // Everybody needs good neighbours

  // binned system variables
  TimeType bin_size;
  SiteType total_no_of_bins;
  SiteType* projected_bin_no;
  std::vector<SiteType> bin_freq;
  std::valarray<TimeType> projected_binned_state;
  std::valarray<TimeType> projected_binned_state_times_N;

  // worm variables
  list<element_type>::iterator worm_it;  // iterator pointing at the next element the worm will encounter
  list<element_type>::iterator*  site_it;  // iterator for site
  list<element_type>::iterator*  dummy_it; // iterator pointing at dummy elements useful for measuring diag prop

  element_type worm_head;                     // in Russian called ira
  element_type worm_tail;                     // in Russian called masha
  int worm_dir;                               // direction of the worm, to higher of lower imag times
  int worm_dir_init;
  bool worm_rising;                           // convention : density is increased or decreased between tail and head
  bool worm_rising_init;
  int start_bos;

  // measurement variables
  //TimeType* dns;
  Matrix<uint32_t, double> cdist;

  TimeType* av_dns;
  TimeType* av_dns_inf;
/*
  TimeType* av_state;
  TimeType* av_state_sq;
  TimeType* av_state_projected;
  TimeType* av_state_sq_projected;
  TimeType mZ_state;
*/
  TimeType mZ_dns;

  TimeType* mWinding;
  Matrix<int, TimeType> winding_element;
  Matrix<SiteType, SiteType> bond;

  uint64_t sweeps;
  uint64_t thermalization_sweeps;
  uint64_t total_sweeps;
};

inline bool lowa::shift_right_kink(TimeType& p, double dE, TimeType& t2)
{
  // see if the next interaction can be reached when moving to the right
  t2 = worm_it->time() - worm_head.time();          // t2 = time between worm head and the next interaction
  if (t2 <= 0.) t2 += beta;
  if (dE == 0.) return (true);
  TimeType t1 = p/dE;
  while (abs(t1 - t2) < tol) {
        p = -log(random_real());
        t1 = p/dE;
  }                                                 // p is automatically chosen again if it is not suitable, ie. t1 = p/dE
  if (t1 > t2) return (true);
  return (false);
}

inline bool lowa::shift_left_kink(TimeType& p, double dE, TimeType& t2)
{
  // see if the next interaction can be reached when moving to the left
  t2 = worm_head.time() - worm_it->time();
  if (t2 <= 0.) t2 += beta;
  if (dE == 0.) return (true);
  TimeType t1 = p/dE;
  while (abs(t1 - t2) < tol) {
        p = -log(random_real());
        t1 = p/dE;
  }
  if (t1 > t2) return (true);   // take care when t1 and t2 get close...
  return (false);
}


inline bool lowa::t_between(const TimeType t0, const TimeType t1, const TimeType t2) const
// check if t0 is between t1 and t2, given that t1 happens before t2 (mod beta)
{
  if (t2 > t1)
  {
    if ((t2>=t0) && (t0 > t1)) return (true);
  }
  else
  {
    if ((t0 > t1) || (t0 <= t2)) return (true);
  }
  return (false);
}

inline SiteType lowa::heatbath_search(std::vector<TimeType>& w)
// w is cumulative weight, last element = 1
{
  TimeType q = random_real();
  SiteType i = 0;
  for (std::vector<TimeType>::iterator it = w.begin(); it != w.end(); ++it, ++i)
  {
    if (q < *it) return (i);
  }
  std::cout << "\nError in heatbath_search.";
}

inline void lowa::site_right(const SiteType s)
  // move iterator to the right on site s
{
  ++(site_it[s]);
  if (site_it[s] == operator_string[s].end()) site_it[s] = operator_string[s].begin();
}

inline void lowa::site_left(const SiteType s)
  // move iterator to the left on site s
{
  if (site_it[s] == operator_string[s].begin()) site_it[s] = operator_string[s].end();
  --site_it[s];
}

#ifndef TRAPPEDSYSTEM
inline SiteType lowa::dist(const SiteType s1, const SiteType s2)
  // LATTICE requirements. returns the site that behaves to 0 as s2 behaves to s1
  // This is different for open and periodic boundary conditions (PBC)
{
  if (dim == 1) {
    SiteType s0;
    s0 = (s2 - s1 < 0 ? s2 - s1 + Ls[0] : s2 - s1);
    return s0;
  }
  else if (dim ==2) {
    SiteType i1, j1;
    SiteType i2, j2;
    SiteType i0, j0, s0;
    i1 = s1 % Ls[0];
    j1 = (s1 - i1) / Ls[0];
#ifdef DEBUGMODE
    if (s1 != i1 + j1*Ls[0]) {
      cout << "\n Error in dist1 " << s1 << "\t" << i1 << "\t" << j1;
      char ch; cin >> ch;
    }
    //cout << "\ndist : " << s1 << "\t" << i1 << "\t" << j1;
#endif
    i2 = s2 % Ls[0];
    j2 = (s2 - i2) / Ls[0];
#ifdef DEBUGMODE
    if (s2 != i2 + j2*Ls[0]) {
      cout << "\n Error in dist2 " << s2 << "\t" << i2 << "\t" << j2;
      char ch; cin >> ch;
    }
    //cout << "\ndist : " << s2 << "\t" << i2 << "\t" << j2;
#endif
    i0 = (i2-i1 >= 0 ? i2 - i1 : i2 -i1 + Ls[0]);
    j0 = (j2-j1 >= 0 ? j2 - j1 : j2 -j1 + Ls[1]);
    s0 = i0 + j0*Ls[0];
    //cout << "\ndist : " << s0 << "\t" << i0 << "\t" << j0;
    //char ch; cin >> ch;
    return (s0);
  }
  else if (dim == 3) {
    SiteType i1, j1, k1, t1;
    SiteType i2, j2, k2, t2;
    SiteType i0, j0, k0;
    SiteType s0;
    t1 = s1 % (Ls[0] * Ls[1]);
    k1 = (s1 - t1) / (Ls[0]*Ls[1]);
    i1 = t1 % Ls[0];
    j1 = (t1 - i1) / Ls[0];
#ifdef DEBUGMODE
    if (s1 != i1 + j1*Ls[0] + k1*(Ls[0]*Ls[1])) {
      cout << "\n Error in dist " << s1 << "\t" << i1 << "\t" << j1 << "\t" << k1;
      char ch; cin >> ch;
    }
#endif
    t2 = s2 % (Ls[0] * Ls[1]);
    k2 = (s2 - t2) / (Ls[0]*Ls[1]);
    i2 = t2 % Ls[0];
    j2 = (t2 - i2) / Ls[0];
#ifdef DEBUGMODE
    if (s1 != i1 + j1*Ls[0] + k1*(Ls[0]*Ls[1])) {
      cout << "\n Error in dist " << s1 << "\t" << i1 << "\t" << j1 << "\t" << k1;
      char ch; cin >> ch;
    }
#endif
    i0 = (i2-i1 >= 0 ? i2 - i1 : i2 -i1 + Ls[0]);
    j0 = (j2-j1 >= 0 ? j2 - j1 : j2 -j1 + Ls[1]);
    k0 = (k2-k1 >= 0 ? k2 - k1 : k2 -k1 + Ls[2]);
    s0 = i0 + j0*Ls[0] + k0*(Ls[0]*Ls[1]);
    return (s0);
  }
}
#else
inline SiteType lowa::dist(const SiteType s1, const SiteType s2)
{
  if (dim == 1) {
    SiteType s0;
    s0 = abs(s2 - s1);
    return s0;
  }
  else if (dim ==2) {
    SiteType i1, j1;
    SiteType i2, j2;
    SiteType i0, j0, s0;
    i1 = s1 % Ls[0];
    j1 = (s1 - i1) / Ls[0];
#ifdef DEBUGMODE
    if (s1 != i1 + j1*Ls[0]) {
      cout << "\n Error in dist1 " << s1 << "\t" << i1 << "\t" << j1;
      char ch; cin >> ch;
    }
    //cout << "\ndist : " << s1 << "\t" << i1 << "\t" << j1;
#endif
    i2 = s2 % Ls[0];
    j2 = (s2 - i2) / Ls[0];
#ifdef DEBUGMODE
    if (s2 != i2 + j2*Ls[0]) {
      cout << "\n Error in dist2 " << s2 << "\t" << i2 << "\t" << j2;
      char ch; cin >> ch;
    }
    //cout << "\ndist : " << s2 << "\t" << i2 << "\t" << j2;
#endif
    i0 = abs(i2-i1);
    j0 = abs(j2-j1);
    s0 = i0 + j0*Ls[0];
    //cout << "\ndist : " << s0 << "\t" << i0 << "\t" << j0;
    //char ch; cin >> ch;
    return (s0);
  }
  else if (dim == 3) {
    SiteType i1, j1, k1, t1;
    SiteType i2, j2, k2, t2;
    SiteType i0, j0, k0;
    SiteType s0;
    t1 = s1 % (Ls[0] * Ls[1]);
    k1 = (s1 - t1) / (Ls[0]*Ls[1]);
    i1 = t1 % Ls[0];
    j1 = (t1 - i1) / Ls[0];
#ifdef DEBUGMODE
    if (s1 != i1 + j1*Ls[0] + k1*(Ls[0]*Ls[1])) {
      cout << "\n Error in dist " << s1 << "\t" << i1 << "\t" << j1 << "\t" << k1;
      char ch; cin >> ch;
    }
#endif
    t2 = s2 % (Ls[0] * Ls[1]);
    k2 = (s2 - t2) / (Ls[0]*Ls[1]);
    i2 = t2 % Ls[0];
    j2 = (t2 - i2) / Ls[0];
#ifdef DEBUGMODE
    if (s1 != i1 + j1*Ls[0] + k1*(Ls[0]*Ls[1])) {
      cout << "\n Error in dist " << s1 << "\t" << i1 << "\t" << j1 << "\t" << k1;
      char ch; cin >> ch;
    }
#endif
    i0 = abs(i2-i1);
    j0 = abs(j2-j1);
    k0 = abs(k2-k1);
    s0 = i0 + j0*Ls[0] + k0*(Ls[0]*Ls[1]);
    return (s0);
  }
}
#endif


typedef alps::scheduler::SimpleMCFactory<lowa> lowa_Factory;
#endif
