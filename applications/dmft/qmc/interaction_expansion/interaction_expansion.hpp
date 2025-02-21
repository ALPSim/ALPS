/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef DMFT_QMC_WEAK_COUPLING_H
#define DMFT_QMC_WEAK_COUPLING_H

#include <iostream>
#include <vector>
#include <fstream>
#include <complex>
#include <cassert>
#include <cmath>
#include <string.h>
#include <alps/numeric/matrix.hpp>
#include "types.h"
#include "green_function.h"
#include "solver.h"
#include "alps_solver.h"
#include "fouriertransform.h"
#include "U_matrix.h"
#include "operator.hpp"
#include "green_matrix.hpp"


/*types*/
class c_or_cdagger;
class vertex;
typedef class histogram simple_hist;
typedef std::vector<vertex> vertex_array;



enum measurement_methods {
  selfenergy_measurement_matsubara, //measurement using self energy method
  selfenergy_measurement_itime_rs, //measurement using self energy method in imag time, real space
};



class histogram
{
public:
  
  histogram(unsigned int N):hist_(N, 0){}
  
  unsigned long &operator[](unsigned int n){return hist_[n];}
  const unsigned long &operator[](unsigned int n) const{return hist_[n];}
  unsigned int size() const{return hist_.size();}
  
  unsigned int max_index() const
  { 
    unsigned int max_index=0; 
    double max=0; 
    for(unsigned int i=0;i<hist_.size();++i){
      if(max<hist_[i]){
        max=hist_[i];
        max_index=i;
      }
    }
    return max_index;
  }
  
  unsigned int top_index() const
  { 
    unsigned int top_index=0;
    for(unsigned int i=0;i<hist_.size();++i){
      if(hist_[i]!=0){
        top_index=i;
      }
    }
    return top_index;
  }
  
  double max(const unsigned int index) const
  { 
    double max=0; 
    for(unsigned int i=0;i<index;++i){
      if(max<hist_[i]){
        max=hist_[i];
      }
    }
    return max;
  }
  
  double average(const unsigned int index) const{ 
    double average=0; 
    for(unsigned int i=0;i<index;++i){
      average+=hist_[i];
    }
    return average/index;
  }
  
  bool is_flat(const unsigned int index)const{return max(index)*0.8<average(index);}
  
  void clear()
  {
    for(unsigned int i=0;i<hist_.size();++i){
      hist_[i]=0;
    }
  }

private:
  
  std::vector<unsigned long> hist_;
};



typedef class vertex
{        
public:
  vertex(const spin_t &flavor1, const site_t &site1, const unsigned int &c_dagger_1, const unsigned int &c_1, 
         const spin_t &flavor2, const site_t &site2, const unsigned int &c_dagger_2, const unsigned int &c_2, 
         const double &abs_w)
  {
    z1_=flavor1;
    z2_=flavor2;
    s1_=site1;
    s2_=site2;
    c1dagger_=c_dagger_1;
    c2dagger_=c_dagger_2;
    c1_=c_1;
    c2_=c_2;
    abs_w_=abs_w;
  }
  
  inline const double &abs_w() const {return abs_w_;}
  inline const unsigned int &flavor1() const {return z1_;}
  inline const unsigned int &flavor2() const {return z2_;}
  inline const unsigned int &site1() const {return s1_;}
  inline const unsigned int &site2() const {return s2_;}
  inline void set_site1(site_t site1) {s1_=site1;}
  inline void set_site2(site_t site2) {s2_=site2;}
  inline const unsigned int &c_dagger_1() const {return c1dagger_;}
  inline const unsigned int &c_dagger_2() const {return c2dagger_;}
  inline const unsigned int &c_1() const {return c1_;}
  inline const unsigned int &c_2() const {return c2_;}
  inline unsigned int &flavor1() {return z1_;}
  inline unsigned int &flavor2() {return z2_;}
  inline unsigned int &c_dagger_1() {return c1dagger_;}
  inline unsigned int &c_dagger_2() {return c2dagger_;}
  inline unsigned int &c_1() {return c1_;}
  inline unsigned int &c_2() {return c2_;}
private:
  unsigned int z1_, z2_;
  unsigned int s1_, s2_;
  unsigned int c1_, c2_;
  unsigned int c1dagger_, c2dagger_;
  double abs_w_;
} vertex;



class inverse_m_matrix
{
public:
  alps::numeric::matrix<double> &matrix() { return matrix_;}
  alps::numeric::matrix<double> const &matrix() const { return matrix_;}
  std::vector<creator> &creators(){ return creators_;}
  const std::vector<creator> &creators() const{ return creators_;}
  std::vector<annihilator> &annihilators(){ return annihilators_;}
  const std::vector<annihilator> &annihilators()const{ return annihilators_;}
  std::vector<double> &alpha(){ return alpha_;}
  const std::vector<double> &alpha() const{ return alpha_;}
private:
  alps::numeric::matrix<double> matrix_;
  std::vector<creator> creators_;         //an array of creation operators c_dagger corresponding to the row of the matrix
  std::vector<annihilator> annihilators_; //an array of to annihilation operators c corresponding to the column of the matrix
  std::vector<double> alpha_;             //an array of doubles corresponding to the alphas of Rubtsov for the c, cdaggers at the same index.
};

std::ostream & operator<<(std::ostream &os, inverse_m_matrix const& M);



class InteractionExpansionSim: public alps::scheduler::MCSimulation, public alps::MatsubaraImpurityTask
{
public:

  InteractionExpansionSim(const alps::ProcessList &w, const boost::filesystem::path &p) : alps::scheduler::MCSimulation(w,p) {}
  
  InteractionExpansionSim(const alps::ProcessList &w, const alps::Parameters &p) : alps::scheduler::MCSimulation(w,p) {p_=p;}
  
  std::pair<matsubara_green_function_t,itime_green_function_t> get_result(); 

  void evaluate_selfenergy_measurement_matsubara(const alps::ObservableSet &gathered_measurements, 
                                                 matsubara_green_function_t &green_matsubara_measured,
                                                 const matsubara_green_function_t &bare_green_matsubara, 
                                                 std::vector<double>& densities, const double &beta, 
                                                 const int n_site, const int n_flavors, const int n_matsubara) const;

  void evaluate_selfenergy_measurement_itime_rs(const alps::ObservableSet &gathered_measurements, itime_green_function_t &green_result,
                                                const itime_green_function_t &green0, const double &beta, const int n_site, 
                                                const int n_flavors, const int n_tau, const int n_self) const;

  double green0_spline(const itime_green_function_t &green0, const itime_t delta_t, const int s1, const int s2, 
                       const spin_t flavor, int n_tau, double beta) const;
  
private:

  alps::Parameters p_;
};



class InteractionExpansionRun: public alps::scheduler::MCRun
{
public:

  InteractionExpansionRun(const alps::ProcessList &, const alps::Parameters &parms, int);
  ~InteractionExpansionRun() {}
  bool is_thermalized() const {return true;} //thermalization is done in the constructor. It's not a big deal here.
  void dostep();
  double work_done() const;
  
  // not yet implemented, see rubtosv_auxiliary.cpp
  void save(alps::ODump &) const;
  void load(alps::IDump &);
  
protected:
  
  /*functions*/
  /*io & initialization*/
  void initialize_simulation(const alps::Parameters&); // called by constructor
  // in file io.cpp
  void read_bare_green(std::ifstream &G0_omega, std::ifstream &G0_tau);
  void print(std::ostream &os);
  
  /*green's function*/
  // in file spines.cpp
  double green0_spline(const c_or_cdagger &cdagger, const c_or_cdagger &c) const;
  double green0_spline(const itime_t delta_t, const spin_t flavor, const site_t site1, const site_t site2) const;
  double green0_spline(const itime_t delta_t, const spin_t flavor) const;
  
  /*the actual solver functions*/
  // in file solver.cpp
  void interaction_expansion_step(void);
  void reset_perturbation_series(void);
  
  // in file fastupdate.cpp:
  double fastupdate_up(const int operator_nr, bool compute_only_weight);
  double fastupdate_down(const int operator_nr, const int flavor, bool compute_only_weight);
  
  /*measurement functions*/
  // in file measurements.cpp
  void measure_observables(void);
  void initialize_observables(void);
  
  void compute_W_matsubara();
  void compute_W_itime();
  void measure_Wk(std::vector<std::vector<std::valarray<std::complex<double> > > >& Wk, const unsigned int nfreq);
  void measure_densities();
  
  /*abstract virtual functions. Implement these for specific models.*/
  virtual double try_add()=0;
  virtual void perform_add()=0;
  virtual void reject_add()=0;
  virtual double try_remove(unsigned int vertex_nr)=0;
  virtual void perform_remove(unsigned int vertex_nr)=0;
  virtual void reject_remove()=0;
  
  /*private member variables, constant throughout the simulation*/
  const unsigned int max_order;                        
  const spin_t n_flavors;                                //number of flavors (called 'flavors') in InteractionExpansion
  const site_t n_site;                                //number of sites
  const frequency_t n_matsubara;        //number of matsubara freq
  const frequency_t n_matsubara_measurements;        //number of measured matsubara freq
  const itime_index_t n_tau;                        //number of imag time slices
  const itime_t n_tau_inv;                        //the inverse of n_tau
  const frequency_t n_self;                        //number of self energy (W) binning points
  const boost::uint64_t mc_steps;                        
  const unsigned long therm_steps;                
  const unsigned int nruns;                      
  const double max_time_in_seconds;
  
  const double beta;                                
  const double temperature;                        //only for performance reasons: avoid 1/beta computations where possible        
  const double onsite_U;                        
  const double alpha;                                
  const U_matrix U;
  
  
  const unsigned int recalc_period;                
  const unsigned int measurement_period;        
  const unsigned int convergence_check_period;        
  
  /*InteractionExpansion's roundoff threshold*/
  const double almost_zero;                        
  /*PRNG seed*/
  const int seed;                                
  
  /*private member variables*/
  matsubara_green_function_t green_matsubara;
  matsubara_green_function_t bare_green_matsubara;
  itime_green_function_t bare_green_itime;
  itime_green_function_t green_itime;
  std::vector<green_matrix> g0;
  boost::shared_ptr<FourierTransformer> fourier_ptr;
  
  vertex_array vertices;
  std::vector<inverse_m_matrix> M;
    
  double weight;
  double sign;
  unsigned int measurement_method;
  bool thermalized;
  
  simple_hist pert_hist;
  unsigned int hist_max_index;
  simple_hist **vertex_histograms;
  unsigned int vertex_histogram_size;
  
  unsigned long step;        
  time_t start_time;
  clock_t update_time;
  clock_t measurement_time;

};



/*aux functions*/
std::ostream& operator << (std::ostream& os, const std::vector<double>& v);
std::ostream& operator << (std::ostream &os, const vertex_array &vertices);
std::ostream& operator << (std::ostream &os, const vertex &v);
std::ostream& operator << (std::ostream &os, const c_or_cdagger &c);
std::ostream& operator << (std::ostream& os, const simple_hist &h);



//Use this for the most simple single site Hubbard model.
class HalfFillingHubbardInteractionExpansionRun: public InteractionExpansionRun{
public:
  HalfFillingHubbardInteractionExpansionRun(const std::vector<alps::Process, std::allocator<alps::Process> >& p, 
                                    const alps::Parameters& param, int& i)
    :InteractionExpansionRun(p, param, i)
  {
    if(n_flavors !=1){
      std::cerr<<"you need a different model for n_flavors!=1."<<std::endl;
      exit(1);
    }
  }
  double try_add();
  void perform_add();
  void reject_add();
  double try_remove(unsigned int vertex_nr);
  void perform_remove(unsigned int vertex_nr);
  void reject_remove();
};



class HubbardInteractionExpansionRun: public InteractionExpansionRun{
public:
  HubbardInteractionExpansionRun(const std::vector<alps::Process, std::allocator<alps::Process> >& p, const alps::Parameters& param, int& i)
    :InteractionExpansionRun(p, param, i){
    if(n_flavors !=2){
      std::cerr<<"you need a different model for n_flavors!=2."<<std::endl;
      exit(1);
    }
  }
  double try_add();
  void perform_add();
  void reject_add();
  double try_remove(unsigned int vertex_nr);
  void perform_remove(unsigned int vertex_nr);
  void reject_remove();
};



//Use this for multiple bands where you have terms Un_i n_j
class MultiBandDensityHubbardInteractionExpansionRun: public InteractionExpansionRun{
public:
  MultiBandDensityHubbardInteractionExpansionRun(const std::vector<alps::Process, std::allocator<alps::Process> >& p, const alps::Parameters& param, int& i)
    :InteractionExpansionRun(p, param, i) {
    if(n_site != 1){
      std::cerr<<"Please choose n_site=1."<<std::endl;
      exit(1);
    }
  }
  
  double try_add();
  void perform_add();
  void reject_add();
  double try_remove(unsigned int vertex_nr);
  void perform_remove(unsigned int vertex_nr);
  void reject_remove();
};

#endif
