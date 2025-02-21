/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2005 by Matthias Troyer <troyer@comp-phys.org>,
*                       Andreas Streich <astreich@student.ethz.ch>
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

/* $Id$ */

#ifndef ALPS_APPLICATIONS_MC_SPIN_FITTER_H_
#define ALPS_APPLICATIONS_MC_SPIN_FITTER_H_

#include <alps/scheduler.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cstdio>

#include "factory.h"
#include "helper.h"
#include <alps/plot.h>

#define SQUARE(x) ((x)*(x))
#define ABS(X) (X<0 ? -X : X)

template<class T>
bool is_contained(const std::vector<T>& vec, T elem)
{
  for (int i=0; i<vec.size(); i++)
    if (vec[i] == elem) return true;
  return false;
}

template <class T> inline T square(T x) { return x*x;}

typedef struct {
  alps::ParameterList parms;
  char type;
  bool is_genuine;
} FitterParamList;

alps::oxstream& operator<< (alps::oxstream& oxs, const FitterParamList list);
std::ifstream& operator>> (std::ifstream& is, const FitterParamList& list);

typedef struct {
  double T, mean, error;
} ResidualType;

alps::oxstream& operator<< (alps::oxstream& oxs, const ResidualType residual);
std::ifstream& operator>> (std::ifstream& is, ResidualType& residual);

typedef struct {
  std::vector<ResidualType> res_vals;
  double res_size;
  double avg_res_size;
} ResHistElementType;
  
alps::oxstream& operator<< (alps::oxstream& oxs, const ResHistElementType res); 
std::ifstream& operator>> (std::ifstream& is, ResHistElementType& reshist);

alps::oxstream& operator<< (alps::oxstream& oxs, 
               const alps::scheduler::ResultType);
std::ifstream& operator>> (std::ifstream& is, alps::scheduler::ResultType& res);

alps::oxstream& operator<< (alps::oxstream& oxs, 
               const alps::scheduler::ResultsType);
std::ifstream& operator>> (std::ifstream& is,
               alps::scheduler::ResultsType& results);

typedef struct {
  double min, max;
  bool validmin, validmax;
} Limit;

alps::oxstream& operator<< (alps::oxstream& oxs, const Limit);
std::ifstream& operator>> (std::ifstream& is, Limit&);

typedef struct {
  std::string name;
  Limit limit;
} FittingParam;

alps::oxstream& operator<< (alps::oxstream& oxs, const FittingParam);
std::ifstream& operator>> (std::ifstream& is, FittingParam&);

class AbstractFitter 
{
public:
  AbstractFitter(const std::string& sim_file_name);
  virtual void initialize();
  virtual void initialize(std::ifstream& is);

  virtual ~AbstractFitter() { };

  virtual void add_sim_result(const alps::scheduler::ResultsType sim_res);
  virtual void make_new_jobfile();
  std::string get_current_filename();
  virtual bool is_done(double* time_left);

  virtual void print_results() const;
  virtual void print_results(const std::string filename) const;
  virtual alps::oxstream& write_to_xml(alps::oxstream& oxs) const;
  virtual std::string get_fitter_name() const =0;

  alps::plot::Set<double> get_exp_results_plot_set() const;
  alps::plot::Set<double> get_sim_results_plot_set() const;
  alps::plot::Set<double> get_exp_error_plot_set() const;
  alps::plot::Set<double> get_res_abs_plot_set() const;
  alps::plot::Set<double> get_exp_percent_plot_set() const;
  virtual alps::plot::Set<double> get_res_percent_plot_set() const;

  int get_step_count() const { return step; }
  std::string get_filenamebase() const { return filenamebase; }
  std::string get_fit_observable() const { return fit_obs_name; }

protected:
  void read_target(const std::string file_name);
  std::string print_to_file(const alps::ParameterList params) const;

  std::vector<ResidualType> get_residuum() const;
  std::vector<ResidualType> get_residuum(unsigned int pos) const;
  std::vector<ResidualType> get_residuum(unsigned int pos,
               const alps::scheduler::ResultsType exp_res) const;
  double get_residuum_size(const std::vector<ResidualType>) const;
  double get_avg_res_size(const std::vector<ResidualType>) const;
  double get_last_res_size() const;
  alps::scheduler::ResultsType get_eff_sim_res(alps::scheduler::ResultsType,
               const alps::ParameterList) const;

protected:
  alps::ParameterList adjust_parameters(const alps::ParameterList list, 
               const alps::scheduler::ResultsType& exp_res);

  virtual void find_new_parameters()=0; 
  std::vector<FittingParam> find_variable_parms(alps::Parameters parms) const;
  alps::ParameterList update_parameters(const std::string par_name, 
               double new_value, alps::ParameterList list) const;

  alps::scheduler::ResultsType all_exp_res;
  std::vector<alps::scheduler::ResultsType> sim_results;
  std::vector<FitterParamList> param_history;
  FitterParamList current_params;
  std::list<FitterParamList> next_parameters;
  std::vector<ResHistElementType> residuum_history;
  int last_genuine;
  std::vector<FittingParam> variable_params;
  int num_variable_params;
  double tolerated_error_factor;

  std::string fit_obs_name;  
  std::string filenamebase;
  std::string results_namebase;
  std::fstream tracefile;
 
  std::string current_filename;
  boost::posix_time::ptime start_time, end_time;
  int step;
  int max_steps;
  double max_avg_residuum;
 
};

class DummyFitter : public AbstractFitter 
{
public:
  DummyFitter(const std::string& filename);
  virtual std::string get_fitter_name() const;
  virtual bool is_done(double* seconds);

protected:
  virtual void find_new_parameters();
};

class BaseIncrFitter : public AbstractFitter
{
public:
  BaseIncrFitter(const std::string& filename);
  virtual alps::plot::Set<double> get_res_percent_plot_set() const;

protected:
  void initialize(std::ifstream& is);
  alps::oxstream& write_to_xml(alps::oxstream& oxs) const;

  void determine_curvatures();
  void choose_exp_results(int expected, int min, int max);

  bool enlarged_basis;
  bool should_enlarge_basis();
  alps::ParameterList enlarge_basis(const alps::ParameterList list);

  void add_exp_results(int expected, int min, int max);
  bool is_contained(const alps::scheduler::ResultType& element) const;

  alps::scheduler::ResultsType exp_results;
  std::vector<double> curvature;
  double curvature_sum;
  double remaining_curv_sum;

  void choose_exp_results(int expected, double rem_curv_sum);
  void choose_exp_results(int expected, int min);
  void save_choose(int expected, int min, int max, int nbiggest);
  void simple_add_exp_results(int exp, double rem_sum,
               std::vector<unsigned int>& tmp) const;
  void simple_add_exp_results(int exp, double rem_sum,
               const std::vector<unsigned int>& positions,
               std::vector<unsigned int>& tmp) const;
  void save_add_exp_results(int exp, int nbiggest,
               std::vector<unsigned int>& tmp) const;

  double enlarge_trigger_factor;
  int steps_between_baseincr;
  unsigned int last_base_increment;
};
  
class EstGradFitter : public BaseIncrFitter
{
 
public:
  EstGradFitter(const std::string filename);
  virtual void initialize();
  virtual void initialize(std::ifstream& is);

  virtual ~EstGradFitter();
  virtual std::string get_fitter_name() const;

  void make_new_jobfile();
  virtual void add_sim_result(const alps::scheduler::ResultsType sim_res);
  virtual bool is_done(double* time_left);
  virtual alps::oxstream& write_to_xml(alps::oxstream& oxs) const;

protected:
  alps::ParameterList init_param_list(
               const alps::ParameterList list, double fac);
  virtual void find_new_parameters();
  void prepare_jacobi_approx(const alps::ParameterList list);
  void compute_gauss_newton_step();
  void do_step(alps::ParameterList& plist);
  void find_variable_model_params( const alps::Parameters& sim_param);

  bool use_curie_term;
  std::vector<unsigned int> varPar_history_map;

  // do we use the Levenberg Marquart algorithm?
  bool use_levenberg_marquart;
  virtual double determine_lambda(double size);

  // do we accept new genuine parameter sets only if they yield smaller residui?
  bool enforce_amelioration;
  unsigned int decrease_step_counter;
  unsigned int decrease_step_limit;
  void escape_from_decreasing();
  double stepsize_prefactor;

  std::vector<double> last_step;
  double res_to_underbid;
  unsigned int res_to_underbid_index;
  int best_residui_size[3];

  bool decrease_stepsize();
  double step_decrease_factor;
  // do we check whether the error in the simulation is sufficiently small?
  bool check_error_level;
  double error_decrement;
  int erroneous_accept_count;
  double erroneous_accept_ratio;

  int check_errors();
  alps::ParameterList thighten_error_limit(const alps::ParameterList list);

  bool have_tightened;
  std::vector<unsigned int> thightened_indices;
};

#endif
