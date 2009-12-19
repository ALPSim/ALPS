/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2005 by Matthias Troyer <troyer@comp-phys.org>,
*                       Walter Gander <gander@inf.ethz.ch>,
*                       Andreas Streich <astreich@student.ethz.ch>
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

/* $Id$ */

#include "fitter.h"
#include "lapack.h"

#define STEP_DECREASE_FACTOR 0.316227766 // = 1./sqrt(10)
#define STEP_INCREASE_FACTOR 5.0

/**
 * The constructor of the class. Read in parameters or set default values.
 *
 * @param filename The name of the file constaining the fitting options.
 */
EstGradFitter::EstGradFitter(const std::string filename)
   : BaseIncrFitter(filename)
{
  alps::Parameters sim_param;
  {
    std::ifstream in(filename.c_str());
    in >> sim_param;
  }

  tracefile << "Constructor of EstGradFitter:\n"
            << "=============================\n"
            << "read simulation parameters:\n" << sim_param << "\n";

  // read the parameters determining the fitting algorithm
  if (sim_param.defined("CHECK_ERROR_LEVEL") 
               && sim_param["CHECK_ERROR_LEVEL"]=="no") {
    check_error_level = false;
    tracefile << "Error level will not be checked\n";
  } else {
    check_error_level = true;
    tracefile << "Error Level will be checked\n";
  } 
  tolerated_error_factor = sim_param.value_or_default("TOLERATED_ERROR_FACTOR",
               0.3333333333333333);
  erroneous_accept_count = (int)(sim_param.value_or_default(
               "ERRONEOUS_ACCEPT_COUNT", 1));
  erroneous_accept_ratio = (double)(sim_param.value_or_default(
               "ERRONEOUS_ACCEPT_RATIO", 0.1));
  error_decrement = sim_param.value_or_default("ERROR_DECREMENT",
               0.25*tolerated_error_factor);

  if (sim_param.defined("USE_LEVENBERG_MARQUART") 
               && sim_param["USE_LEVENBERG_MARQUART"]=="yes") {
    use_levenberg_marquart = true;
    tracefile << "Levenberg-Marquart algorithm will be used\n";
  } else {
    use_levenberg_marquart = false;
    tracefile << "Levenberg-Marquart algorithm will not be used\n";
  }

  if (sim_param.defined("ENFORCE_AMELIORATION") 
               && sim_param["ENFORCE_AMELIORATION"]=="yes") {
    enforce_amelioration = true;
    tracefile << "Amelioration will be enforced\n";

  } else {
    enforce_amelioration = false;
    tracefile << "Amelioration will not be enforced\n";
  }
  step_decrease_factor = sim_param.value_or_default("STEP_DECREASE_FACTOR",
               0.316227766); // = 1/sqrt(10)
  decrease_step_limit = sim_param.value_or_default("DECREASE_STEP_LIMIT", 4);

  alps::ParameterList init_params;
  {
    std::ifstream param_stream(filenamebase.c_str());
    param_stream >> init_params;
  }
  FitterParamList initial;
  initial.parms = AbstractFitter::adjust_parameters(init_params, exp_results);
  initial.is_genuine = true;
  initial.type = 'g';
  FittingParam myfp;
  myfp.limit.validmin = false;
  myfp.limit.validmax = false;

  if (sim_param.defined("FIT_CONSTANT")) {
    initial.parms = update_parameters("FIT_CONSTANT",
               sim_param["FIT_CONSTANT"], initial.parms);
  }
  if (sim_param.defined("CONST_PREFAC")) {
    initial.parms = update_parameters("CONST_PREFAC",
               sim_param["CONST_PREFAC"], initial.parms);
  }
  if (sim_param.defined("POS_PREFAC")) {
    initial.parms = update_parameters("POS_PREFAC",
               sim_param["POS_PREFAC"], initial.parms);
  }

  if (sim_param.defined("C_OVER_T")) {
    initial.parms = update_parameters("C_OVER_T",
               sim_param["C_OVER_T"], initial.parms);
  }

  if (sim_param.defined("USE_CW_LAW") &&
                (sim_param["USE_CW_LAW"]=="yes")) {
    use_curie_term = true;

    initial.parms = update_parameters("CURIE_CONST",
               sim_param["CURIE_CONST"], initial.parms);
    initial.parms = update_parameters("CURIE_TEMP",
               sim_param["CURIE_TEMP"], initial.parms);
  }
  // possibly more to come

  else
    use_curie_term = false;

  // determine which model parameters are variable
  find_variable_model_params(sim_param);

  next_parameters.push_back(initial);
  prepare_jacobi_approx(initial.parms);
}

/**
 * Determine the model parameters which are variable and must be fited.
 *
 * @param sim_param The parameters of the fitting.
 */
void EstGradFitter::find_variable_model_params(
               const alps::Parameters& sim_param)
{
  if (!sim_param.defined("VAR_MODEL_PARAMS"))
    return;
  std::string tmp_str = sim_param["VAR_MODEL_PARAMS"];
  Helper<std::string> mh(tmp_str);
  int count = mh.elemcount();
  std::vector<std::string> vmp;
  vmp.resize(count);
  for (int i=0; i<count; i++)
    vmp[i] = mh.getNextString();

  FittingParam fp;
  fp.limit.validmin = false;
  fp.limit.validmax = false;

  // check all possible model parameters whether they are variable
  if (::is_contained(vmp, std::string("FIT_CONSTANT"))) {
    fp.name = "FIT_CONSTANT";
    variable_params.push_back(fp);
  }
  if (::is_contained(vmp, std::string("CONST_PREFAC"))) {
    fp.name = "CONST_PREFAC";
    variable_params.push_back(fp);
  }
  if (::is_contained(vmp, std::string("POS_PREFAC"))) {
    fp.name = "POS_PREFAC";
    variable_params.push_back(fp);
  }
  if (::is_contained(vmp, std::string("C_OVER_T"))) {
    fp.name = "C_OVER_T";
    variable_params.push_back(fp);
  }
  if ((::is_contained(vmp, std::string("CURIE_CONST"))) && use_curie_term) {
    fp.name = "CURIE_CONST";
    variable_params.push_back(fp);
  }
  if ((::is_contained(vmp, std::string("CURIE_TEMP"))) && use_curie_term) {
    fp.name = "CURIE_TEMP";
    variable_params.push_back(fp);
  }
}

/**
 * Initialize the class reading information from a in file stream in XML format.
 *
 * @param is The in file stream.
 */
void EstGradFitter::initialize(std::ifstream& is)
{
  BaseIncrFitter::initialize(is);

  // read parameters from the checkpoint file
  alps::XMLTag tag = alps::parse_tag(is, true);
  enlarge_trigger_factor = alps::evaluate<double>(alps::parse_content(is));
  std::cout << "enlarge_trigger_factor set to " << enlarge_trigger_factor 
            << "\n";
  tag = alps::parse_tag(is, true);
  tag = alps::parse_tag(is, true);

  last_base_increment = 
               (unsigned int)(alps::evaluate<float>(alps::parse_content(is)));
  std::cout << "Last base increment was done at step " << last_base_increment 
            << "\n";
  tag = alps::parse_tag(is, true);
  tag = alps::parse_tag(is, true);

  enlarged_basis = (alps::parse_content(is) == "1");
  tag = alps::parse_tag(is, true);
  tag = alps::parse_tag(is, true);

  decrease_step_counter=(int)(alps::evaluate<float>(alps::parse_content(is)));
  std::cout << "Stepsize was decreased " << decrease_step_counter 
            << " times in a row.\n";
  tag = alps::parse_tag(is, true);
  tag = alps::parse_tag(is, true);

  stepsize_prefactor = (alps::evaluate<float>(alps::parse_content(is)));
  std::cout << "Stepsize decrease factor is " << stepsize_prefactor << "\n";
  tag = alps::parse_tag(is, true);
  tag = alps::parse_tag(is, true);

  res_to_underbid = alps::evaluate<double>(alps::parse_content(is));
  std::cout << "Residuum to underbid is " << res_to_underbid << "\n";
  tag = alps::parse_tag(is, true);
  tag = alps::parse_tag(is, true);

  res_to_underbid_index=(int)(alps::evaluate<float>(alps::parse_content(is)));
  tag = alps::parse_tag(is, true);
  std::cout << "Its index is " << res_to_underbid_index << "\n";
  tag = alps::parse_tag(is, true);

  have_tightened = (alps::parse_content(is) == "1");
  std::cout << "have_tightened set to " << have_tightened << "\n";
  tag = alps::parse_tag(is, true);
  tag = alps::parse_tag(is, true);

  num_variable_params = (int)(alps::evaluate<float>(alps::parse_content(is)));
  std::cout << "variable parameters:\n";
  tag = alps::parse_tag(is, true);

  tag = alps::parse_tag(is, true);

  FittingParam fp;
  if (tag.type == alps::XMLTag::SINGLE) {
    std::cout << "Invalid checkpoint file: No variable parameters.\n";
    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  } else {
    variable_params.clear();
    tag = alps::parse_tag(is,true);
    while (tag.name != "/FITTING_PARAMETER") {
      is >> fp;
      variable_params.push_back(fp);
      tag = alps::parse_tag(is, true);
      tag = alps::parse_tag(is, true);
      tag = alps::parse_tag(is, true);
    }
  }

  tag = alps::parse_tag(is, true);
  if (tag.type == alps::XMLTag::SINGLE)
    std::cout << "No Values for last step given.\n";
  else {
    last_step.clear();
    std::cout << "Values for the last step:\n";
    tag = alps::parse_tag(is,true);
    while (tag.name != "/FITTING_PARAMETER") {
      double component = alps::evaluate<double>(alps::parse_content(is));
      std::cout << component  << "\n";
      last_step.push_back(component);
      tag = alps::parse_tag(is, true);
      tag = alps::parse_tag(is, true);
    }
  }

  tag = alps::parse_tag(is, true);
  if (tag.name == "FITTING_PARAMETER") {
    for (int i =0; i<3; i++) {
      tag = alps::parse_tag(is,true);
      int current = (int)(alps::evaluate<double>(alps::parse_content(is)));
      best_residui_size[i] = current;
      tag = alps::parse_tag(is, true);
    }
    tag = alps::parse_tag(is, true);
    tag = alps::parse_tag(is, true);
  }

  thightened_indices.clear();
  if (tag.type == alps::XMLTag::SINGLE)
    std::cout << "No thightening indices given in the checkpoint file.\n";
  else {
    tag = alps::parse_tag(is,true);
    while (tag.name != "/FITTING_PARAMETER") {   
      int pos = (int)(alps::evaluate<float>(alps::parse_content(is)));
      thightened_indices.push_back(pos);
      tag = alps::parse_tag(is, true);
      tag = alps::parse_tag(is, true);
    }
  }
  std::cout << "Initialization from checkpoint file successfully ended.\n";
}

/**
 * Initialize the class with default values.
 */
void EstGradFitter::initialize()
{
  step = -1;
  last_base_increment = 0;
  have_tightened = false;
  enlarged_basis = false;
  stepsize_prefactor = 1.0;

  tracefile << "\nThe variable parameters are:\n";
  int i;
  for (i=0; i<num_variable_params; i++)
    tracefile << " " << variable_params[i].name
              << " (derivation numerically approximated)\n";
  for (; i<variable_params.size(); i++)
    tracefile << " " << variable_params[i].name
              << " (analytic derivation)\n";

  // set the ranking list elements to -1, i.e. not set yet
  best_residui_size[0] = -1;
  best_residui_size[1] = -1;
  best_residui_size[2] = -1;
}

EstGradFitter::~EstGradFitter()
{
  tracefile.flush();
  tracefile.close();
}

/**
 * Initialize the parameter list for a first use, ie. by setting the summary
 * variable and the the error limit.
 * 
 * @param list The list of parameters to be initialized
 * @param fac  The error limit is set to fac*mean_value
 */
alps::ParameterList EstGradFitter::init_param_list(
               const alps::ParameterList list, double fac)
{
  alps::ParameterList res = list;
  unsigned int size = res.size();
  // ensure the summary variable is the fitted observable
  alps::Parameter pSUM("SUMMARY_VARIABLE",fit_obs_name.c_str());
  for (int i=0; i<size; i++) {
    alps::Parameter pERR("ERROR_LIMIT",fac*exp_results[i].mean);
    res[i] << pERR;
    res[i] << pSUM;
  }
  return res;  
}

/**
 * Check the error level of the current residual.
 */
int EstGradFitter::check_errors()
{
  int count = 0;
  if (residuum_history.empty())
    return 0;
  unsigned int run_number = residuum_history.size()-1;
  unsigned int size = residuum_history[run_number].res_vals.size();

  tracefile << "\nChecking error limits of the residui:\n";
  for (int i=0; i<size; i++) {
    double mean = residuum_history[run_number].res_vals[i].mean;
    mean = std::abs(mean);
    double error = residuum_history[run_number].res_vals[i].error;
    error = std::abs(error);

    tracefile << " i = " << i << " mean = " << mean << " error = " << error;
    if (error > tolerated_error_factor*mean) {
      count++;
      tracefile << " : not accepted\n";
    } else {
      tracefile << " : accepted\n";
    }
  }
  tracefile << count << " of " << size << " not accepted\n";
  return count;
}

/**
 * Add a simulation results.
 *
 * @param sim_res The simulation result.
 */
void EstGradFitter::add_sim_result(const alps::scheduler::ResultsType sim_res)
{
  tracefile.flush();
  if (!have_tightened) {
    sim_results.push_back(sim_res);
    ResHistElementType cur_res;
    cur_res.res_vals = get_residuum(step, exp_results);
    cur_res.res_size = get_residuum_size(cur_res.res_vals);
    cur_res.avg_res_size = cur_res.res_size/cur_res.res_vals.size();
    residuum_history.push_back(cur_res);
                                                                                
    // simulation results and parameters only make sence together
    // the parameters are also put into history here
    param_history.push_back(current_params);
  } else {
    // copy the previous results and replace the simulations with thightened
    // error limit by the more exact ones of this run
    alps::scheduler::ResultsType temp = sim_results[sim_results.size()-1];
    for (int i =0; i<thightened_indices.size(); i++)
      temp[thightened_indices[i]] += sim_res[i];
    sim_results.push_back(temp);
    
    // update residuum history
    ResHistElementType cur_res;
    cur_res.res_vals = get_residuum(step, exp_results);
    cur_res.res_size = get_residuum_size(cur_res.res_vals);
    cur_res.avg_res_size = cur_res.res_size/cur_res.res_vals.size();
    residuum_history.push_back(cur_res);

    // simulation results and parameters only make sence together
    // the parameters are also put into history here

    // re-assemble the complete parameter set
    FitterParamList p_temp = param_history[param_history.size()-1];
    for (int i=0; i<thightened_indices.size(); i++)
      p_temp.parms[thightened_indices[i]] = current_params.parms[i];
    p_temp.is_genuine = current_params.is_genuine;
    p_temp.type = 't';
    param_history.push_back(p_temp);
    current_params.parms = p_temp.parms;
  }

  tracefile << "\n\nStep " << std::setw(4) << step << "\n=========\n"
            << "Added a simulation result:\n";
  alps::scheduler::ResultsType temp = sim_results[sim_results.size()-1];  
  for (int i=0; i<temp.size(); i++)
    tracefile << " " << temp[i].T << " " << temp[i].mean << " "
              << temp[i].error << " " << temp[i].count <<"\n";
  tracefile << "there are now " << sim_results.size() << " sim_results.\n";
  tracefile.flush();

  // update the residui ranking list
  if ((!current_params.is_genuine) && (!have_tightened)) {
    double this_res = residuum_history[step].avg_res_size; 
    if (best_residui_size[0] == -1)
      best_residui_size[0] = step;
    else if (residuum_history[best_residui_size[0]].avg_res_size > this_res) {
      best_residui_size[2] = best_residui_size[1];
      best_residui_size[1] = best_residui_size[0];
      best_residui_size[0] = step;
    } else if (best_residui_size[1] == -1) {
      best_residui_size[2] = best_residui_size[1];
      best_residui_size[1] = step;
    } else if (residuum_history[best_residui_size[1]].avg_res_size > this_res) {
      best_residui_size[2] = best_residui_size[1];
      best_residui_size[1] = step;
    } else if (best_residui_size[2] == -1)
      best_residui_size[2] = step;
    else if (residuum_history[best_residui_size[2]].avg_res_size > this_res)
      best_residui_size[2] = step;
  }
}

/**
 * Set a tighter error limit to the parameter list
 * 
 * @param list The list which error limits have to be tightened.
 */
alps::ParameterList EstGradFitter::thighten_error_limit(
               const alps::ParameterList list)
{
  alps::ParameterList new_list;
  unsigned int size = list.size();
  unsigned int run_number = residuum_history.size()-1;
  tracefile << "thightening error limit\n";
  ResHistElementType res = residuum_history[run_number];
  thightened_indices.clear(); 
  
  tracefile << " The following simulations will be re-run with a thighter "
            << "error limit:\n ";
  for (unsigned int i=0; i<size; i++) {
    double mean = residuum_history[run_number].res_vals[i].mean;
    mean = std::abs(mean);
    double error = residuum_history[run_number].res_vals[i].error;
    error = ABS(error);
    if (error > tolerated_error_factor*mean) {
      // not accepted, needs to be thightened
      alps::Parameter p("ERROR_LIMIT",error_decrement*mean);
      alps::Parameters tmp = list[i];
      tmp << p;
      new_list.push_back(tmp); 
      thightened_indices.push_back(i);
      tracefile << i << " ";
    }
  }
  tracefile << "\n";
  return new_list;
}

/**
 * Checks whether the fitting is done.
 *
 * @param seconds_left If the fitting is not done yet, the number of seconds
 *                     left is indicated.
 */
bool EstGradFitter::is_done(double* seconds_left)
{
  tracefile << "\nChecking whether it's done:\n";
  if (end_time > start_time) {
    boost::posix_time::time_duration time_left = end_time 
                 - boost::posix_time::second_clock::local_time();
    *seconds_left = time_left.total_seconds();
    if (*seconds_left <= 0) {
      tracefile << " Time limit exceeded, ending\n";
      std::cout << " Time limit exceeded, ending\n";
      return true;
    } else {
      tracefile << " Still " << *seconds_left << " seconds left.\n";
    }
  } else
    *seconds_left = -1;

  tracefile << " Did " << step << " steps, maximal number of steps is "
            << max_steps << "\n";

  if (step >= max_steps) {
    std::cout << "maximal number of steps has been reached, ending";
    tracefile << " maximal number of steps has been reached, ending.\n";
    return true;
  }
  if (step < 0)
    return false;

  double last_res;

  if (residuum_history.empty()) {
    tracefile << "No residuum computed so far.\n";
    last_res = max_avg_residuum + 1.;
    return false;
  } else {
    last_res = residuum_history.back().avg_res_size;
    tracefile << " Last residuum size is " << last_res << ", limit is "
              << max_avg_residuum << "\n";
  }

  if (last_res < max_avg_residuum) {
    std::cout << "\nGot a sufficiently small residuum.\n";
    tracefile << " Got a sufficiently small residuum.\n";

    if (all_exp_res.size() == exp_results.size()) {
      tracefile << " All experimental data taken into account, ending\n";
      return true;
    } else {
      std::cout << " Not all experimental data taken into account.\n"
                << " Only " << exp_results.size() << " of " <<all_exp_res.size()
                << " experimental data points are currently used for fitting."
                << "\nFitting basis will be enlarged at the next occasion.\n";
      tracefile << " Not all experimental data taken into account.\n"
                << " Only " << exp_results.size() << " of " <<all_exp_res.size()
                << " experimental data points are currently used for fitting."
                << "\n Fitting basis will be enlarged at the next occasion.\n";

      return false;
    }
  }
  return false;
}

/**
 * Escape from a series of unsuccesful step size decrements, either by enlarging
 * the fitting base or by using a non-genuine parameter set as new starting
 * point for an approximation of the Jacobian matrix.
 */
void EstGradFitter::escape_from_decreasing()
{
  tracefile<< "\nEscaping from a series of ineffective step size decrements:\n";
  tracefile.flush();

  // if possible: Increase the basis
  unsigned int size_now = exp_results.size();
  add_exp_results(size_now, size_now-2, size_now+2);

  // take the best available non-genuine parameter set and use it for an
  // approximation of the Jacobian matrix
  if (best_residui_size[0] == -1) {
    // no better non-genuine parameter set is available
    // accept the step anyway
    // this case should actually not happen if there is at least one
    // parameter with numerical differentiation
    last_step.clear();
    prepare_jacobi_approx(param_history[last_genuine].parms);
    tracefile << " no non-genuine parameter sets available, using the "
              << "parameters\n with the current step size anyway.\n";
  } else {
    tracefile << " using one of the parameter sets used for the numerical "
              <<" derivation as new starting point.\n The parameter set number "
              << best_residui_size[0] << " will be used.\n";
    last_genuine = best_residui_size[0];

    // there's a good non-genuine parameter set available
    // if the size of the corresponding residui is the same as now:
    // set the index as last_genuine, make jacobi approximation
    if (sim_results[best_residui_size[0]].size() == exp_results.size()) {
      prepare_jacobi_approx(param_history[best_residui_size[0]].parms);
      tracefile << " The fitting base has not changed since, use the then "
                << "simulation results,\n directly prepare jacobi "
                << "approximation.\n";
    } else {
      // the fitting basis has changed, ie. make a new simulation.
      FitterParamList myList;
      myList.parms =adjust_parameters(param_history[best_residui_size[0]].parms,
               exp_results);
      myList.is_genuine = true;
      myList.type = 'g';
      next_parameters.push_back(myList);
      prepare_jacobi_approx(myList.parms);
      tracefile << " The fitting base has changed since, recompute the "
                << "simulation with the then parameters.\n A Jacobi "
                << "approximation will be done afterwards\n";
    }
    // delete the used best non-genuine parameter set
    best_residui_size[0] = best_residui_size[1];
    best_residui_size[1] = best_residui_size[2];
    best_residui_size[2] = -1;
    varPar_history_map.clear();
    tracefile.flush();
  }
}

/**
 * Generate a new jobfile and the corresponding task files. New parameters are
 * computed if no precomputed parameter values are available.
 */
void EstGradFitter::make_new_jobfile()
{
  // make_new_jobfile is the beginning of a new fitting step
  step++;
  bool read_next_parameter_set;
  if ((current_params.parms.size() > 0) && (check_error_level)) {
    int deficient_count = check_errors();
    tracefile << "there are " << deficient_count 
              << " simulations with too high error level.\n";

    if ((deficient_count > erroneous_accept_count) 
               // take all simulations into account, also the ones that are 
               // not re-computed (e.g. after an error thightening)
               &&(deficient_count>erroneous_accept_ratio*exp_results.size())) {
      if (have_tightened) {
        std::cerr << "\nWARNING : Requested error limit not achieved, too low "
                  <<              "sweep limit"
                  << "\n          Continuing with data with high relative "
                                  "errors."
                  << "\n          Please increase the number of sweeps\n\n";
        tracefile << "\nWARNING : Requested error limit not achieved, too low "
                  <<              "sweep limit"
                  << "\n          Continuing with data with high relative "
                                  "errors."
                  << "\n          Please increase the number of sweeps\n";

        read_next_parameter_set = true;
        have_tightened = false;
      } else {
        tracefile << "\nToo high error level in residuum, thighening error "
                  << "limit\n";
        current_params.parms = thighten_error_limit(current_params.parms);
        have_tightened = true;
        read_next_parameter_set = false;
      }
    } else {
      read_next_parameter_set = true;
      have_tightened = false;
      tracefile<< "\nError level of the current simulation is ok, continuing\n";
    }
  } else
    // no error level check was done
    read_next_parameter_set = true;

  if (read_next_parameter_set) {
    // accept the results, housekeeping
    if (current_params.is_genuine) {
      last_genuine = param_history.size() - 1;
    } else if (param_history.size() > 0) 
      // update the lookup-table vector
      varPar_history_map.push_back(param_history.size() - 1);

    if (next_parameters.empty())
      find_new_parameters();

    current_params = next_parameters.front();
    next_parameters.pop_front();
  }
  current_filename = print_to_file(current_params.parms);
}

void EstGradFitter::find_new_parameters()
{
  FitterParamList next_params;
  double this_residuum;
  if (last_step.empty()) {
    // make a Gauss-Newton step
    compute_gauss_newton_step();
    if (should_enlarge_basis())
      next_params.parms = enlarge_basis(param_history[last_genuine].parms);
    else
      next_params.parms = param_history[last_genuine].parms;
    do_step(next_params.parms);
    next_params.is_genuine = true;
    next_params.type = 'g';
    next_parameters.push_back(next_params);
  } else if (enforce_amelioration) {
    this_residuum = residuum_history[step-1].avg_res_size;
    tracefile << "\nChecking whether the stepsize is accepted:\n"
              << " this_residuum = "<< this_residuum << ", res_to_underbid = "
              << res_to_underbid << "\n";
    if (this_residuum < res_to_underbid) {
      // the current parameter set is accepted, prepare another approximation
      // of the Jacobian matrix
      last_step.clear();
      tracefile << " step accepted, preparing Jacobi matrix approximation\n";
      decrease_step_counter = 0;
      prepare_jacobi_approx(param_history[last_genuine].parms);
    } else if (decrease_step_counter > decrease_step_limit) {
      // we got a too high residuum and have been decreasing many times.
      // try to escape from the decreasings
      escape_from_decreasing();
      last_step.clear();
      decrease_step_counter = 0;
    } else {
      // residuum is too high, continue decreasing the step size
      tracefile << " step not accepted, decrease_stepsize\n";
                                                                                
      decrease_stepsize();
      decrease_step_counter++;
      next_params.parms = param_history[res_to_underbid_index].parms;
      if (next_params.parms.size() != current_params.parms.size())
        next_params.parms = adjust_parameters(next_params.parms, exp_results);
      do_step(next_params.parms);
      next_params.is_genuine = true;
      next_params.type = 'd';
      next_parameters.push_back(next_params);
    }
  } else {
    // no enforcement, just continue with a jacobi approximation
    last_step.clear();
    prepare_jacobi_approx(param_history[last_genuine].parms);
  }
}

/**
 * Return the name of the fitter.
 */
std::string EstGradFitter::get_fitter_name() const
{
  std::string theName("Estimated gradient fitter");
  if (use_levenberg_marquart)
    theName += " with Levenberg Marquart stabilization";
  if (enforce_amelioration)
    theName += "\nSize of new step is decreased until a shorter residuum than the current one is found";
  if (check_error_level)
    theName += "\nError level of the residuum is checked";
  return theName;
}

/**
 * Prepare an approximation of the Jacobian matrix by adding a small delta to 
 * all variable parameters.
 * 
 * @param list The current parameters.
 */
void EstGradFitter::prepare_jacobi_approx(alps::ParameterList list)
{
  // for every variable parameter:
  // add a small 'delta' to it, in order to approximate the gradient
  tracefile << "\nStarting prepare_jacobi_approx()\n";
  if (num_variable_params == 0) {
    tracefile << "\nNO VARIABLE PARAMETERS FOUND !!!!\n\n";
    std::cerr << "\nno variable parameters found !!!!\n\n";
    boost::throw_exception(std::runtime_error("no variable parameters found"));
  }
  tracefile << " There are " << num_variable_params << " variable parameters\n";
  FitterParamList tmp;
  std::string par_name;
  double old_value, new_value;
  for (int i = 0; i<num_variable_params; i++) {
    par_name = variable_params[i].name;
    old_value = list[0][par_name];
    if (old_value != 0.0)
      new_value = 1.05*(double)(list[0][par_name]);
    else 
      new_value = 0.001;
    tmp.parms = update_parameters(par_name, new_value, list);
    tracefile << " Parameter " << par_name << " will be set to " << new_value
              << ".\n";
    tmp.is_genuine = false;
    tmp.type = 'a';
    next_parameters.push_back(tmp);
  }
}

/**
 * Determine the lambda factor for the Levenberg-Marquart stabilization.
 */
double EstGradFitter::determine_lambda(double size)
{
  if (enforce_amelioration) {
    tracefile << "\nComputing Lambda for the Levenberg-Marquart stabilization"
              << " (with enforcement):\n Lambda set to 0.0001.\n";
    return 0.0001;
  }

  // determine c
  double c;
  if (size > 10.0)
    c = 10.0;
  else if (size > 1.0)
    c = 1.0;
  else
    c = 0.01;

  double lambda = c*size;
  tracefile << "\nComputing Lambda for the Levenberg-Marquart stabilization:\n"
            << " Residuum size is " << size << ", c set to " << c
            << "\n Lambda set to " << lambda << "\n";
  return lambda;
}

/**
 * Decrease the size of the last Gauss-Newton step.
 */
bool EstGradFitter::decrease_stepsize()
{
  unsigned int size = last_step.size();
  for (int i=0; i<size; i++)
    last_step[i] *= step_decrease_factor;
//  stepsize_prefactor *= STEP_DECREASE_FACTOR;
  return true;
}

/**
 * Compute a Gauss-Newton step.
 */
void EstGradFitter::compute_gauss_newton_step()
{
  tracefile << "\nComputing a Gauss-Newton step:\n";
  
  unsigned int result_count = exp_results.size();

  // we make one row per variable parameter, no matter whether the derivation is
  // is approximated of computed analytically.
  unsigned int col_count = variable_params.size();

  unsigned int row_count = result_count;
  if (use_levenberg_marquart) 
    // add rows for the multiple of the identity matrix and the zeros
    row_count += col_count;

  res_to_underbid = residuum_history[last_genuine].avg_res_size;
  res_to_underbid_index = last_genuine;
  
  // prepare the Jacobian matrix (containing the gradients of the residuum)
  boost::numeric::ublas::matrix< double, boost::numeric::ublas::column_major, 
                boost::numeric::ublas::unbounded_array<double> > 
                jacobian(row_count,col_count);
  
  // fill the Jacobian matrix with the finite difference approximation to the
  // derivative
  std::string parm_name;
  tracefile << " Parameters of the simulation number " << last_genuine 
            << " are the current position.\n The results of the following"
            << " simulations are used for numerical derivation:\n";

  for (int i=0; i<num_variable_params; i++)
    tracefile << "  Simulation number " << varPar_history_map[i] 
              << " for the derivation with respect to " 
              << variable_params[i].name << ".\n";

  double prefac;
  if (current_params.parms[0].defined("CONST_PREFAC") 
               || current_params.parms[0].defined("POS_PREFAC")) {
    prefac = alps::evaluate<double>
               (current_params.parms[0].value_or_default("POS_PREFAC",0.0));
    prefac *= prefac;
    prefac += alps::evaluate<double>
               (current_params.parms[0].value_or_default("CONST_PREFAC",0.0));
  } else
    prefac = 1.0;

  double r_here, r_there, r_delta, param_here, param_there, param_delta;
  int col;
  int sim_index;
  for (col=0; col<num_variable_params; col++) {
    sim_index = varPar_history_map[col];
    parm_name = variable_params[col].name;
    param_here = (param_history[last_genuine].parms)[0][parm_name];
    param_there = (param_history[sim_index].parms)[0][parm_name];
    param_delta = param_there - param_here;
    if (param_delta == 0.0) {
      std::cerr << "cannot compute gradient with no difference in the "
                << "parameter values\n";
      boost::throw_exception(std::runtime_error(
               "invalid parameter values for gradient approximation"));
    }
    for (int row=0; row<result_count; row++) {
      r_here = residuum_history[last_genuine].res_vals[row].mean;
      r_there = residuum_history[sim_index].res_vals[row].mean;

      r_delta = r_there - r_here;
      jacobian(row,col) = prefac*r_delta/param_delta;
    }
  }

  // the rest of the Jacobian matrix is filled with the analytic derivatives
  if ((col<variable_params.size())
               && (variable_params[col].name == "FIT_CONSTANT")) {
    for (int row=0; row<result_count; row++)
      jacobian(row,col) = -1.0;
    col++;
  }

  alps::scheduler::ResultsType eff_res;
  if ((col<variable_params.size())
               && (variable_params[col].name == "CONST_PREFAC")) {
    for (int row=0; row<result_count; row++)
      jacobian(row,col) = - sim_results[last_genuine][row].mean;
    col++;
  }

  if ((col<variable_params.size())
               && (variable_params[col].name == "POS_PREFAC")) {
    double g = (param_history[last_genuine].parms)[0]["POS_PREFAC"];
    for (int row=0; row<result_count; row++)
      jacobian(row,col) = - 2.0*g*sim_results[last_genuine][row].mean;
    col++;
  }

  if ((col<variable_params.size())
               && (variable_params[col].name == "C_OVER_T")) {
    for (int row=0; row<result_count; row++)
      jacobian(row,col) = - 1.0/sim_results[last_genuine][row].T;
    col++;
  }

  double C, T, T_C;
  if ((col<variable_params.size())
               && (variable_params[col].name == "CURIE_CONST")) {
    C = (param_history[last_genuine].parms)[0]["CURIE_CONST"];
    T_C = (param_history[last_genuine].parms)[0]["CURIE_TEMP"];
    for (int row=0; row<result_count; row++) {
      T = (param_history[last_genuine].parms)[row]["T"];
      // derivation after C ("CURIE_CONST")
      jacobian(row,col) = -2.0*C/(T-T_C);
    }
    col++;
  }

  if ((col<variable_params.size())
               && (variable_params[col].name == "CURIE_TEMP")) {
    C = (param_history[last_genuine].parms)[0]["CURIE_CONST"];
    T_C = (param_history[last_genuine].parms)[0]["CURIE_TEMP"];
    for (int row=0; row<result_count; row++) {
      T = (param_history[last_genuine].parms)[row]["T"];
       // derivation after T_C ("CURIE_TEMP")
      jacobian(row,col) = SQUARE(C/(T-T_C));
    }
    col++;
  }

  // possibly more to come for other analytically differentiable parts of the
  // model

  // by now, all columns should have been filled
  if (col != col_count)
    std::cerr << "WARNNG : something went wrong while filling the jacobian "
              << "matrix.\n         col = " << col << ", col_count = " 
              << col_count << "\n\n";

  if (use_levenberg_marquart) {
    double lambda = determine_lambda(residuum_history[last_genuine].res_size);

    for (int row=result_count; row<row_count; row++) {
      for (int col=0; col<col_count; col++)
        jacobian(row,col) = 0.0;
      jacobian(row, row-result_count) = lambda;
    }
  }

  // compute the absolute value of the resiuum, put it into the right hand
  // side vector
  boost::numeric::ublas::vector<double> r_vec(row_count);
  int row;
  for (row = 0; row<result_count; row++)
    r_vec(row) = residuum_history[last_genuine].res_vals[row].mean;

  // if use_levenberg_marquart : fill the rest of the vector with 0 
  // else, row_count == result_count, i.e. nothing is done here
  for (; row<row_count; row++)
    r_vec(row) = 0.0;

  tracefile << "\nThe Jacobian matrix:\n ";
  for (int row=0; row<row_count; row++) {
    for (int col=0; col<col_count; col++) 
      tracefile << jacobian(row,col) << " ";
    tracefile << "\n ";
  }

  // solve the equation
  boost::numeric::ublas::vector<double> h(col_count);
  double cond;
  cond = solve_llsp(jacobian,r_vec,h);
  last_step.clear();
  varPar_history_map.clear();
  tracefile << "\nThe condition number of the matrix was " << cond << "\n";

  tracefile << "\nThe right hand side vector (residuum):\n";
  for (int i=0; i<row_count; i++)
    tracefile << " " << r_vec(i) << "\n";

  tracefile << "\nValues for the step:\n ";
  for (int i=0; i<col_count; i++) {
    tracefile << h(i) << " ";
    last_step.push_back(h(i));
  }
  tracefile << "\n";
}

/**
 * Add the current step size to the parameter list.
 *
 * @param plist A list parameters, the step is added to every single element
 *              of the list and stored into next_parameters.
 */
void EstGradFitter::do_step(alps::ParameterList& plist)
{
  for (int i=0; i<variable_params.size(); i++) {
    std::string param_name = variable_params[i].name;
    double new_param = (double)plist[0][param_name] - (double)last_step[i];
    if ( ((variable_params[i].limit.validmax) 
               && (new_param > variable_params[i].limit.max)) ||
          ((variable_params[i].limit.validmin)
               && (new_param < variable_params[i].limit.min))) {
      std::cerr << "Parameter " << param_name << " should be set to "
                << new_param << ".\nThis value lies outside the given "
                << "interval (" << variable_params[i].limit.min << ","
                << variable_params[i].limit.max << ") (validity:("
                << std::boolalpha << variable_params[i].limit.validmin
                << "," << std::boolalpha << variable_params[i].limit.validmax
                << ")).\nThe simulation is ended.\n";
      boost::throw_exception(std::runtime_error("value out of range"));
    } else
      plist = update_parameters(param_name, new_param, plist);
  }
}

/**
 * Write the current state of the fitter into an XML file.
 */
alps::oxstream& EstGradFitter::write_to_xml(alps::oxstream& oxs) const
{
  BaseIncrFitter::write_to_xml(oxs);

  // write special parameters of this type of fitter
  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","enlarge_trigger_factor")
      << alps::no_linebreak << enlarge_trigger_factor 
      << alps::end_tag("FITTING_PARAMETER");
  oxs << alps::start_tag("FITTING_PARAMETER") 
      << alps::attribute("name","last_base_increment")
      << alps::no_linebreak << last_base_increment 
      << alps::end_tag("FITTING_PARAMETER");
  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","enlarged_basis")
      << alps::no_linebreak << enlarged_basis
      << alps::end_tag("FITTING_PARAMETER");

  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","decrease_step_counter")
      << alps::no_linebreak << decrease_step_counter
      << alps::end_tag("FITTING_PARAMETER");

  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","stepsize_prefactor")
      << alps::no_linebreak << stepsize_prefactor
      << alps::end_tag("FITTING_PARAMETER");

  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","res_to_underbid")
      << alps::no_linebreak << res_to_underbid 
      << alps::end_tag("FITTING_PARAMETER");

  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","res_to_underbid_index")
      << alps::no_linebreak << res_to_underbid_index 
      << alps::end_tag("FITTING_PARAMETER");

  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","have_tightened")
      << alps::no_linebreak << have_tightened 
      << alps::end_tag("FITTING_PARAMETER");

  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","num_variable_params")
      << alps::no_linebreak << num_variable_params
      << alps::end_tag("FITTING_PARAMETER");

  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","variable_params");
//  int i=0;
  for (int i=0; i<variable_params.size(); i++)
    oxs << alps::start_tag("VECTOR_ELEMENT") << alps::no_linebreak
        << variable_params[i] << alps::end_tag("VECTOR_ELEMENT");
  oxs  << alps::end_tag("FITTING_PARAMETER");

  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","last_step");
  for (int i=0; i<last_step.size(); i++)
    oxs << alps::start_tag("VECTOR_ELEMENT") << alps::no_linebreak
        << last_step[i] << alps::end_tag("VECTOR_ELEMENT");
  oxs  << alps::end_tag("FITTING_PARAMETER");

  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","best_residui_size");
  for (int i=0; i<3; i++)
    oxs << alps::start_tag("VECTOR_ELEMENT") << alps::no_linebreak
        << best_residui_size[i] << alps::end_tag("VECTOR_ELEMENT");
  oxs  << alps::end_tag("FITTING_PARAMETER");

  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","thightened_indices");
  for (int i=0; i<thightened_indices.size(); i++)
    oxs << alps::start_tag("VECTOR_ELEMENT") << alps::no_linebreak 
        << thightened_indices[i] << alps::end_tag("VECTOR_ELEMENT");
  oxs  << alps::end_tag("FITTING_PARAMETER");

  return oxs;
}
