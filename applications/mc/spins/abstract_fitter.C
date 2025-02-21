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

#include "fitter.h"
#include "fitting_scheduler.h"

alps::oxstream& operator<< (alps::oxstream& oxs, const FitterParamList list) {
  oxs << list.parms;
  oxs << alps::start_tag("TYPE") << alps::attribute("name", "is_genuine")
      << alps::no_linebreak << list.is_genuine << alps::end_tag("TYPE");
  oxs << alps::start_tag("TYPE") << alps::attribute("name", "type")
      << alps::no_linebreak << list.type << alps::end_tag("TYPE");
  return oxs;
}

std::ifstream& operator>> (std::ifstream& is, FitterParamList& list) {
  using namespace alps;
  alps::XMLTag tag = alps::parse_tag(is,true);
  if (tag.name != "PARAMETERLIST") {
    std::cerr << "An error occured while parsing the checkpoint file\n";
    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  }
  tag = alps::parse_tag(is,true);
  while (tag.name != "/PARAMETERLIST") {
    alps::Parameters param;
    param.read_xml(tag,is,true);
    list.parms.push_back(param);
    tag = alps::parse_tag(is,true);
  }
  tag = alps::parse_tag(is,true);
  if (tag.name == "TYPE")
    list.is_genuine = (alps::parse_content(is) == "1");
  tag = alps::parse_tag(is,true);
  tag = alps::parse_tag(is,true);
  if (tag.name == "TYPE")
    list.type = ((alps::parse_content(is)).c_str())[0];
  return is;
}
   
alps::oxstream& operator<< (alps::oxstream& oxs, const ResidualType residual) {
  oxs << alps::start_tag("RESIDUAL");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","T") 
      << alps::no_linebreak << residual.T << alps::end_tag("ELEMENT");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","mean")
      << alps::no_linebreak << residual.mean << alps::end_tag("ELEMENT");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","error")
      << alps::no_linebreak << residual.error << alps::end_tag("ELEMENT");
  oxs << alps::end_tag("RESIDUAL");
  return oxs;
}

std::ifstream& operator>> (std::ifstream& is, ResidualType& residual) {
  alps::XMLTag tag = alps::parse_tag(is,true);
  if (tag.name != "RESIDUAL") {
    std::cerr << "An error occured while parsing the checkpoint file\n";
    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  }
  // read temperature
  tag = alps::parse_tag(is,true);
  residual.T = alps::evaluate<double>(alps::parse_content(is));
  tag = alps::parse_tag(is,true);
  // read mean
  tag = alps::parse_tag(is,true);
  residual.mean =  alps::evaluate<double>(alps::parse_content(is));
  tag = alps::parse_tag(is,true);
  // read error
  tag = alps::parse_tag(is,true);
  residual.mean =  alps::evaluate<double>(alps::parse_content(is));
  tag = alps::parse_tag(is,true);
  return is;
}
   
alps::oxstream& operator<< (alps::oxstream& oxs, const ResHistElementType res) {
  oxs << alps::start_tag("RESIDUUM_HISTORY_ELEMENT");
  for (int i=0; i<res.res_vals.size(); i++)
    oxs << res.res_vals[i];
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name", "res_size")
      << alps::no_linebreak << res.res_size << alps::end_tag("ELEMENT");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name", "avg_res_size")
      << alps::no_linebreak << res.avg_res_size << alps::end_tag("ELEMENT");
  oxs << alps::end_tag("RESIDUUM_HISTORY_ELEMENT");
  return oxs;
}

std::ifstream& operator>> (std::ifstream& is, ResHistElementType& reshist) {
  alps::XMLTag tag = alps::parse_tag(is,true);
  if (tag.name != "RESIDUAL") {
    std::cerr << "An error occured while parsing the checkpoint file\n";
    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  }
  while (tag.name != "ELEMENT") {
    ResidualType rt;
    tag = alps::parse_tag(is,true);
    rt.T = alps::evaluate<double>(alps::parse_content(is));
    tag = alps::parse_tag(is,true);
    tag = alps::parse_tag(is,true);
    rt.mean = alps::evaluate<double>(alps::parse_content(is));
    tag = alps::parse_tag(is,true);
    tag = alps::parse_tag(is,true);
    rt.error = alps::evaluate<double>(alps::parse_content(is));
    reshist.res_vals.push_back(rt);
    tag = alps::parse_tag(is,true);
    tag = alps::parse_tag(is,true);
    tag = alps::parse_tag(is,true);
  }
  // read res_size
  reshist.res_size = alps::evaluate<double>(alps::parse_content(is));
  tag = alps::parse_tag(is,true);
  // read avg_res_size
  tag = alps::parse_tag(is,true);
  reshist.avg_res_size = alps::evaluate<double>(alps::parse_content(is));
  tag = alps::parse_tag(is,true);
  return is;
}

alps::oxstream& operator<< (alps::oxstream& oxs,
               const alps::scheduler::ResultType res) {
  oxs << alps::start_tag("RESULT");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","T")
      << alps::no_linebreak << res.T << alps::end_tag("ELEMENT");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","mean")
      << alps::no_linebreak << res.mean << alps::end_tag("ELEMENT");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","error")
      << alps::no_linebreak << res.error << alps::end_tag("ELEMENT");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","count")
      << alps::no_linebreak << res.count << alps::end_tag("ELEMENT");
  oxs << alps::end_tag("RESULT");
  return oxs;
}

std::ifstream& operator>> (std::ifstream& is, alps::scheduler::ResultType& res) {
  alps::XMLTag tag =  alps::parse_tag(is,true);
  if (tag.name != "RESULT") {
    std::cerr << "An error occured while parsing the checkpoint file\n";
    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  }
  tag = alps::parse_tag(is,true);
  res.T =  alps::evaluate<double>(alps::parse_content(is));
  tag = alps::parse_tag(is,true);
  tag = alps::parse_tag(is,true);
  res.mean = alps::evaluate<double>(alps::parse_content(is));
  tag = alps::parse_tag(is,true);
  tag = alps::parse_tag(is,true);
  res.error = alps::evaluate<double>(alps::parse_content(is));
  tag = alps::parse_tag(is,true);
  tag = alps::parse_tag(is,true);
  res.count = alps::evaluate<double>(alps::parse_content(is));
  tag = alps::parse_tag(is,true);
  return is;
}

alps::oxstream& operator<< (alps::oxstream& oxs, 
               const alps::scheduler::ResultsType results)
{
  oxs << alps::start_tag("RESULTS");
  for (int i=0; i<results.size(); i++)
    oxs << results[i];
  oxs << alps::end_tag("RESULTS");
  return oxs;
}

std::ifstream& operator>> (std::ifstream& is,
               alps::scheduler::ResultsType& results) {
  alps::XMLTag tag =  alps::parse_tag(is,true);
  if (tag.name != "RESULTS") {
    std::cerr << "An error occured while parsing the checkpoint file\n";
    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  }
  tag = alps::parse_tag(is,true);
  while (tag.name != "/RESULTS") {
    alps::scheduler::ResultType res;
    tag = alps::parse_tag(is,true);
    res.T = alps::evaluate<double>(alps::parse_content(is));
    tag = alps::parse_tag(is,true);
    tag = alps::parse_tag(is,true);
    res.mean = alps::evaluate<double>(alps::parse_content(is));
    tag = alps::parse_tag(is,true);
    tag = alps::parse_tag(is,true);
    res.error = alps::evaluate<double>(alps::parse_content(is));
    tag = alps::parse_tag(is,true);
    tag = alps::parse_tag(is,true);
    res.count = alps::evaluate<double>(alps::parse_content(is));
    tag = alps::parse_tag(is,true);
    tag = alps::parse_tag(is,true);
    results.push_back(res);
    tag = alps::parse_tag(is,true);
  }
  return is;
}

alps::oxstream& operator<< (alps::oxstream& oxs, const Limit limit)
{
  oxs << alps::start_tag("LIMIT");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","min")
      << alps::no_linebreak << limit.min << alps::end_tag("ELEMENT");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","max")
      << alps::no_linebreak << limit.max << alps::end_tag("ELEMENT");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","validmin")
      << alps::no_linebreak << limit.validmin << alps::end_tag("ELEMENT");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","validmax")
      << alps::no_linebreak << limit.validmax << alps::end_tag("ELEMENT");
  oxs << alps::end_tag("LIMIT");
  return oxs;
}

std::ifstream& operator>> (std::ifstream& is, Limit& limit)
{
  alps::XMLTag tag = alps::parse_tag(is,true);
  if (tag.name != "LIMIT") {
    std::cerr << "An error occured while parsing the checkpoint file\n";
    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  }
  tag = alps::parse_tag(is,true);
  // read min
  limit.min = alps::evaluate<double>(alps::parse_content(is));
  tag = alps::parse_tag(is,true);
  tag = alps::parse_tag(is,true);
  // read max
  limit.max = alps::evaluate<double>(alps::parse_content(is));
  tag = alps::parse_tag(is,true);
  tag = alps::parse_tag(is,true);
  // read valid flags
  limit.validmin = (alps::parse_content(is) == "1");
  tag = alps::parse_tag(is,true);
  tag = alps::parse_tag(is,true);
  limit.validmax = (alps::parse_content(is) == "1");
  tag = alps::parse_tag(is,true);
  tag = alps::parse_tag(is,true);
  return is;
}

alps::oxstream& operator<< (alps::oxstream& oxs, const FittingParam fp)
{
  oxs << alps::start_tag("FITTING_PARAM");
  oxs << alps::start_tag("ELEMENT") << alps::attribute("name","name")
      << fp.name << alps::end_tag("ELEMENT");
  oxs << fp.limit;
  oxs << alps::end_tag("FITTING_PARAM");
  return oxs;
}

std::ifstream& operator>> (std::ifstream& is, FittingParam& fp)
{
  alps::XMLTag tag = alps::parse_tag(is,true);
  if (tag.name != "FITTING_PARAM") {
    std::cerr << "An error occured while parsing the checkpoint file\n";
    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  }
  tag = alps::parse_tag(is,true);
  fp.name = alps::parse_content(is);
  tag = alps::parse_tag(is,true);
  is >> fp.limit;
  return is;
} 

/**
 * Constructor for the absract class AbstractFitter.
 *
 * @param fit_file_name The name of the file containing the parameters for the
 *                      fitting
 */
AbstractFitter::AbstractFitter(const std::string& fit_file_name)
{
  // reset everything
  all_exp_res.clear();
  sim_results.clear();
  param_history.clear();
  current_params.parms.clear();
  next_parameters.clear();
  residuum_history.clear();
  last_genuine = -1;
  variable_params.clear();
  results_namebase = fit_file_name;

  // read the simulation parameter file
  alps::Parameters fit_param;
  {
    std::ifstream in(fit_file_name.c_str());
    in >> fit_param;
  }
  // handle the time limit
  long time_limit = static_cast<boost::int64_t>
               (fit_param.value_or_default("TIME_LIMIT",-1l));
  start_time = boost::posix_time::second_clock::local_time();
  end_time = start_time + boost::posix_time::seconds(time_limit);

  std::string fname = fit_file_name + "_LOG";
  int counter;
  char new_name[256];
  if (boost::filesystem::exists(fname)) {
    counter = 1;
    sprintf(new_name, "%s_restart%i_LOG", fit_file_name.c_str(), counter);
    while (boost::filesystem::exists(new_name)) {
      counter++;
      sprintf(new_name, "%s_restart%i_LOG", fit_file_name.c_str(), counter);
    }
    tracefile.open(new_name, std::fstream::out);
  } else
    tracefile.open(fname.c_str(), std::fstream::out);
                                                                                
  // check for the observable to be fitted for
  if (!fit_param.defined("FIT_OBSERVABLE")) {
    std::cerr << "No observable to fit defined, returning\n";
    boost::throw_exception(std::runtime_error("no observable to fit"));
  }
  fit_obs_name = fit_param["FIT_OBSERVABLE"];

  // check for variable parameters
  variable_params = find_variable_parms(fit_param);

  num_variable_params = variable_params.size();
  if (num_variable_params == 0) {
    std::cerr << "No variable parameters found.\n"
              << "Cannot do a fitting without variable parameters!!\n";
    boost::throw_exception(std::runtime_error("no variable parameters"));
  }

  // check for the target file name
  if (!fit_param.defined("TARGET_FILE_NAME")) {
    std::cerr << "No target file name is given, aborting\n";
    boost::throw_exception(std::runtime_error("no target file name"));
  }

  // read for the target file
  read_target(fit_param["TARGET_FILE_NAME"]);

  // read the initial parameter file, store the values into the future
  // parameters list
  if (!fit_param.defined("PARAMETER_FILE")) {
    // later : produce a default init file
    std::cerr << "No initial parameter file given, aborting\n";
    boost::throw_exception(std::runtime_error("no parameter file"));
  }
  filenamebase = fit_param["PARAMETER_FILE"];
  srand(time(NULL));

  // read the values for the stop criterion
  max_steps = fit_param.value_or_default("MAX_STEPS",50);
  max_avg_residuum = fit_param.value_or_default("MAX_AVG_RESIDUUM",0.0);
  // by default, there's no upper limit for the residuum
}

/**
 * Initialize the state of the fitting, reading parameter values and 
 * status data from the stream.
 *
 * @param is The stream containing the data.
 */
void AbstractFitter::initialize(std::ifstream& is) 
{
  // read parameters from the checkpoint file
  alps::XMLTag tag = alps::parse_tag(is, true);
  if (tag.type == alps::XMLTag::SINGLE)
    std::cout << "No parameter history in checkpoint file.\n";
  else {
    tag = alps::parse_tag(is, true);
    while (tag.name != "/PARAMETER_HISTORY") {
      FitterParamList fpl;
      is >> fpl;
      if (fpl.is_genuine)
        last_genuine = param_history.size();

      param_history.push_back(fpl);
      tag = alps::parse_tag(is, true);
      tag = alps::parse_tag(is, true);
      tag = alps::parse_tag(is, true);
    }
    std::cout << "Read " << param_history.size() << " previous parameters.\n";
  }
  
  step = param_history.size();
  std::cout << "Step set to " << step << "\n";

  tag = alps::parse_tag(is, true);
  if (tag.type == alps::XMLTag::SINGLE) {
    std::cerr << "Invalid checkpoint file: No current parameters.\n";
    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  } else {
    is >> current_params;
    tag = alps::parse_tag(is, true);
    std::cout << "read current parameters.\n";
  }
  tag = alps::parse_tag(is, true);
  tag = alps::parse_tag(is, true);

  next_parameters.clear();
  if (tag.type == alps::XMLTag::SINGLE)
    std::cout << "No next parameters defined in the checkpoint file.\n";
  else {
    tag = alps::parse_tag(is, true); // NEXT_PARAM opening tag
    while (tag.name != "/NEXT_PARAMETERS") {
      FitterParamList fpl;
      is >> fpl;
      next_parameters.push_back(fpl);

      tag = alps::parse_tag(is, true);
      tag = alps::parse_tag(is, true);
      tag = alps::parse_tag(is, true);
    }
    std::cout << "read " << next_parameters.size() << " next parameter sets.\n";
  }


  tag = alps::parse_tag(is, true);
  all_exp_res.clear();
  if (tag.type == alps::XMLTag::SINGLE) {
    std::cerr << "Invalid checkpoint file: No experimental results.\n";
    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  } else {
    is >> all_exp_res;
    std::cout << "read full set of experimental results.\n";
    tag = alps::parse_tag(is, true);
  }  

  tag = alps::parse_tag(is,true);
  if (tag.type == alps::XMLTag::SINGLE) {
    std::cerr << "No simulation results in checkpoint file.\n";
  } else {
    tag = alps::parse_tag(is, true);

    while (tag.name != "/SIMULATION_RESULTS") {
      alps::scheduler::ResultsType rt;
      is >> rt;
      sim_results.push_back(rt);
      tag = alps::parse_tag(is, true);
      tag = alps::parse_tag(is, true);
    }
    std::cout << "read " << sim_results.size() << " simulation results.\n";
  }

  tag = alps::parse_tag(is, true);
  if (tag.type == alps::XMLTag::SINGLE) {
    std::cerr << "No residuum history given.\n";
  } else {
    tag = alps::parse_tag(is,true); // starting tag for RES_HIST_ELEMENT
    while (tag.name != "/RESIDUUM_HISTORY") {
      ResHistElementType rhet;
      is >> rhet;
      residuum_history.push_back(rhet);
      tag = alps::parse_tag(is, true); // end RES.HIST.ELEMENT
      tag = alps::parse_tag(is, true); // start " or end RESIDUUM_HISTORY
    }
    std::cout << "Read " << residuum_history.size() << " residui sets.\n";
  }

  tag = alps::parse_tag(is, true);
  if (tag.name == "FITTING_PARAMETER") {
    filenamebase = alps::parse_content(is);
    std::cout << "Filenamebase set to " << filenamebase << "\n";
    tag = alps::parse_tag(is,true);
  }
  tag = alps::parse_tag(is,true);

  if (tag.name == "FITTING_PARAMETER") {
    current_filename = alps::parse_content(is);
    tag = alps::parse_tag(is,true);
    std::cerr << "Current job file name set to " << current_filename << "\n";
  }
}

/**
 * Initialize the class if no checkpoint file is available. Nothing is to be
 * done.
 */
void AbstractFitter::initialize() 
{
}

/**
 * Read in the measurements of a physical experience and stores the values.
 * In the file, every line describes an measurements. For every measurement, the
 * name of the observable, the temperature, the mean and optionally the error is
 * given (in this order). If no error is given in the file, the error is set to
 * 0 (zero).
 * 
 * @param file_name the name of the input file.
 */
void AbstractFitter::read_target(const std::string file_name)
{ 
  FILE* myFile;
  alps::scheduler::ResultType cur;
  myFile = fopen(file_name.c_str(),"r"); 
  if (NULL == myFile) {
    std::cerr << "File could not be opened, returning ...\n";
    boost::throw_exception(std::runtime_error("could not read file"));
  }
  char currentLine[256];
  std::string name_string;
  std::string T_string, mean_string, error_string;
  int pos=0;
  fgets(currentLine, 256, myFile);
  while (!feof(myFile)) {
    pos = 0;
    while ((pos<256) && (currentLine[pos] == ' ')&&(currentLine[pos] != '\n'))
      pos++;
    while ((pos<256) && (currentLine[pos] != ' ')&&(currentLine[pos] != '\n')) {
      T_string.append(1,currentLine[pos]);
      pos++;
    }
    while ((pos<256) && (currentLine[pos] == ' ')&&(currentLine[pos] != '\n')) 
       pos++;
    while ((pos<256) && (currentLine[pos] != ' ')&&(currentLine[pos] != '\n')) {
      mean_string.append(1,currentLine[pos]);
      pos++;
    }
    while ((pos<256) && (currentLine[pos] == ' ')&&(currentLine[pos] != '\n'))
      pos++;
    while ((pos<256) && (currentLine[pos] != ' ')&&(currentLine[pos] != '\n')) {
      error_string.append(1,currentLine[pos]);
      pos++;
    }
    while ((pos<256) && (currentLine[pos] == ' ')&&(currentLine[pos] != '\n'))
      pos++;

    // parse the different strings
    cur.T = alps::evaluate<double>(T_string);
    cur.mean = alps::evaluate<double>(mean_string);
    cur.error = alps::evaluate<double>(error_string);
    if (cur.T != 0.0) {
      // if T == 0.0, we get a division by 0 when calculating beta, which 
      // returns nan and returns errors when the simulation results are handled.
      all_exp_res.push_back(cur);
    }

    // prepare for the next line
    T_string.clear();
    mean_string.clear();
    error_string.clear();
    fgets(currentLine,256,myFile);
  }
  fclose(myFile);
}

/**
 * Return a set containing the data to plot the exponential results.
 */
alps::plot::Set<double> AbstractFitter::get_exp_results_plot_set() const
{
  alps::plot::Set<double> res(alps::plot::xydy);
  for (int i=0; i<all_exp_res.size(); i++) {
    res << all_exp_res[i].T;
    res << all_exp_res[i].mean;
    res << all_exp_res[i].error;
  }
  // set the label of the plot set
  res << "Experimental data";
  return res;
}
   
/**
 * Returns a set containing the data to produce a plot of the results of the
 * last simulation.
 */
alps::plot::Set<double> AbstractFitter::get_sim_results_plot_set() const
{
  alps::plot::Set<double> res(alps::plot::xydy);
  res.clear();
  if (sim_results.empty())
    return res;

  unsigned int pos = sim_results.size() - 1;
  alps::scheduler::ResultsType eff_sim_res;
  eff_sim_res = get_eff_sim_res(sim_results[pos],current_params.parms);

  for (int i=0; i<eff_sim_res.size(); i++) {
    res << eff_sim_res[i].T;
    res << eff_sim_res[i].mean;
    res << eff_sim_res[i].error;
  }
  // set the label of the plot set
  res << "Simulation data";
  return res;
}

/**
 * Returns a set containing the data to generate a plot showing the error of
 * the experimental results.
 */
alps::plot::Set<double> AbstractFitter::get_exp_error_plot_set() const
{
  alps::plot::Set<double> res(alps::plot::xy);
  res.clear();

  for (int i =0; i<all_exp_res.size(); i++) {
    res << all_exp_res[i].T;
    res << all_exp_res[i].error;
  }
  // label
  res << "Error of experimental data";
  return res;
}

/**
 * Returns a set containing the data to produce a plot showing the residual
 * values.
 */
alps::plot::Set<double> AbstractFitter::get_res_abs_plot_set() const
{ 
  alps::plot::Set<double> res(alps::plot::xydy);
  res.clear();
  if (residuum_history.empty())
    return res;
  
  std::vector<ResidualType> res_vals = residuum_history.back().res_vals;
  for (int i =0; i<res_vals.size(); i++) {
    res << res_vals[i].T;
    res << res_vals[i].mean;
    res << res_vals[i].error;
  }
  // set label
  res << "Residuum size";
  return res;
}

/**
 * Returns a set with the percentual error of the experimental data.
 */
alps::plot::Set<double> AbstractFitter::get_exp_percent_plot_set() const
{
  alps::plot::Set<double> res(alps::plot::xy);
  res.clear();
  for (unsigned int i=0; i<all_exp_res.size(); i++) {
    res << all_exp_res[i].T;
    if (all_exp_res[i].mean == 0.0)
      res << 0.0;
    else
      res << 100*all_exp_res[i].error/all_exp_res[i].mean;
  }
  // label
  res << "Percentual error in measurement data";
  return res;
}

/**
 * Returns a set containing the data for a plot of the size of the residual
 * in percent of the mean value of the experimental results.
 */
alps::plot::Set<double> AbstractFitter::get_res_percent_plot_set() const
{
  alps::plot::Set<double> res(alps::plot::xydy);
  res.clear();
  if (residuum_history.empty()) 
    return res;

  std::vector<ResidualType> res_vals = residuum_history.back().res_vals;
  for (unsigned int i=0; i<all_exp_res.size(); i++) 
    if (all_exp_res[i].mean != 0) {
      res << all_exp_res[i].T;
      res << 100.0*res_vals[i].mean/all_exp_res[i].mean;
      res << 100.0*res_vals[i].error/all_exp_res[i].mean;
    }
  res << "Persental size of Residuum";
  return res;
}

/**
 * Takes the parameters and writes the job and task files.
 * The name of the job file is <basename>.step<s>.in.xml, the task files will
 * have the names <basename>.step<s>.task<t>.in.xml, where <basename> is the
 * name of the initial parameter file given in the simulation file, s is the
 * step count and t is the number of the task.
 */
std::string AbstractFitter::print_to_file(alps::ParameterList list) const
{
  if (list.size() == 0) return std::string();

  std::string basename = filenamebase+".step"
                +boost::lexical_cast<std::string,int>(step);
  std::cout << "Writing job file " << basename << ".in.xml\n";

  int bits = 31;
  for (int n = 1; n < list.size(); n<<=1, --bits);

  boost::uint32_t baseseed;
  if (list[0].defined("BASESEED")) {
	  baseseed = boost::lexical_cast<boost::uint32_t>(list[0]["BASESEED"]);
  } else {
    baseseed =
      boost::posix_time::microsec_clock::local_time().time_of_day().
      total_microseconds() << 24;
  }
  baseseed &= ((1<<30)|((1<<30)-1));

  alps::oxstream out(boost::filesystem::path(basename+".in.xml"));
  out << alps::header("UTF-8")
      << alps::stylesheet(alps::xslt_path("ALPS.xsl"))
      << alps::start_tag("JOB")
      << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
      << alps::attribute("xsi:noNamespaceSchemaLocation",
                         "http://xml.comp-phys.org/2003/8/job.xsd")
      << alps::start_tag("OUTPUT")
      << alps::attribute("file", basename+".out.xml")
      << alps::end_tag("OUTPUT");

  for (int i = 0; i < list.size(); ++i) {
    std::string taskname =
      basename+".task"+boost::lexical_cast<std::string,int>(i+1);

    if (!list[i].defined("SEED")) {
      unsigned seed = baseseed ^ (i << bits);
      list[i]["SEED"] = seed;
    }

    out << alps::start_tag("TASK") << alps::attribute("status","new")
        << alps::start_tag("INPUT")
        << alps::attribute("file", taskname + ".in.xml")
        << alps::end_tag("INPUT")
        << alps::start_tag("OUTPUT")
        << alps::attribute("file", taskname + ".out.xml")
        << alps::end_tag("OUTPUT")
        << alps::end_tag("TASK");

    alps::oxstream task(boost::filesystem::path(taskname+".in.xml"));
    task << alps::header("UTF-8")
         << alps::stylesheet(alps::xslt_path("ALPS.xsl"));
    task << alps::start_tag("SIMULATION")
         << alps::xml_namespace("xsi",
                                "http://www.w3.org/2001/XMLSchema-instance")
         << alps::attribute("xsi:noNamespaceSchemaLocation",
                            "http://xml.comp-phys.org/2002/10/QMCXML.xsd");
    task << list[i];
    task << alps::end_tag("SIMULATION");
  }
  out << alps::end_tag("JOB");
  return basename.append(".in.xml");
}

/**
 * Returns the simulation results, taking into account all model parameters.
 *
 * @param sim_res The pure simulation results
 * @param param   A list of parameters defining the fitting model.
 */
alps::scheduler::ResultsType AbstractFitter::get_eff_sim_res(
               const alps::scheduler::ResultsType sim_res,
               const alps::ParameterList param) const
{
  alps::scheduler::ResultsType res;
  unsigned int length = sim_res.size();
  res.resize(length);
  double fit_constant = alps::evaluate<double>
               (param[0].value_or_default("FIT_CONSTANT",0.0));

  double prefac;
  if (param[0].defined("CONST_PREFAC") || param[0].defined("POS_PREFAC")) {
    prefac = alps::evaluate<double>
               (param[0].value_or_default("POS_PREFAC",0.0));
    prefac *= prefac;
    prefac += alps::evaluate<double>
               (param[0].value_or_default("CONST_PREFAC",0.0));
  } else
    prefac = 1.0;

  double c_over_t = alps::evaluate<double>
               (param[0].value_or_default("C_OVER_T",0.0));
  double curie_const = alps::evaluate<double>
               (param[0].value_or_default("CURIE_CONST",0.0));
  curie_const *= curie_const;
  double curie_temp = alps::evaluate<double>
               (param[0].value_or_default("CURIE_TEMP",0.0));

  for (int i=0; i<length; i++) {
    res[i].T = sim_res[i].T;
    res[i].mean = fit_constant + prefac*sim_res[i].mean
               + curie_const/(sim_res[i].T - curie_temp);
    if (c_over_t != 0.0)
      res[i].mean += c_over_t/res[i].T;
    res[i].error = prefac*sim_res[i].error;
  }
  return res;
}

/** 
 * Adds a simulation result ("summary") to the internal collection of 
 * results.
 * 
 * @param sim_res the result of the simulation
 */
void AbstractFitter::add_sim_result(const alps::scheduler::ResultsType sim_res)
{
  sim_results.push_back(sim_res);

  // update residuum history
  ResHistElementType cur_res;
  cur_res.res_vals = get_residuum(step,all_exp_res);
  cur_res.res_size = get_residuum_size(cur_res.res_vals);
  cur_res.avg_res_size = cur_res.res_size/cur_res.res_vals.size();
  residuum_history.push_back(cur_res);
 
  // simulation results and parameters only make sence together
  // the parameters are also put into history here
  param_history.push_back(current_params);
}

/**
 * Returns the residuum of the simulation number pos with respect to the 
 * experimental data given as second parameter.
 * 
 * @param pos The number of the simulation the residuum has to be computed
 *      of.
 * @param exp_res The experimental result used to compute the residuum
 */
std::vector<ResidualType> AbstractFitter::get_residuum(unsigned int pos,
               const alps::scheduler::ResultsType exp_res) const
{
  if (pos >= sim_results.size()) {
    boost::throw_exception(std::runtime_error(
               "invalid index to compute residuum"));
  }
  if (exp_res.empty()) {
    std::cerr << "No experimental data available, thus no residuum can be "
              << "computed.\ncomputation is aborted\n";
    boost::throw_exception(std::runtime_error("no experimental data given"));
  }

  alps::scheduler::ResultsType eff_res = 
               get_eff_sim_res(sim_results[pos],current_params.parms);
  unsigned int size = eff_res.size();
  std::vector<ResidualType> res;
  res.resize(size);

  for (unsigned int i=0; i<size; i++) {
    if (exp_res[i].T != eff_res[i].T) {
      std::cerr << "\ncannot compute residuum for measurements with different "
                << "temperatures!!"
                << "\nTemperature of oberved data is " << exp_res[i].T
                << "\nTemperature of simulation data is " << eff_res[i].T
                << "\nNo residuum will be computed\n\n";
    } else
      res[i].T = exp_res[i].T;
      res[i].mean = exp_res[i].mean - eff_res[i].mean;
      if (exp_res[i].error == 0)
        res[i].error = eff_res[i].error;
      else
        res[i].error = sqrt(SQUARE(eff_res[i].error)
                           +SQUARE(exp_res[i].error));
  }
  return res;
}

/**
 * Compute the lenght of the given residual vector
 *
 * @param residuum The residual vector.
 */
double AbstractFitter::get_residuum_size(
               const std::vector<ResidualType> residuum) const
{
  double res = 0.0;
  unsigned int size = residuum.size();
  if (size==0) return 0.0;

  for (int i=0; i<size; i++)
    res += square(residuum[i].mean);
  return sqrt(res);
}

/**
 * Compute the length of the given residual vector divided by the number of
 * components of the vector
 *
 * @param residuum The residual vector which size per component is to be 
 *                 computed.
 */
double AbstractFitter::get_avg_res_size(
               const std::vector<ResidualType> residuum) const
{
  unsigned int size = residuum.size();
  if (size==0) return 0.0;
  return get_residuum_size(residuum)/size;
}

/**
 * Returns the length of the last residual.
 */
double AbstractFitter::get_last_res_size() const
{ return (residuum_history.back()).res_size; }

/**
 * Produce a new job file along with the corresponding task files.
 * If necessary, new parameter sets are generated, the concrete algorithms
 * therefore are to be implemented in the subclasses.
 */
void AbstractFitter::make_new_jobfile()
{
  // make_new_jobfile is the beginning of a new fitting step
  step++;
  if (current_params.is_genuine)
    last_genuine = param_history.size() - 1;

  if (next_parameters.empty())
    find_new_parameters();

  current_params = next_parameters.front();
  next_parameters.pop_front();
  current_filename = print_to_file(current_params.parms);
}

/**
 * Get the name of the current job file.
 */
std::string AbstractFitter::get_current_filename() {
  return current_filename;
}

/**
 * Print out several results of the fitting into a series of files.
 */
void AbstractFitter::print_results() const
{
  std::string results_file = results_namebase + "_RESULTS";
  print_results(results_file);
}

/**
 * Print out several results of the fitting into a series of files.
 *
 * @param filename The base to be used for the file names.
 */
void AbstractFitter::print_results(const std::string filename) const
{
  if ((current_params.parms.size()==0) && (param_history.size() == 0))
    return;

  using namespace std;

  std::string mean_name, error_name, rel_name, count_name, res_name, 
               allres_name, ressign_name;
  std::fstream mean_stream, error_stream, rel_stream, count_stream, res_stream,
               allres_stream, ressign_stream;

  mean_name = filename + ".mean";
  error_name = filename + ".error";
  rel_name = filename + ".rel_error";
  count_name = filename + ".count";
  res_name = filename + ".residual";
  allres_name = filename + ".all_res";
  ressign_name = filename + ".res_sign";
  mean_stream.open(mean_name.c_str(), std::fstream::out);
  error_stream.open(error_name.c_str(), std::fstream::out);
  rel_stream.open(rel_name.c_str(), std::fstream::out);
  count_stream.open(count_name.c_str(), std::fstream::out);
  res_stream.open(res_name.c_str(), std::fstream::out);
  allres_stream.open(allres_name.c_str(), std::fstream::out);
  ressign_stream.open(ressign_name.c_str(), std::fstream::out);

  mean_stream << "results of the simulation: Values for the observable "
              << fit_obs_name << "\n==========================\n\n";
  error_stream << "results of the simulation: Values for the observable "
               << fit_obs_name << "\n==========================\n\n";
  rel_stream << "results of the simulation: Values for the observable "
             << fit_obs_name << "\n==========================\n\n";
  count_stream << "results of the simulation: Values for the observable "
               << fit_obs_name << "\n==========================\n\n";
  res_stream << "residui for the fitting : obervable " << fit_obs_name
             << "\n=======================\n\n";
  allres_stream << "residui for the fitting : obervable " << fit_obs_name
             << "\n=======================\n\n";
  ressign_stream << "error/mean of the residui for the fitting : obervable " 
                 << fit_obs_name << "\n=========================\n\n";

  std::cout << "Printing results of the simulation\n";
  int sim_count = sim_results.size();
  int T_count;
  double avg_res_size, least_avg_res_size;
  if (sim_count == 0)
    T_count = 0;
  else {
    T_count = sim_results[sim_count-1].size();
    least_avg_res_size = residuum_history[0].res_size;
  }
  unsigned int var_par_count = variable_params.size();
  alps::ParameterList print_param;

  mean_stream << "step   ";
  error_stream << "step   ";
  rel_stream << "step   ";
  count_stream << "step   ";
  res_stream << "step   ";
  allres_stream << "step   ";
  ressign_stream << "step   ";

  // print out the names of all variable parameters
  for (unsigned int i=0; i<var_par_count; i++) {
    mean_stream  << std::setw(13) << variable_params[i].name << " ";
    error_stream << std::setw(13) << variable_params[i].name << " ";
    rel_stream << std::setw(13) << variable_params[i].name << " ";
    count_stream << std::setw(13) << variable_params[i].name << " ";
    res_stream << std::setw(13) << variable_params[i].name << " ";
    allres_stream << std::setw(13) << variable_params[i].name << " ";
    ressign_stream << std::setw(13) << variable_params[i].name << " ";
  }

  res_stream << setw(12) << "|r|/n";
  res_stream << setw(13) << "n";
  allres_stream << setw(12) << "|r|/n";

  // print out all the temperature values
  for (unsigned int i=0; i<T_count; i++) {
    mean_stream << std::setw(12) << std::showpoint 
                << sim_results[sim_count-1][i].T << " ";
    error_stream << std::setw(12) << std::showpoint 
                 << sim_results[sim_count-1][i].T << " ";
    rel_stream << std::setw(12) << std::showpoint 
               << sim_results[sim_count-1][i].T << " ";
    count_stream << std::setw(12) << std::showpoint
                 << sim_results[sim_count-1][i].T << " ";
    if (i<9) 
      allres_stream << std::setw(11) << "r(" << i+1 << ") ";
    else
      allres_stream << std::setw(10) << "r(" << i+1 << ") ";
    if (i<9)
      ressign_stream << std::setw(11) << "r(" << i+1 << ") ";
    else
      ressign_stream << std::setw(10) << "r(" << i+1 << ") ";

  }
  mean_stream << "\n";
  error_stream << "\n";
  rel_stream << "\n";
  count_stream << "\n";
  res_stream << "\n";
  allres_stream << "\n";
  ressign_stream << std::setw(10) << "count\n";

  unsigned int least_avg_res_index = 0;
  for (int sim=0; sim < sim_count; sim++) {
    char type = param_history[sim].type;

    // print step number
    mean_stream << std::setw(4) << sim << " " << type << " ";
    error_stream << std::setw(4) << sim << " " << type << " ";
    rel_stream << std::setw(4) << sim << " " << type << " ";
    count_stream << std::setw(4) << sim << " " << type << " ";
    res_stream << std::setw(4) << sim << " " << type << " ";
    allres_stream << std::setw(4) << sim << " " << type << " ";
    ressign_stream << std::setw(4) << sim << " " << type << " ";

    print_param = param_history[sim].parms;
    double value;
    for (int i=0; i<var_par_count; i++) {
      value = (double)print_param[0][variable_params[i].name];
      if (value >=0.0) {
        mean_stream << "+";
        error_stream << "+";
        rel_stream << "+";
        count_stream << "+";
        res_stream << "+";
        allres_stream << "+";
        ressign_stream << "+";
      }

      mean_stream << setw(8) << scientific << value << " ";
      error_stream << setw(8) << scientific << value << " ";
      rel_stream << setw(8) << scientific << value << " ";
      count_stream << setw(8) << scientific << value << " ";
      res_stream << setw(8) << scientific << value << " ";
      allres_stream << setw(8) << scientific << value << " ";
      ressign_stream << setw(8) << scientific << value << " ";
    }

    avg_res_size = residuum_history[sim].avg_res_size;
    res_stream << std::setw(12) << std::scientific << avg_res_size << " ";
    allres_stream << std::setw(12) << std::scientific << avg_res_size << " ";

    if (avg_res_size < least_avg_res_size) {
      least_avg_res_size = avg_res_size;
      least_avg_res_index = sim;
    }
    
    // the number of temperatures may change in every simulation
    T_count = sim_results[sim].size();
    int erroneous_count = 0;

    for (int i=0; i<T_count; i++) {
      value = (double)sim_results[sim][i].mean;
      if (value >=0.0) mean_stream << "+";
      mean_stream << std::setw(7) << std::scientific << value << " ";

      value = (double)sim_results[sim][i].error;
      if (value >=0.0) error_stream << "+";
      error_stream << std::setw(7) << std::scientific << value << " ";

      value=(double)sim_results[sim][i].error/(double)sim_results[sim][i].mean;
      if (value >=0.0) rel_stream << "+";
      rel_stream << std::setw(7) << std::scientific << value << " ";

      value = (double)sim_results[sim][i].count;
      count_stream << std::setw(10) << std::scientific << value << " ";

      value = (double)residuum_history[sim].res_vals[i].mean;
      if (value >=0.0) allres_stream << "+";
      allres_stream << std::setw(7) << std::scientific << value << " ";

      value = residuum_history[sim].res_vals[i].error/value;
      value = std::abs(value);
      ressign_stream << std::setw(8) << std::scientific << value << " ";
      if (value > tolerated_error_factor) 
        erroneous_count++;

    }

    mean_stream << "\n";
    error_stream << "\n";
    rel_stream << "\n";
    count_stream << "\n";
    res_stream << std::setw(12) << residuum_history[sim].res_vals.size() <<"\n";
    allres_stream << "\n";
    ressign_stream << std::setw(12) << erroneous_count << " out of " 
                   << residuum_history[sim].res_vals.size() << "\n";
    erroneous_count = 0;
  }
  std::string fitter_name = get_fitter_name();
  mean_stream << "\nDescription of the fitter:\n" << fitter_name;
  error_stream << "\nDescription of the fitter:\n" << fitter_name;
  rel_stream << "\nDescription of the fitter:\n" << fitter_name;
  count_stream << "\nDescription of the fitter:\n" << fitter_name;
  res_stream << "\nDescription of the fitter:\n" << fitter_name;
  allres_stream << "\nDescription of the fitter:\n" << fitter_name;
  ressign_stream << "\nDescription of the fitter:\n" << fitter_name;

  char time_measure[256];
  boost::posix_time::time_duration time_used;
  time_used = boost::posix_time::second_clock::local_time() - start_time;
  sprintf(time_measure,"\nRan %i hours, %i minutes %i seconds.", 
               time_used.hours(), time_used.minutes(), 
               time_used.seconds());

  mean_stream << "\n" << time_measure << std::endl;
  error_stream << "\n" << time_measure << std::endl;
  rel_stream << "\n" << time_measure << std::endl;
  count_stream << "\n" << time_measure << std::endl;
  
  if (sim_count > 0) {
    res_stream << "\n\nLeast residuum size achieved was " << least_avg_res_size
               << " in step " << least_avg_res_index << "." << time_measure;
    allres_stream << "\n\nLeast residuum size achieved was " 
                 << least_avg_res_size << " in step " << least_avg_res_index 
                 << "." << time_measure;
    ressign_stream << "\n" << time_measure << std::endl;
  }

  mean_stream.close();
  error_stream.close();
  rel_stream.close();
  count_stream.close();
  res_stream.close();
  allres_stream.close();
  ressign_stream.close();
}

/**
 * Write the status of the fitting into a stream in XML format (this is
 * used for checkpointing).
 *
 * @param oxs The stream to be written into.
 */
alps::oxstream& AbstractFitter::write_to_xml(alps::oxstream& oxs) const
{
  oxs << alps::start_tag("PARAMETER_HISTORY");
  for (int i=0; i<param_history.size(); i++)
    oxs << alps::start_tag("HISTORY_ELEMENT") << param_history[i]
        << alps::end_tag("HISTORY_ELEMENT");
  oxs << alps::end_tag("PARAMETER_HISTORY");

  oxs << alps::start_tag("CURRENT_PARAMS");
  oxs << current_params;
  oxs << alps::end_tag("CURRENT_PARAMS");

  oxs << alps::start_tag("NEXT_PARAMETERS");
  std::list<FitterParamList>::const_iterator it = next_parameters.begin();
  for (; it != next_parameters.end(); it++)
    oxs << alps::start_tag("NEXT_PARAM") << *it << alps::end_tag("NEXT_PARAM");
  oxs << alps::end_tag("NEXT_PARAMETERS");

  oxs << alps::start_tag("ALL_EXP_RESULTS");
  oxs << all_exp_res;
  oxs << alps::end_tag("ALL_EXP_RESULTS");

  oxs << alps::start_tag("SIMULATION_RESULTS");
  for (int i=0; i<sim_results.size(); i++)
    oxs << alps::start_tag("SIM_RESULTS") << sim_results[i]
        << alps::end_tag("SIM_RESULTS");
  oxs << alps::end_tag("SIMULATION_RESULTS");

  oxs << alps::start_tag("RESIDUUM_HISTORY");
  for (int i=0; i<residuum_history.size(); i++)
    oxs << residuum_history[i];
  oxs << alps::end_tag("RESIDUUM_HISTORY");

  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","filenamebase")
      << alps::no_linebreak << filenamebase
      << alps::end_tag("FITTING_PARAMETER");
  oxs << alps::start_tag("FITTING_PARAMETER")
      << alps::attribute("name","current_filename")
      << alps::no_linebreak << current_filename
      << alps::end_tag("FITTING_PARAMETER");
  return oxs;
}

/**
 * Checks whether the fitting process is done.
 *
 * @param seconds_left The number of seconds left with respect to the time 
 *                     limit given to the fitter.
 */
bool AbstractFitter::is_done(double* seconds_left)
{
  using namespace boost::posix_time;
  if (end_time > start_time) {
    time_duration time_left = end_time - second_clock::local_time();
    *seconds_left = time_left.total_seconds();
    if (*seconds_left <= 0) {
      std::cout << "\nTime limit exceeded, ending\n";
      return true;
    }
  }
  std::cout << "Did " << step << " steps, maximal number of steps is "
            << max_steps << "\n";

  if (step >= max_steps) {
    std::cout << "maximal number of steps has been reached, ending ...\n";
    return true;
  }
  if (step < 0)
    return false;

  double last_res = residuum_history.back().avg_res_size;
  std::cout << "last residuum size is " << last_res << ", limit is " 
            << max_avg_residuum << "\n";

  if (last_res < max_avg_residuum) {
    std::cout << "\nGot a sufficiently small residuum, ending ...\n";
    return true;
  }
  return false;
}

/**
 * Find the variable parameters for the fitting in the set of parameters and set
 * the limitations for them.
 *
 * @param parms The fitting parameters
 */
std::vector<FittingParam> AbstractFitter::find_variable_parms
        (alps::Parameters parms) const
{
  Helper<std::string> myHelper(parms["VARIABLE_PARAMETERS"],parms);
  int size = myHelper.elemcount();
  std::vector<FittingParam> res;
  res.resize(size);
  for (int i=0; i<size; i++) {
    res[i].name = myHelper.getNextString();
    res[i].limit.validmin = false;
    res[i].limit.validmax = false;
  }
  
  // read the limits
  Helper<double> limitHelper;
  std::string tmp;
  int limits_size;
  if (parms.defined("LIMITS")) {
    limitHelper.setString(parms["LIMITS"]);
    limits_size = limitHelper.elemcount();
    if ((limits_size % 3) != 0) {
      std::cerr << parms["LIMITS"] << " is an invalid limit set.\n";
      boost::throw_exception(std::runtime_error("Invalid parameter input"));
    } else {
      for (int i=0; i<limits_size; i+=3) {
        tmp = limitHelper.getNextString();
        // find the position of the variable parameter
        int pos=0;
        while ((pos<size) && (res[pos].name != tmp)) pos++;
        if (pos==size) {
          std::cerr << "variable parameter " << tmp << " not found.\n"
                    << "Ignoring limits " << limitHelper.getNextElement()
                    << " and " << limitHelper.getNextElement() << "\n";
        } else {
          tmp = limitHelper.getNextString();
          if (tmp != "-") {
            res[pos].limit.min = alps::evaluate<double>(tmp,parms);
            res[pos].limit.validmin = true;
          }
          tmp = limitHelper.getNextString();
          if (tmp != "-") {
            res[pos].limit.max = alps::evaluate<double>(tmp,parms);
            res[pos].limit.validmax = true;
          }
          if (res[pos].limit.validmax && res[pos].limit.validmin
               && (res[pos].limit.min > res[pos].limit.max)) {
            std::cerr << "Invalid limit value for " << res[pos].name << "!\n"
                      << "minimum " << res[pos].limit.min
                      << " larger than maximum " << res[pos].limit.max << "\n";
            boost::throw_exception(std::runtime_error("Invalid limits"));
          }
        }
      }
    }
  }
  return res;
}

/**
 * Add a new parameter with the given value to all parameter sets in the 
 * parameter list. If a parameter with the given name already exists, its value
 * is overwritten.
 *
 * @param par_name The name of the parameter to be set
 * @param value    The value the new parameter is to be set to
 * @param list     The list of parameter sets to be updated.
 */
alps::ParameterList AbstractFitter::update_parameters(
               const std::string par_name, double value, 
               alps::ParameterList list) const
{
  alps::Parameter new_par(par_name, value);
  alps::ParameterList res;
  res = list;
  for (int j=0; j<list.size(); j++)
    res[j] << new_par;
  return res;
}

/** 
 * Adjust the temperature values for the Monte-Carlo simulations to the
 * temperature values of the chosen experimental results.
 * 
 * @param list The parameter list that is to be adjusted.
 */
alps::ParameterList AbstractFitter::adjust_parameters(
               const alps::ParameterList list,
               const alps::scheduler::ResultsType& exp_res)
{
  if (list.size() == 0) {
    std::cerr << "Initial parameter list has no elements, aborting.\n";
    boost::throw_exception(std::runtime_error("invalid parameter file"));
  }

  alps::ParameterList res;
  alps::Parameters current_params;
  current_params = list[0];
  
  unsigned int exp_size = exp_res.size();
  res.resize(exp_size);

  for (int i=0; i<exp_size; i++) {
    alps::Parameter param("T",exp_res[i].T);
    current_params << param;
    res[i] =  current_params;
  }
  return res;
}
