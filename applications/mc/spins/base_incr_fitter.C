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

#define CURVATURE_PREFACTOR 0.0000000001
#define INV_CURVATURE_PREFAC 1e10
#define INITIAL_ERROR_FACTOR 50.0

/**
 * The constructor, reads in the parameter file and internalizes the values.
 * 
 * @param filename The name of the file containing the parameter for the
 * fitting.
 */
BaseIncrFitter::BaseIncrFitter(const std::string& filename)
    : AbstractFitter(filename)
{
  exp_results.clear();
  alps::Parameters sim_param;
  {
    std::ifstream in(filename.c_str());
    in >> sim_param;
  }
  enlarge_trigger_factor = sim_param.value_or_default("ENLARGE_TRIGGER_FACTOR",
               25.0);
  steps_between_baseincr = sim_param.value_or_default("STEPS_BETWEEN_BASEINCR",
               10);

  if (sim_param.defined("INITIAL_BASIS_SIZE")) {
    int size = sim_param["INITIAL_BASIS_SIZE"];
    choose_exp_results(size,size,size);
  } else {
    int avg = sim_param.value_or_default("AVG_INITIAL_BASE_SIZE",
               num_variable_params+5);
    int min = sim_param.value_or_default("MIN_INITIAL_BASE_SIZE",avg-1);
    int max = sim_param.value_or_default("MAX_INITIAL_BASE_SIZE",avg+1);
    choose_exp_results(avg,min,max);
  }

}

/** 
 * Determine the absolute value of the curvature of the experimental data.
 */
void BaseIncrFitter::determine_curvatures()
{
  unsigned int exp_count = all_exp_res.size();
  curvature.resize(exp_count);
  curvature_sum = 0.0;
                                                                                
  double first_der_1, first_der_0, second_der, temp, curv;
  // the first and the last experiment get curvature 0.5 (default value)
                                                                                
  first_der_0 = (all_exp_res[1].mean - all_exp_res[0].mean)
               /(all_exp_res[1].T - all_exp_res[0].T);
  for (int i=1; i<exp_count-1; i++) {
    first_der_1  = (all_exp_res[i+1].mean - all_exp_res[i].mean)
                  /(all_exp_res[i+1].T - all_exp_res[i].T);
    second_der = (first_der_1 - first_der_0)
               / (all_exp_res[i].T - all_exp_res[i-1].T);
    second_der *= CURVATURE_PREFACTOR;
    temp = 1+first_der_1*first_der_0;
    curv = second_der / ( temp*sqrt(temp) );
    if (curv < 0.0)
      curvature[i] = -curv;
    else 
      curvature[i] = curv;
    first_der_0 = first_der_1;
  }
  for (int i=1; i<exp_count-1; i++)
    curvature_sum += curvature[i];

  // set the points at the end to the double value of the average
  double average = curvature_sum/exp_count;
  curvature[0] = 2*average;
  curvature[exp_count-1] = 2*average;
  curvature_sum += 4*average;
}

/**
 * Randomly choose a subset of the experimental data, with probability pro -
 * portional to the curvature.
 *
 * @param expected_number  The number of expected experimental results.
 * @param rem_curv_sum     The sum of the curvatures of all elements not in the 
 *                         fitting basis yet.
 */
void BaseIncrFitter::choose_exp_results(int expected_number,double /*rem_curv_sum*/)
{
  if (all_exp_res.empty()) {
    std::cerr << "no exponential results given, returning\n";
    return;
  }
  srand(time(NULL));
  exp_results.clear();
  if (curvature.empty())
    determine_curvatures();

  // every experimental result is chosen with probability
  // p_acc = scale/curvature
  int count = all_exp_res.size();

  // curvature_sum is the unscaled expected number of chosen
  // experiments
  double scale = (double)expected_number/curvature_sum;
  remaining_curv_sum = curvature_sum;
  for(unsigned int i=0; i<count; i++) {
    double r = (double)rand()/RAND_MAX;
    if (r < scale*curvature[i]) {
      // choose the i-th element
      exp_results.push_back(all_exp_res[i]);
      remaining_curv_sum -= curvature[i];
    }
  }
}

/**
 * Choose a number of experimental results. The (nbiggest) elements with the
 * largest curvature are chosen for sure, the remaining elements are chosen
 * randomly with probability proportional to their absolute curvature.
 *
 * @param expected The expected number of elements to be chosen.
 * @param min      The minimal number to be chosen.
 * @param max      The maximal number to be chosen.
 * @param nbiggest The number of data points to be chosen deterministically
 *                 according to the largest curvature.
 */
void BaseIncrFitter::save_choose(int expected, int min, int max, int nbiggest)
{
  if (expected < nbiggest) {
    std::cerr << "invalid input in save_choose : truncating the number of ranks to the number of expected exponential results for the basis\n";
    expected = nbiggest;
  }
  exp_results.clear();
  remaining_curv_sum = curvature_sum;
  // search the nbiggest elements with the biggest curvature
  int* ranking = new int[nbiggest];
  for (int i=0; i<nbiggest; i++)
    ranking[i]=-1;
  ranking[0]=0; // 0 is the best
  for (int pos=1; pos<all_exp_res.size(); pos++) {
    int i=0;
    while ((i<nbiggest) && (ranking[i] >=0) 
               && (curvature[pos] < curvature[ranking[i]])) i++;
    if (i<nbiggest) {
      for (int j=nbiggest-2; j>=i; j--) ranking[j+1] = ranking[j];
      ranking[i] = pos;
    }
  }
  for (int i=0; i<nbiggest; i++) {
    exp_results.push_back(all_exp_res[ranking[i]]);
    remaining_curv_sum -= curvature[i];
  }
  delete ranking;
  // choose the remaining elements by adding elements
  if (min > nbiggest)
    add_exp_results(expected-nbiggest, min-nbiggest, max-nbiggest);
}

/**
 * Choose experimental results. If possible, the number of added data points
 * lies between min and max
 * 
 * @param expected_number The expected number of newly added elements to the 
 *                        fitting base
 * @param min             The minimal number of base elements to be added
 * @param max             The maximal number of base elements to be added.
 */
void BaseIncrFitter::choose_exp_results(int expected_number, int min, int max)
{
  if (all_exp_res.empty()) {
    std::cerr << "no exponential results given, returning\n";
    return;
  }
  if (min > max) {
    std::cerr << "invalid input in choose_exp_results: min = " << min
              << ", max = " << max;
    return;
  }
  choose_exp_results(expected_number,remaining_curv_sum);
  int counter=-10;
  while (((exp_results.size()<min) || (exp_results.size()>max)) 
               && (counter <= 0)) {
    choose_exp_results(expected_number,remaining_curv_sum);
    counter++;
  }
  while((exp_results.size()<min) || (exp_results.size()>max)) {
    save_choose(expected_number, min, max, counter);
    counter++;
  }
}

/**
 * Decide whether the fitting basis should be enlarged.
 */
bool BaseIncrFitter::should_enlarge_basis()
{
  unsigned int size_now = exp_results.size();
  double last_res = residuum_history.back().avg_res_size;
  tracefile << "\nChecking whether the fitting basis should be enlarged:\n";

  if (size_now == all_exp_res.size()) {
    tracefile << " no, as all available data points are already used.\n";
    return false;
  }
  if (last_base_increment+steps_between_baseincr < step) {
    tracefile << " yes, no basis increment was done within the last "
              << steps_between_baseincr << " steps.\n";
    return true;
  }
  if (last_res < enlarge_trigger_factor*max_avg_residuum) {
    tracefile << " yes, got a small residuum close to the requested value:\n"
              << " last residuum size is " << last_res << ".\n Limit is "
              << enlarge_trigger_factor << "*" << max_avg_residuum << " = "
              << enlarge_trigger_factor*max_avg_residuum << "\n";
    return true;
  }
  tracefile << " no, there's no reason to do so.\n";
  return false;
}

/**
 * Enlarge the fitting base and adjust the parameter list accordingly.
 *
 * @param list The list of parameters to be adjusted to the newly chosen
 *             fitting base.
 */
alps::ParameterList BaseIncrFitter::enlarge_basis(
               const alps::ParameterList list) 
{
  unsigned int size_now = exp_results.size();
  alps::ParameterList res = list;
                                                                                
  add_exp_results(size_now, size_now-2, size_now+2);
  for (unsigned int i = size_now; i < exp_results.size(); i++) {
    res.push_back(res[0]);
    // set temperature and initial error limit
    alps::Parameter pE("ERROR_LIMIT",INITIAL_ERROR_FACTOR*exp_results[i].mean);
    alps::Parameter pT("T",exp_results[i].T);
    res[i] << pE;
    res[i] << pT;
  }
  tracefile << "\nThe basis for the fitting was enlarged to "
            << exp_results.size() << ".\n Previously, the size was "
            << size_now << "\n";
  last_base_increment = step;
  return res;
}

/**
 * Checks whether a result type is contained in the fitting base.
 *
 * @param element The result element to be looked for in the current fitting
 *                base.
 */
bool BaseIncrFitter::is_contained(const alps::scheduler::ResultType& element)
               const
{
  for (int i=0; i<exp_results.size(); i++)
    if (exp_results[i].T == element.T)
      return true;
  return false;
}

/**
 * Add a number of expperimental results to the fitting base. If possible,
 * the data points are chosen with probability proportional to the curvature.
 * If this strategy is not successful, a few data points with the biggest
 * curvature are chosen determinisitally and the rest of the points is chosen
 * randomly. After the end of the function, the fitting base has been increased
 * by a number of elements in the given range.
 *
 * @param expected The expected number of newly added data points.
 * @param min      The minimal number of newly added data points.
 * @param max      The maximal number of data points to add to the basis.
 */
void BaseIncrFitter::add_exp_results(int expected, int min, int max)
{
   if (min > max) {
    std::cerr << "invalid input in add_exp_results: min = " << min
              << ", max = " << max;
    return;
  }
  unsigned int all_count = all_exp_res.size();
  int remaining_elements = all_count - exp_results.size();
  if (remaining_elements == 0) {
    std::cout << "all exponential results are already chosen, returning\n";
    return;
  }
  if (expected >= remaining_elements) {
    // add all remaining elements to exp_results
    for (int i = 0; i<all_exp_res.size(); i++)
      if (!is_contained(all_exp_res[i]))
        exp_results.push_back(all_exp_res[i]);
    std::cout << "all exponential results are now chosen, returning\n";
    return;
  }
  // the normal case
  std::vector<unsigned int> tmp;
  int counter = -10;
  simple_add_exp_results(expected, remaining_curv_sum,tmp);
  while (((tmp.size() < min) || (tmp.size() > max)) && (counter<=0)) {
    tmp.clear();
    simple_add_exp_results(expected, remaining_curv_sum,tmp);
    counter++;
  }
  while ((tmp.size() < min) || (tmp.size() > max)) {
    tmp.clear();
    save_add_exp_results(expected,counter,tmp);
    counter++;
  }

  // add the chosen elements
  for (int i=0; i<tmp.size(); i++) {
    exp_results.push_back(all_exp_res[tmp[i]]);
    remaining_curv_sum -= curvature[tmp[i]];
  }
}

/**
 * Randomly choose a number of experimental results with probability
 * proportional to the curvature.
 *
 * @param expected The expected number of data points to be added to the base
 * @param rem_sum  The sum of the curvatures of all remaining data points
 * @param tmp      A vector to return the chosen data points.
 */
void BaseIncrFitter::simple_add_exp_results(int expected, double rem_sum,
               std::vector<unsigned int>& tmp) const
{
  tmp.clear();
  double scale = (double)expected/rem_sum;
  unsigned int all_count = all_exp_res.size();
  for (unsigned int i = 0; i<all_count; i++) {
    double r = (double)rand()/RAND_MAX;
    if ((r < scale*curvature[i]) && (!is_contained(all_exp_res[i])))
      // choose the i-th element
      tmp.push_back(i);
  }
}

/**
 * Randomly choose a number of experimental results with probability
 * proportional to the curvature.
 *
 * @param expected The expected number of data points to be added to the base
 * @param rem_sum  The sum of the curvatures of all remaining data points
 * @param positions A vector of already chosen data points.
 * @param tmp      The vector to return the chosen points.
 */
void BaseIncrFitter::simple_add_exp_results(int expected, double rem_sum,
               const std::vector<unsigned int>& positions,
               std::vector<unsigned int>& tmp) const
{
  tmp.clear();
  double scale = (double)expected/rem_sum;
  unsigned int all_count = all_exp_res.size();
  for (unsigned int i = 0; i<all_count; i++) {
    double r = (double)rand()/RAND_MAX;
    if ((r < scale*curvature[i]) && (!is_contained(all_exp_res[i]))
               && (!::is_contained(positions,i)))
      // choose the i-th element
      tmp.push_back(i);
  }
}

/**
 * Add a number of points with the largest curvature to the fitting base. 
 * Chosse the remaining number of points randomly with probability proportional
 * to the curvature.
 *
 * @param expected The expected number of points to add to the base.
 * @param nbiggest The number of base points with the largest curvature to be
 *                 taken for sure.
 * @param tmp      The vector to return the indices of the chosen data points.
 */
void BaseIncrFitter::save_add_exp_results(int expected, int nbiggest,
               std::vector<unsigned int>& tmp) const
{
  if (expected < nbiggest) {
    std::cerr << "invalid input in save_add_exp+results : \n"
              << "truncating the number of ranks to the number of expected "
              << "exponential results for the basis\n";
    expected = nbiggest;
  }
  tmp.clear();
  if (nbiggest <= 0) return;

  unsigned int all_count = all_exp_res.size();
  int remaining_elements = all_count - exp_results.size();
  if (remaining_elements == 0) {
    std::cout << "all exponential results are already chosen, returning\n";
  }
  if (expected >= remaining_elements) {
    // add all remaining elements to exp_results
    for (unsigned int i = 0; i<all_exp_res.size(); i++)
      if (!is_contained(all_exp_res[i]))
        tmp.push_back(i);
    std::cout << "all exponential results are now chosen, returning\n";
  }

  // choose the first nbiggest elements as the elements with the biggest
  // curvature
  double rcs = remaining_curv_sum;
  //double scale = (double)expected/remaining_curv_sum;
  int* ranking = new int[nbiggest];
  for (unsigned int i=0; i<nbiggest; i++)
    ranking[i]=-1;
  for (int pos=1; pos<all_exp_res.size(); pos++) {
    for (int i=0; i<nbiggest; i++)
    if (!is_contained(all_exp_res[pos])) {
      int i=0;
      while ((i<nbiggest) && (ranking[i] >=0)
                 && (curvature[pos] < curvature[ranking[i]])) i++;
      if (i<nbiggest) {
        for (int j=nbiggest-2; j>=i; j--) ranking[j+1] = ranking[j];
        ranking[i] = pos;
      }
    }
  }
  for (int i=0; i<nbiggest; i++) {
    tmp.push_back(ranking[i]);
    rcs -= curvature[ranking[i]];
  }
  std::vector<unsigned int> tmp2;
  if (expected>nbiggest) {
    simple_add_exp_results(expected-nbiggest,rcs,tmp,tmp2);
    for (unsigned int i=0; i<tmp2.size();i++)
      tmp.push_back(tmp2[i]);
  }
}

/**
 * Return a set to plot the ratio of the residual in percents of the mean
 * value of the experimental data.
 */
alps::plot::Set<double> BaseIncrFitter::get_res_percent_plot_set() const
{
  alps::plot::Set<double> res(alps::plot::xydy);
  res.clear();
  if (residuum_history.empty())
    return res;
                                                                                
  std::vector<ResidualType> res_vals = residuum_history.back().res_vals;
  for (unsigned int i=0; i<exp_results.size(); i++) {
    if (exp_results[i].mean != 0) {
      res << exp_results[i].T;
      res << 100.0*res_vals[i].mean/exp_results[i].mean;
      res << 100.0*res_vals[i].error/exp_results[i].mean;
    }
  }
  res << "Percentual size of Residuum";
  return res;
}

/**
 * Write the status information to an xml stream
 * 
 * @param oxs The xml stream
 */
alps::oxstream& BaseIncrFitter::write_to_xml(alps::oxstream& oxs) const
{
  AbstractFitter::write_to_xml(oxs);

  oxs << alps::start_tag("EXP_RESULTS");
  oxs << exp_results;
  oxs << alps::end_tag("EXP_RESULTS");
                                                                                
  oxs << alps::start_tag("CURVATURE");
  for (int i=0; i<curvature.size(); i++)
    oxs << alps::start_tag("CURVATURE") << alps::no_linebreak
        << curvature[i] << alps::end_tag();
  oxs << alps::end_tag("CURVATURE");
                                                                                
  oxs << alps::start_tag("CURVATURE_SUM") << alps::no_linebreak
      << curvature_sum << alps::end_tag("CURVATURE_SUM");
                                                                                
  oxs << alps::start_tag("REMAINING_CURV_SUM") << alps::no_linebreak
      << remaining_curv_sum << alps::end_tag("REMAINING_CURV_SUM");
  return oxs;
}

/**
 * Initialize the BaseIncrFitter reading parameters from the in stream in 
 * XML format.
 * 
 * @param is The infile stream containing the parameters in XML format.
 */
void BaseIncrFitter::initialize(std::ifstream& is)
{
  AbstractFitter::initialize(is);
  alps::XMLTag tag = alps::parse_tag(is, true);
  
  exp_results.clear();
  if (tag.type == alps::XMLTag::SINGLE) {
    std::cerr << "Invalid checkpoint file: No experimental results.\n";
    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  } else {
    is >> exp_results;
    std::cout << "read full set of experimental results.\n";
    tag = alps::parse_tag(is, true);
  }
                                                                                
  tag = alps::parse_tag(is, true);
  curvature.clear();
  if (tag.type == alps::XMLTag::SINGLE) {
    std::cout << "Corrupted checkpoint file: Curvature not stored.\n";
    std::cout << "Recovering ...";
    determine_curvatures();
    std::cout << " ok, continuing.\n";
  } else {
    tag = alps::parse_tag(is,true);
    while (tag.name != "/CURVATURE") {
      double component = alps::evaluate<double>(alps::parse_content(is));
      curvature.push_back(component);
      tag = alps::parse_tag(is,true);
      tag = alps::parse_tag(is,true);
    }
  }
  tag = alps::parse_tag(is,true);
                                                                                
  if (tag.type == alps::XMLTag::SINGLE) {
    std::cerr << "Curvatures sum not stored, using recomputed sum\n";
  } else {
    curvature_sum = alps::evaluate<double>(alps::parse_content(is));
    tag = alps::parse_tag(is,true);
  }
                                                                                
  tag = alps::parse_tag(is,true);
  if (tag.type == alps::XMLTag::SINGLE) {
    std::cerr<<"Invalid checkpoint file: Remaining curvature sum is missing.\n";    boost::throw_exception(std::runtime_error("Invalid checkpoint file"));
  } else {
    remaining_curv_sum = alps::evaluate<double>(alps::parse_content(is));
    tag = alps::parse_tag(is,true);
  }
}
