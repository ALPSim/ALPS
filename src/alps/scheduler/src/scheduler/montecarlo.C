/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/scheduler/montecarlo.h>
#include <alps/scheduler/types.h>
#include <alps/scheduler/scheduler.h>
#include <alps/alea/observables_io.hpp>
#include <alps/parser/xslt_path.h>
#include <alps/osiris/mpdump.h>
#include <alps/osiris/std/vector.h>
#include <filesystem>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>

namespace alps {
namespace scheduler {

namespace {

// Where a run keeps what it measured, in its HDF5 file: the results for
// reading, and the measurements to continue from.
const std::string run_results = "/simulation/realizations/0/clones/0/results";
const std::string run_measurements = "/simulation/realizations/0/clones/0/measurements";

// Digits that show the mean to about a thousandth of its error
int mean_digits(double mean, double error, int otherwise)
{
  const double digits = 4. - std::log10(std::abs(error / mean));
  return digits >= 3. && digits < 20. ? static_cast<int>(digits) : otherwise;
}

bool error_underflow(double mean, double error)
{
  return error != 0. && mean != 0.
    && std::abs(mean) * 10. * std::sqrt(std::numeric_limits<double>::epsilon()) > std::abs(error);
}

const char* convergence_text(alps::alea::error_convergence c)
{
  return c == alps::alea::CONVERGED ? "yes" : c == alps::alea::MAYBE_CONVERGED ? "maybe" : "no";
}

struct average_values
{
  explicit average_values(const alps::alea::observable_result& r)
    : count(r.count()), mean(r.mean_vector()), error(r.error_vector()), variance(r.variance_vector()),
      tau(r.tau_vector()), converged(r.converged_vector()), binning(r.size())
  {
    for (std::size_t i = 0; i != r.size(); ++i)
      binning[i] = r.binning_errors(i);
    // A binning analysis belongs to a measured observable. Other results with
    // bins are functions of observables, evaluated by jackknife.
    method = !binning.empty() && !binning[0].empty() ? "binning" : r.has_bins() ? "jackknife" : "simple";
  }

  std::uint64_t count;
  std::vector<double> mean, error, variance, tau;
  std::vector<alps::alea::error_convergence> converged;
  std::vector<std::vector<double> > binning;
  std::string method;
};

// The elements inside one <SCALAR_AVERAGE>, for component i
void write_component(oxstream& out, const average_values& v, std::size_t i, int digits)
{
  const bool jackknife = v.method == "jackknife";
  out << start_tag("COUNT") << no_linebreak << v.count << end_tag("COUNT");
  out << start_tag("MEAN") << no_linebreak << attribute("method", jackknife ? "jackknife" : "simple")
      << precision(v.mean[i], digits) << end_tag("MEAN");
  out << start_tag("ERROR") << attribute("converged", convergence_text(v.converged[i]));
  if (error_underflow(v.mean[i], v.error[i]))
    out << attribute("underflow", "true");
  out << attribute("method", v.method) << no_linebreak << precision(v.error[i], 3) << end_tag("ERROR");
  if (v.variance[i] == v.variance[i])
    out << start_tag("VARIANCE") << no_linebreak << attribute("method", "simple")
        << precision(v.variance[i], 3) << end_tag("VARIANCE");
  if (v.tau[i] == v.tau[i])
    out << start_tag("AUTOCORR") << no_linebreak << attribute("method", "binning")
        << precision(v.tau[i], 3) << end_tag("AUTOCORR");
  for (std::size_t level = 0; level != v.binning[i].size(); ++level)
    out << start_tag("BINNED") << attribute("size", boost::lexical_cast<std::string>(std::uint64_t(1) << level))
        << no_linebreak << start_tag("COUNT") << (v.count >> level) << end_tag("COUNT")
        << start_tag("ERROR") << attribute("method", "simple") << no_linebreak
        << precision(v.binning[i][level], 3) << end_tag("ERROR") << end_tag("BINNED");
}

void write_xml_average(oxstream& out, const std::string& name, const alps::alea::observable_result& r)
{
  if (r.count() == 0)
    return;
  if (r.is_histogram()) {
    const std::vector<std::uint64_t>& counts = r.histogram();
    out << start_tag("HISTOGRAM") << attribute("name", name) << attribute("nvalues", counts.size());
    for (std::size_t i = 0; i != counts.size(); ++i)
      out << start_tag("ENTRY") << attribute("indexvalue", i)
          << start_tag("COUNT") << no_linebreak << r.count() << end_tag("COUNT")
          << start_tag("VALUE") << no_linebreak << counts[i] << end_tag("VALUE")
          << end_tag("ENTRY");
    out << end_tag("HISTOGRAM");
    return;
  }
  const average_values v(r);
  if (r.is_scalar()) {
    out << start_tag("SCALAR_AVERAGE") << attribute("name", name);
    if (!r.sign().empty())
      out << attribute("signed", "true") << attribute("sign", r.sign());
    write_component(out, v, 0, mean_digits(v.mean[0], v.error[0], 8));
    out << end_tag("SCALAR_AVERAGE");
    return;
  }
  out << start_tag("VECTOR_AVERAGE") << attribute("name", name) << attribute("nvalues", r.size());
  if (!r.sign().empty())
    out << attribute("signed", "true") << attribute("sign", r.sign());
  for (std::size_t i = 0; i != r.size(); ++i) {
    const std::string label = i < r.labels().size() && !r.labels()[i].empty()
      ? r.labels()[i] : boost::lexical_cast<std::string>(i);
    out << start_tag("SCALAR_AVERAGE") << attribute("indexvalue", label);
    write_component(out, v, i, v.count == 1 ? 19 : mean_digits(v.mean[i], v.error[i], 8));
    out << end_tag("SCALAR_AVERAGE");
  }
  out << end_tag("VECTOR_AVERAGE");
}

} // namespace

void write_xml_averages(oxstream& out, const alps::alea::result_set& results, int id)
{
  out << start_tag("AVERAGES");
  if (id >= 0)
    out << attribute("id", id);
  for (const std::string& name : results.names())
    write_xml_average(out, name, results[name]);
  out << end_tag("AVERAGES");
}

namespace {

double xml_number(const std::string& text, const std::string& name)
{
  const char* begin = text.c_str();
  char* end = nullptr;
  const double value = std::strtod(begin, &end);   // also nan and inf
  if (end == begin || *end != 0)
    boost::throw_exception(std::runtime_error("Cannot read the number '" + text + "' in <AVERAGES> of '" + name + "'"));
  return value;
}

alps::alea::error_convergence xml_convergence(const XMLAttributes& attributes)
{
  const std::string converged = attributes.value_or_default("converged", "maybe");
  return converged == "yes" ? alps::alea::CONVERGED
    : converged == "no" ? alps::alea::NOT_CONVERGED : alps::alea::MAYBE_CONVERGED;
}

} // namespace

AveragesXMLHandler::AveragesXMLHandler(std::vector<alps::alea::result_set>& sets)
  : XMLHandlerBase("AVERAGES"), sets_(sets), histogram_(false), histogram_min_(0)
{
}

void AveragesXMLHandler::start_element(const std::string& name, const XMLAttributes& attributes,
                                       xml::tag_type type)
{
  if (type != xml::element)
    return;
  text_.clear();
  if (path_.empty() && name == "AVERAGES") {
    sets_.push_back(alps::alea::result_set());
    path_.push_back(name);
    return;
  }
  const std::string parent = path_.empty() ? std::string() : path_.back();
  path_.push_back(name);
  if (parent == "AVERAGES") {
    name_ = attributes.value_or_default("name", "");
    numbers_ = alps::alea::stored_numbers();
    histogram_ = name == "HISTOGRAM";
    counts_.clear();
    histogram_min_ = 0;
    if (name == "SCALAR_AVERAGE" || name == "VECTOR_AVERAGE") {
      numbers_.scalar = name == "SCALAR_AVERAGE";
      if (attributes.value_or_default("signed", "") == "true")
        numbers_.sign = attributes.value_or_default("sign", "Sign");
      if (numbers_.scalar)
        finish_component();
    }
  } else if (parent == "VECTOR_AVERAGE" && name == "SCALAR_AVERAGE") {
    finish_component();
    numbers_.labels.push_back(attributes.value_or_default("indexvalue", ""));
  } else if (parent == "SCALAR_AVERAGE" && name == "ERROR") {
    numbers_.converged.back() = xml_convergence(attributes);
  }
}

// a new component, with what is unknown until it is read
void AveragesXMLHandler::finish_component()
{
  const double nan = std::numeric_limits<double>::quiet_NaN();
  numbers_.mean.push_back(nan);
  numbers_.error.push_back(nan);
  numbers_.variance.push_back(nan);
  numbers_.tau.push_back(nan);
  numbers_.converged.push_back(alps::alea::MAYBE_CONVERGED);
  numbers_.binning_errors.push_back(std::vector<double>());
}

void AveragesXMLHandler::text(const std::string& text)
{
  text_ += text;
}

void AveragesXMLHandler::end_element(const std::string& name, xml::tag_type type)
{
  if (type != xml::element || path_.empty())
    return;
  path_.pop_back();
  const std::string parent = path_.empty() ? std::string() : path_.back();
  std::string content = text_;
  text_.clear();
  content.erase(0, content.find_first_not_of(" \t\r\n"));
  content.erase(content.find_last_not_of(" \t\r\n") + 1);

  if (parent == "SCALAR_AVERAGE" && !content.empty()) {
    const double value = xml_number(content, name_);
    if (name == "COUNT") {
      if (numbers_.count == 0)
        numbers_.count = static_cast<std::uint64_t>(value);
    } else if (name == "MEAN") {
      numbers_.mean.back() = value;
    } else if (name == "ERROR") {
      numbers_.error.back() = value;
    } else if (name == "VARIANCE") {
      numbers_.variance.back() = value;
    } else if (name == "AUTOCORR") {
      numbers_.tau.back() = value;
    }
  } else if (parent == "BINNED" && name == "ERROR" && !content.empty()) {
    numbers_.binning_errors.back().push_back(xml_number(content, name_));
  } else if (parent == "ENTRY" && name == "VALUE" && !content.empty()) {
    counts_.push_back(static_cast<std::uint64_t>(std::llround(xml_number(content, name_))));
  } else if (parent == "AVERAGES") {
    if (name_.empty())
      boost::throw_exception(std::runtime_error("An element <" + name + "> in <AVERAGES> has no name"));
    if (histogram_) {
      if (!counts_.empty())
        sets_.back().insert(name_, alps::alea::stored_histogram(name_, histogram_min_, counts_));
    } else if (name == "SCALAR_AVERAGE" || name == "VECTOR_AVERAGE") {
      // Without labels, the writer gives the index of each component.
      bool labelled = false;
      for (std::size_t i = 0; i != numbers_.labels.size(); ++i)
        labelled = labelled || numbers_.labels[i] != boost::lexical_cast<std::string>(i);
      if (!labelled)
        numbers_.labels.clear();
      // Variances, autocorrelation times and binning analyses that are
      // missing for every component are unknown.
      auto all_nan = [](const std::vector<double>& v) {
        for (double x : v) if (x == x) return false;
        return true;
      };
      if (all_nan(numbers_.variance)) numbers_.variance.clear();
      if (all_nan(numbers_.tau)) numbers_.tau.clear();
      bool binned = false;
      for (const std::vector<double>& levels : numbers_.binning_errors)
        binned = binned || !levels.empty();
      if (!binned)
        numbers_.binning_errors.clear();
      sets_.back().insert(name_, alps::alea::stored_result(name_, numbers_));
    }
  }
}

alps::alea::result_set read_xml_averages(std::istream& in, const XMLTag& tag)
{
  std::vector<alps::alea::result_set> sets;
  AveragesXMLHandler handler(sets);
  handler.start_element(tag.name, tag.attributes, xml::element);
  if (tag.type != XMLTag::SINGLE) {
    for (int depth = 1; depth > 0;) {
      XMLTag next = parse_tag(in);
      if (next.type == XMLTag::OPENING) {
        handler.start_element(next.name, next.attributes, xml::element);
        ++depth;
        const std::string content = parse_content(in);
        if (!content.empty())
          handler.text(content);
      } else if (next.type == XMLTag::SINGLE) {
        handler.start_element(next.name, next.attributes, xml::element);
        handler.end_element(next.name, xml::element);
      } else if (next.type == XMLTag::CLOSING) {
        handler.end_element(next.name.substr(1), xml::element);
        --depth;
      }
    }
  }
  handler.end_element(tag.name, xml::element);
  return sets.front();
}

/**
 * Returns the summary of the MCSimulation for the Observable given by the
 * name
 *
 * @params name the name of the observable
 */
ResultType MCSimulation::get_summary(const std::string name) const
{
  ResultType res;
  const alps::alea::result_set results = get_measurements();
  res.T = parms["T"];
  res.mean = results[name].mean();
  res.error = results[name].error();
  res.count = results[name].count();
  return res;
}

/**
 * Returns the summary, the name of the Observable is specified in the job file
 */
ResultType MCSimulation::get_summary() const
{
  std::string theName;
    if (parms.defined("SUMMARY_VARIABLE"))
      theName = parms["SUMMARY_VARIABLE"];
    else
      theName = parms["ERROR_VARIABLE"];
  std::cout << "\nMaking summary for the observable " << theName << "\n";
  if (theName.length() == 0) {
    std::cout << "cannot find the tag ERROR_VARIABLE in the parameter set\n"
              << "so summary can be made\n";
    boost::throw_exception(std::runtime_error("no variable name to make summary after"));
  }
  return get_summary(theName);
}

// collect all measurements
alps::alea::result_set MCSimulation::get_measurements() const
{
  // Runs that measured with this version are merged. A run read from a file
  // of an earlier version has results only, which cannot be merged.
  alps::alea::observable_set merged;
  std::vector<alps::alea::result_set> earlier;
  std::size_t number_merged = 0;
  ProcessList where_master;
  for (unsigned int i=0;i<runs.size();i++) {
    if(workerstatus[i]==RemoteRun) {
      if(!runs[i])
        boost::throw_exception(std::runtime_error( "run does not exist in MCSimulation::get_measurements"));
      where_master.push_back( Process(dynamic_cast<const RemoteWorker*>(runs[i])->process()));
    }
    else if(runs[i]) {
      const MCRun& run = dynamic_cast<const MCRun&>(*runs[i]);
      if (run.can_continue()) {
        merged.merge(run.get_results());
        ++number_merged;
      } else
        earlier.push_back(run.get_results());
    }
  }
  if(!where_master.empty()) {
    // broadcast request to all slaves
    OMPDump send;
    send.send(where_master,MCMP_get_measurements);
    // collect results
    for (unsigned int i=0;i<where_master.size();i++) {
      // receive dump from remote process, abort if error
      IMPDump receive(where_master[i],MCMP_measurements);
      // bytes, not a string: a string on a dump ends at the first zero
      std::vector<char> packed;
      receive >> packed;
      merged.merge(alps::alea::unpack(std::string(packed.begin(), packed.end())));
      ++number_merged;
    }
  }

  alps::alea::result_set results;
  if (earlier.empty())
    results = merged.results();
  else if (earlier.size() == 1 && number_merged == 0)
    results = earlier.front();
  else if (!stored_results().names().empty())
    results = stored_results();
  else
    boost::throw_exception(std::runtime_error(
      "The runs of this simulation were written by an earlier version of ALPS and cannot be merged"
      " by this one, and the task file holds no results. Evaluate the simulation with the version"
      " that wrote it."));

  // results of an earlier evaluation, kept with the task
  const alps::alea::result_set& stored = stored_results();
  for (const std::string& name : stored.names())
    if (!results.has(name))
      results.insert(name, stored[name]);
  for (const std::string& name : added_results_.names())
    results.insert(name, added_results_[name]);
  return results;
}


void MCSimulation::add_result(const std::string& name, const alps::alea::observable_result& result)
{
  added_results_.insert(name, result);
}


std::string MCSimulation::worker_tag() const 
{
  return "MCRUN";
}


void MCSimulation::write_xml_body(oxstream& out, std::filesystem::path const& fn, bool writeall) const {
  if (writeall)
    write_xml_averages(out, get_measurements());
  WorkerTask::write_xml_body(out,fn,writeall);
}


void MCRun::print_copyright(std::ostream& out) 
{
  out << "Non-copyrighted Monte Carlo program. Please insert your own copyright statement by overwriting the print_copyright static member function of your MCRun class or the print_copyright virtual function of your factory class.\n\n";
}


MCRun::MCRun(const ProcessList& w,const Parameters&  myparms,int n)
  : Worker(w,myparms,n)
{
}


void MCRun::load_worker(IDump& dump)
{
  Worker::load_worker(dump);
  // Before ALPS 2.1 the measurements were part of this file.
  if(node==0 && dump.version() < 400)
    boost::throw_exception(std::runtime_error(
      "This run was written by a version of ALPS before 2.1, which kept the measurements in the"
      " run file itself. This version reads the measurements of a run from its HDF5 file only."));
  load(dump);
}

void MCRun::save_worker(ODump& dump) const
{
  Worker::save_worker(dump);
  save(dump);
}

alps::alea::result_set MCRun::get_results() const
{
  return earlier_results_ ? *earlier_results_ : measurements.results();
}

void MCRun::save(hdf5::archive & ar) const {
  Worker::save(ar);
  if(node==0) {
    alps::alea::save_results(ar, run_results, get_results());
    if (can_continue())
      alps::alea::save_checkpoint(ar, run_measurements, measurements);
  }
}

void MCRun::load(hdf5::archive & ar) {
  Worker::load(ar);
  if(node!=0)
    return;
  if (ar.is_group(run_measurements))
    alps::alea::load_checkpoint(ar, run_measurements, measurements);
  else if (ar.is_group(run_results)) {
    // A file of an earlier version of ALPS: what it measured can be read as
    // results. The run cannot go on, see run().
    const alps::alea::result_set results = alps::alea::load_results(ar, run_results);
    if (!results.names().empty())
      earlier_results_.reset(new alps::alea::result_set(results));
  }
}

void MCRun::save(ODump&) const
{
}

void MCRun::load(IDump&)
{
}

// start/restart the run
std::string MCRun::work_phase()
{
  return is_thermalized() ? "running" : "equilibrating";
}
                
// do some work, either thermalization or measurements serrps
void MCRun::run()
{
  if (!can_continue())
    boost::throw_exception(std::runtime_error(
      "This run was written by an earlier version of ALPS and cannot be continued by this one."
      " Its results can still be read. To measure more, finish the simulation with the version"
      " that wrote it, or start a new simulation."));
  bool thermalized=is_thermalized();
  // do some work
  Worker::run();
  // check if it just became thermalized
  if(!thermalized&&node==0&&is_thermalized()) {
    measurements.reset();
    change_phase("running");
  }
}

// do one sweep
bool MCRun::is_thermalized() const
{
  // this should be implemented
  boost::throw_exception( std::logic_error("is_thermalized needs to be implemented"));
  return false;
}


void MCRun::write_xml(const std::filesystem::path& name) const
{
  oxstream xml(name.parent_path()/(name.filename().string()+".xml"));
  std::filesystem::path fn_hdf5(name.parent_path()/(name.filename().string()+".hdf"));
  const alps::alea::result_set results = get_results();

  xml << header("UTF-8") << stylesheet(xslt_path("ALPS.xsl"));
  xml << start_tag("SIMULATION") << xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
        << attribute("xsi:noNamespaceSchemaLocation","http://xml.comp-phys.org/2002/10/ALPS.xsd");
  xml << parms;
  write_xml_averages(xml, results);
  xml << start_tag("MCRUN");
  xml << start_tag("CHECKPOINT") << attribute("format","osiris")
      << attribute("file", name.string())
      << end_tag("CHECKPOINT");
  xml << start_tag("CHECKPOINT") << attribute("format","hdf5")
      << attribute("file", fn_hdf5.string())
      << end_tag("CHECKPOINT");

  xml << get_info();
  write_xml_averages(xml, results);
  xml << end_tag("MCRUN") << end_tag("SIMULATION");
}


bool MCRun::handle_message(const Process& runmaster,int32_t tag)
{
  IMPDump message;
  OMPDump dump;
  switch(tag) {
  case MCMP_get_measurements:
    message.receive(runmaster,MCMP_get_measurements);
    {
      const std::string packed = alps::alea::pack(get_results());
      dump << std::vector<char>(packed.begin(), packed.end());
    }
    dump.send(runmaster,MCMP_measurements);
    return true;
  default:
    return Worker::handle_message(runmaster,tag);
  }
}
  

void DummyMCRun::dostep()
{
  boost::throw_exception(std::logic_error("User-level checkpointing needs to be implemented for restarting from a checkpoint\n"));
}

double DummyMCRun::work_done() const 
{
  boost::throw_exception(std::logic_error("User-level checkpointing needs to be implemented for restarting from a checkpoint\n"));
  return 0.;
}

// astreich, 06/23
ResultType DummyMCRun::get_summary() const
{
  boost::throw_exception(std::logic_error("User-level checkpointing needs to be implemented for restarting from a checkpoint\n"));
  ResultType res;
  res.count = 0;
  return res;
}

void MCSimulation::load(hdf5::archive & ar) {
  WorkerTask::load(ar);
  if (ar.is_group("/simulation/results"))
    stored_results_ = alps::alea::load_results(ar, "/simulation/results");
}

void MCSimulation::save(hdf5::archive & ar) const {
  WorkerTask::save(ar);
  alps::alea::save_results(ar, "/simulation/results", get_measurements());
}

// The results in the XML file of the task are used where its HDF5 file has none.
void MCSimulation::handle_tag(std::istream& infile, const XMLTag& tag) 
{
  if (tag.name!="AVERAGES")
    WorkerTask::handle_tag(infile,tag);
  else
    xml_results_ = read_xml_averages(infile, tag);
}

const alps::alea::result_set& MCSimulation::stored_results() const
{
  return stored_results_.names().empty() ? xml_results_ : stored_results_;
}

DummyMCRun::DummyMCRun()
  : MCRun(ProcessList(),Parameters(),0) 
{
}

DummyMCRun::DummyMCRun(const ProcessList& w,const Parameters& p,int n)
: MCRun(w,p,n) 
{
}

} // namespace scheduler
} // namespace alps
