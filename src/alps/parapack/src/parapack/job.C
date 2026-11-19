/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2013 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/
#include <alps/parapack/job.h>
#include <alps/utilities/fs/absolute.hpp>
#include <alps/parapack/clone.h>
#include <alps/parapack/filelock.h>
#include <alps/parapack/measurement.h>
#include <alps/parapack/simulation_p.h>

#include <alps/alea/observables_io.hpp>
#include <alps/hdf5.hpp>

// some file (probably a python header) defines a tolower macro ...
#undef tolower
#undef toupper

#include <boost/config.hpp>
#include <boost/foreach.hpp>
#include <regex>

namespace alps {

namespace {

// Where the results of a task are written: one file, or one per replica
std::filesystem::path results_file(std::filesystem::path const& basedir, std::string const& base,
  std::size_t replicas, std::size_t i) {
  return alps::fs::absolute(std::filesystem::path(replicas == 1 ? base + ".out.h5" :
    base + ".replica" + boost::lexical_cast<std::string>(i+1) + ".h5"), basedir);
}

bool has_results(std::filesystem::path const& file) {
  if (!exists(file)) return false;
  bool found;
  #pragma omp critical (hdf5io)
  {
    hdf5::archive h5(file.string());
    found = h5.is_group("/simulation/results");
  }
  return found;
}

// The results written by the last evaluation of the task, if any
std::vector<alps::alea::result_set> load_results(std::filesystem::path const& basedir,
  std::string const& base) {
  std::filesystem::path single = results_file(basedir, base, 1, 0);
  std::filesystem::path first = results_file(basedir, base, 2, 0);
  bool replicas = has_results(first) &&
    (!has_results(single) || last_write_time(first) > last_write_time(single));
  std::vector<alps::alea::result_set> results;
  for (std::size_t i = 0; ; ++i) {
    std::filesystem::path file = replicas ? results_file(basedir, base, 2, i) : single;
    if (!has_results(file)) break;
    #pragma omp critical (hdf5io)
    {
      hdf5::archive h5(file.string());
      results.push_back(alps::alea::load_results(h5, "/simulation/results"));
    }
    if (!replicas) break;
  }
  return results;
}

} // end namespace

//
// task
//

task::task() : status_(task_status::Undefined) {}

task::task(std::filesystem::path const& file) : status_(task_status::Undefined) {
  basedir_ = file.parent_path();
  file_in_str_ = file.filename().string();
  file_out_str_ = std::regex_replace(file.filename().string(), std::regex("\\.in\\.xml$"), ".out.xml");
  if (file_in_str_ == file_out_str_) {
    file_in_str_ = std::regex_replace(file.filename().string(), std::regex("\\.out\\.xml$"), ".in.xml");
    file_out_str_ = file.filename().string();
  }
  base_ = std::regex_replace(file_out_str_, std::regex("\\.out\\.xml$"), "");
}

bool task::on_memory() const {
  return status_ == task_status::Ready || status_ == task_status::Running ||
    status_ == task_status::Continuing || status_ == task_status::Idling;
}

task::range_type const& task::num_clones() const {
  if (!on_memory()) boost::throw_exception(std::logic_error(
    "task::num_clones() task not loaded"));
  return num_clones_;
}

uint32_t task::num_running() const {
  if (!on_memory()) boost::throw_exception(std::logic_error(
    "task::num_running() task not loaded"));
  return running_.size();
}

uint32_t task::num_suspended() const {
  if (!on_memory()) boost::throw_exception(std::logic_error(
    "task::num_suspended() task not loaded"));
  return suspended_.size();
}

uint32_t task::num_finished() const {
  if (!on_memory()) boost::throw_exception(std::logic_error(
    "task::num_finished() task not loaded"));
  return finished_.size();
}

uint32_t task::num_started() const {
  return num_running() + num_suspended() + num_finished();
}

void task::load() {
  if (on_memory()) boost::throw_exception(std::logic_error("task::load() task already loaded"));
  params_.clear();
  obs_.clear();
  clone_info_.clear();
  std::filesystem::path file_in = alps::fs::absolute(std::filesystem::path(file_in_str_), basedir_);
  std::filesystem::path file_out = alps::fs::absolute(std::filesystem::path(file_out_str_), basedir_);
  simulation_xml_handler handler(params_, clone_info_);
  XMLParser parser(handler);
  if (!exists(file_out)) {
    parser.parse(file_in);
  } else {
    filelock lock(file_out, /* lock_now = */ true, /* wait = */ 60);
    if (!lock.locked())
      boost::throw_exception(std::runtime_error("task::load() lock failed after 60 seconds"));
    parser.parse(file_out);
    obs_ = load_results(basedir_, base_);
  }

  num_clones_ = range_type(params_.value_or_default("NUM_CLONES", "1"), params_);
  clone_status_.clear();
  clone_master_.clear();
  running_.clear();
  suspended_.clear();
  finished_.clear();
  BOOST_FOREACH(clone_info const& info, clone_info_) {
    if (info.progress() < 1){
      clone_status_.push_back(clone_status::Suspended);
      suspended_.insert(info.clone_id());
    } else {
      clone_status_.push_back(clone_status::Finished);
      finished_.insert(info.clone_id());
    }
    clone_master_.push_back(Process());
  }

  status_ = task_status::Ready;
  progress_ = calc_progress();
  boost::tie(weight_, dump_weight_) = calc_weight();
  status_ = calc_status();
}

void task::save(alps::parapack::option const& opt) const {
  if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
  std::filesystem::path file_out = alps::fs::absolute(std::filesystem::path(file_out_str_), basedir_);
  {
    filelock lock(file_out, /* lock_now = */ true, /* wait = */ 60);
    if (!lock.locked())
      boost::throw_exception(std::runtime_error("task::load() lock failed after 60 seconds"));
    if (exists(file_out)) {
      // keep the results that were written with the task
      simulation_xml_writer(file_out, opt.write_xml, true, params_, load_results(basedir_, base_),
        clone_info_);
    } else {
      simulation_xml_writer(file_out, opt.write_xml, true, params_, obs_, clone_info_);
    }
  }
}

void task::save_observable(alps::parapack::option const& opt) const {
  if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
  std::filesystem::path file_out = alps::fs::absolute(std::filesystem::path(file_out_str_), basedir_);
  {
    filelock lock(file_out, /* lock_now = */ true, /* wait = */ 60);
    if (!lock.locked())
      boost::throw_exception(std::runtime_error("task::load() lock failed after 60 seconds"));
    if (exists(file_out)) {
      Parameters params_tmp;
      std::deque<clone_info> clone_info_tmp;
      simulation_xml_handler handler(params_tmp, clone_info_tmp);
      XMLParser parser(handler);
      parser.parse(file_out);
      if (obs_.size() > 1 && !params_tmp.defined("NUM_REPLICAS"))
        params_tmp["NUM_REPLICAS"] = obs_.size();
      simulation_xml_writer(file_out, opt.write_xml, true, params_tmp, obs_, clone_info_tmp);
      for (std::size_t i = 0; i < obs_.size(); ++i) {
        Parameters p = params_tmp;
        if (obs_.size() > 1 && !p.defined("REPLICA")) p["REPLICA"] = i+1;
        #pragma omp critical (hdf5io)
        {
          hdf5::archive h5(results_file(basedir_, base_, obs_.size(), i).string(), "a");
          h5["/parameters"] << p;
          alps::alea::save_results(h5, "/simulation/results", obs_[i]);
        }
      }
    } else {
      boost::throw_exception(std::logic_error("task::save_observable()"));
    }
  }
}

void task::halt() {
  if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
  if (running_.size()) boost::throw_exception(std::logic_error("running clone exists"));

  switch (status_) {
  case task_status::Ready:
    status_ = task_status::NotStarted;
    break;
  case task_status::Running:
    status_ = task_status::Suspended;
    break;
  case task_status::Continuing:
    status_ = task_status::Finished;
    break;
  case task_status::Idling:
    status_ = task_status::Completed;
    break;
  default:
    boost::throw_exception(std::logic_error("unknown task_status"));
  }

  params_.clear();
  obs_.clear();
  clone_status_.clear();
  clone_master_.clear();
  clone_info_.clear();
  running_.clear();
  suspended_.clear();
  finished_.clear();
}

void task::check_parameter(alps::parapack::option const& opt) {
  // change in NUM_CLONES : keep old calculations and add new ones
  // change in SEED : ignored
  // change in other parameters : throw away all the old clones

  std::filesystem::path file_in = alps::fs::absolute(std::filesystem::path(file_in_str_), basedir_);
  std::filesystem::path file_out = alps::fs::absolute(std::filesystem::path(file_out_str_), basedir_);
  if (exists(file_out)) {
    // loading parameters from *.in.xml
    Parameters params_in;
    simulation_parameters_xml_handler handler(params_in);
    XMLParser parser(handler);
    parser.parse(file_in);

    // loading from *.out.xml
    load();
    if (params_.defined("SEED")) params_in["SEED"] = params_["SEED"];

    bool changed = false;
    BOOST_FOREACH(Parameter const& p, params_) {
      if (p.key() != "NUM_CLONES")
        if (!params_in.defined(p.key()) || params_in[p.key()] != p.value()) changed = true;
    }
    BOOST_FOREACH(Parameter const& p, params_in) {
      if (p.key() != "NUM_CLONES")
        if (!params_.defined(p.key())) changed = true;
    }
    if (changed) {
      std::cout << "Info: parameters in " << logger::task(task_id_) << " have been changed. "
                << "All the clones are being thrown away." << std::endl;
      params_ = params_in;
      obs_.clear();
      clone_info_.clear();
      clone_status_.clear();
      clone_master_.clear();
      running_.clear();
      suspended_.clear();
      finished_.clear();
    } else if (params_in.defined("NUM_CLONES")) {
      if ((params_.defined("NUM_CLONES") && params_in["NUM_CLONES"] != params_["NUM_CLONES"]) ||
          !params_.defined("NUM_CLONES")) {
        std::cout << "Info: number of clones in " << logger::task(task_id_) << " has been changed."
                  << std::endl;
        params_["NUM_CLONES"] = params_in["NUM_CLONES"];
        changed = true;
      }
    }

    if (changed) {
      num_clones_ = range_type(params_.value_or_default("NUM_CLONES", "1"), params_);
      progress_ = calc_progress();
      boost::tie(weight_, dump_weight_) = calc_weight();
      status_ = calc_status();
      save(opt);
    }
    halt();
  }
}

bool task::can_dispatch() const {
  return num_suspended() > 0 || num_started() < num_clones_.max BOOST_PREVENT_MACRO_SUBSTITUTION ();
}

void task::info_updated(cid_t cid, clone_info const& info) {
  if (clone_status_[cid] == clone_status::Running) {
    clone_info_[cid] = info;
    if (info.progress() >= 1) clone_status_[cid] = clone_status::Idling;
  }
}

void task::clone_halted(cid_t cid, clone_info const& info) {
  info_updated(cid, info);
  clone_halted(cid);
}

void task::clone_halted(cid_t cid) {
  if (clone_status_[cid] != clone_status::Stopping)
    boost::throw_exception(std::logic_error("clone is not stopping"));
  clone_status_[cid] = clone_status::Finished;
  running_.erase(cid);
  finished_.insert(cid);
  progress_ = calc_progress();
  status_ = calc_status();
  boost::tie(weight_, dump_weight_) = calc_weight();

  // clear cached observables (read from XML)
  obs_.clear();
}

void task::evaluate(alps::parapack::option const& opt) {
  std::cout << "evaluating " << file_out_str() << std::endl;

  if (!on_memory()) load();

  // bool same_weight = params_.defined("EVALUATE_CLONES_WITH_SAME_WEIGHT") &&
  //   static_cast<bool>(alps::evaluate("EVALUATE_CLONES_WITH_SAME_WEIGHT", params_));
  bool only_finished = params_.defined("EVALUATE_ONLY_FINISHED_CLONES") &&
    static_cast<bool>(alps::evaluate("EVALUATE_ONLY_FINISHED_CLONES", params_));

  obs_.clear();
  std::vector<cid_t> clones;
  BOOST_FOREACH(cid_t cid, finished_) clones.push_back(cid);
  if (!only_finished) {
    BOOST_FOREACH(cid_t cid, running_) clones.push_back(cid);
    BOOST_FOREACH(cid_t cid, suspended_) clones.push_back(cid);
  }
  std::sort(clones.begin(), clones.end());

  alps::Parameters p = params_;
  p["DIR_NAME"] = basedir_.string();
  p["BASE_NAME"] = base_;

  boost::shared_ptr<parapack::abstract_evaluator> evaluator
    = parapack::evaluator_factory::instance()->make_evaluator(p);

  std::cout << "  loading clones: ";
  std::vector<alps::alea::observable_set> merged;
  BOOST_FOREACH(cid_t cid, clones) {
    std::cout << (cid+1) << ' ' << std::flush;
    const std::size_t workers = clone_info_[cid].checkpoints().size();
    for (std::size_t w = 0; w < workers; ++w) {
      std::filesystem::path dump_h5 =
        alps::fs::absolute(std::filesystem::path(clone_info_[cid].checkpoints()[w] + ".h5"), basedir_);
      std::vector<alps::alea::observable_set> o;
      bool success;
      #pragma omp critical (hdf5io)
      {
        hdf5::archive h5(dump_h5);
        success = workers == 1 ? load_observable(h5, cid, o) : load_observable(h5, cid, w, o);
      }
      if (success) {
        std::vector<alps::alea::result_set> r;
        BOOST_FOREACH(alps::alea::observable_set const& m, o) r.push_back(m.results());
        evaluator->load(r, merged);
      }
    }
  }
  std::cout << std::endl;
  if (clones.size() > 0) {
    BOOST_FOREACH(alps::alea::observable_set const& m, merged) obs_.push_back(m.results());
    evaluator->evaluate(obs_);
    save_observable(opt);
  }
  halt();
}

void task::write_xml_summary(oxstream& os) const {
  os << start_tag("TASK")
     << attribute("id", task_id_+1)
     << attribute("status", task_status::to_string(status_))
     << attribute("progress", precision(progress() * 100, 3) + '%')
     << attribute("weight", precision(dump_weight_, 3))
     << start_tag("INPUT")
     << attribute("file", file_in_str_)
     << end_tag("INPUT")
     << start_tag("OUTPUT")
     << attribute("file", file_out_str_)
     << end_tag("OUTPUT")
     << end_tag("TASK");
}

void task::write_xml_archive(oxstream& os) const {
  os << alps::start_tag("SIMULATION")
     << attribute("id", task_id_+1)
     << attribute("status", task_status::to_string(status_))
     << attribute("progress", precision(progress() * 100, 3) + '%')
     << params_;
  if (obs_.size() == 1) {
    scheduler::write_xml_averages(os, obs_[0]);
  } else {
    for (unsigned int i = 0; i < obs_.size(); ++i) scheduler::write_xml_averages(os, obs_[i], i+1);
  }
  for (unsigned int i = 0; i < clone_info_.size(); ++i) if (clone_info_[i].clone_id() == i) os << clone_info_[i];
  os << alps::end_tag("SIMULATION");
}

double task::calc_progress() const {
  if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
  return (double)(num_finished()) / num_clones().min BOOST_PREVENT_MACRO_SUBSTITUTION ();
}

std::pair<double, double> task::calc_weight() const {
  if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
  double w;
  if (num_suspended() > 0)
    w = 4.0;
  else if (num_started() == 0)
    // NotStarted
    w = 3.0;
  else if (num_started() < num_clones().min BOOST_PREVENT_MACRO_SUBSTITUTION ())
    // Running
    w = 2.0 - (double)(num_started()) / num_clones().min BOOST_PREVENT_MACRO_SUBSTITUTION ();
  else
    // Continuing
    w = 1.0 - (double)(num_started()) / num_clones().max BOOST_PREVENT_MACRO_SUBSTITUTION ();
  double d = (num_running() + num_suspended() > 0) ? 4.0 : w;
  return std::make_pair(w, d);
}

task_status_t task::calc_status() const {
  if (!on_memory()) boost::throw_exception(std::logic_error("task not loaded"));
  if (num_started() == 0)
    return task_status::Ready;
  else if (num_finished() < num_clones_.min BOOST_PREVENT_MACRO_SUBSTITUTION ())
    return task_status::Running;
  else if (num_finished() < num_clones_.max BOOST_PREVENT_MACRO_SUBSTITUTION ())
    return task_status::Continuing;
  else
    return task_status::Idling;
}

} // end namespace alps
