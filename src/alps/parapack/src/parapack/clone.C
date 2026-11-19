/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2014 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/
#include <alps/parapack/clone.h>
#include <alps/parapack/logger.h>
#include <alps/alea/observables_io.hpp>
#include <filesystem>
#include <alps/utilities/fs/absolute.hpp>

namespace alps {

namespace {

std::string clones_path() { return "simulation/realizations/0/clones"; }

void save_sets(alps::hdf5::archive & ar, std::string const& prefix,
               std::vector<alps::alea::observable_set> const& obs) {
  for (std::size_t m = 0; m < obs.size(); ++m) {
    std::string section = obs.size() == 1 ? "" : "/sections/" + boost::lexical_cast<std::string>(m);
    alps::alea::save_results(ar, prefix + "/results" + section, obs[m].results());
    alps::alea::save_checkpoint(ar, prefix + "/measurements" + section, obs[m]);
  }
}

// `suffix` is appended to the path of the only clone in the file
bool load_sets(alps::hdf5::archive & ar, cid_t cid, std::string const& suffix,
               std::vector<alps::alea::observable_set>& obs) {
  obs.clear();
  std::string prefix = clones_path();
  if (!ar.is_group(prefix)) return false;
  std::vector<std::string> children = ar.list_children(prefix);
  if (children.size() > 1) {
    std::cerr << logger::header() << "warning: more than one clone is found in in prefix = "
              << prefix << ". Only the first one will be loaded." << std::endl;
  }
  if (children.empty()) return false;
  if (children[0] != boost::lexical_cast<std::string>(cid)) {
    std::cerr << logger::header() << "warning: try to find clone id = " << cid << ", but found "
              << children[0] << " in prefix = " << prefix << std::endl;
  }
  std::string clone = prefix + "/" + children[0] + suffix;
  std::string measurements = clone + "/measurements";
  if (!ar.is_group(measurements)) {
    if (ar.is_group(clone + "/results"))
      boost::throw_exception(std::runtime_error(
        "The clone in '" + ar.get_filename() + "' was written by an earlier version of ALPS."
        " This version cannot continue it or evaluate it again. The results of a finished"
        " simulation can still be read from its result files."));
    return false;
  }
  if (!ar.is_group(measurements + "/sections/0")) {
    obs.resize(1);
    alps::alea::load_checkpoint(ar, measurements, obs[0]);
  } else {
    for (int m = 0; ; ++m) {
      std::string p = measurements + "/sections/" + boost::lexical_cast<std::string>(m);
      if (!ar.is_group(p)) break;
      obs.push_back(alps::alea::observable_set());
      alps::alea::load_checkpoint(ar, p, obs[m]);
    }
  }
  return true;
}

} // end namespace

void save_observable(alps::hdf5::archive & ar, cid_t cid,
                     std::vector<alps::alea::observable_set> const& obs) {
  save_sets(ar, clones_path() + "/" + boost::lexical_cast<std::string>(cid), obs);
}

void save_observable(alps::hdf5::archive & ar, cid_t cid, int rank,
                     std::vector<alps::alea::observable_set> const& obs) {
  save_sets(ar, clones_path() + "/" + boost::lexical_cast<std::string>(cid) +
            "/workers/" + boost::lexical_cast<std::string>(rank), obs);
}

bool load_observable(alps::hdf5::archive & ar, cid_t cid,
                     std::vector<alps::alea::observable_set>& obs) {
  return load_sets(ar, cid, "", obs);
}

bool load_observable(alps::hdf5::archive & ar, cid_t cid, int rank,
                     std::vector<alps::alea::observable_set>& obs) {
  return load_sets(ar, cid, "/workers/" + boost::lexical_cast<std::string>(rank), obs);
}


//
// clone
//

clone::clone(std::filesystem::path const& basedir, alps::parapack::option opt, tid_t tid,
  cid_t cid, Parameters const& params, std::string const& base, bool is_new)
  : task_id_(tid), clone_id_(cid), params_(params), basedir_(basedir),
    dump_policy_(opt.dump_policy), timer_(opt.check_interval) {
  params_["DIR_NAME"] = basedir_.string();
  params_["BASE_NAME"] = base;

  params_["TASK_ID"] = task_id_ + 1;
  params_["CLONE_ID"] = clone_id_ + 1;

  info_ = clone_info(clone_id_, params_, base);
  params_["WORKER_SEED"] = info_.worker_seed();
  params_["DISORDER_SEED"] = info_.disorder_seed();

  worker_ = parapack::worker_factory::make_worker(params_);
  if (!is_new) {
    bool exists = 
      std::filesystem::exists(alps::fs::absolute(std::filesystem::path(info_.dumpfile()), basedir_)) &&
      std::filesystem::exists(alps::fs::absolute(std::filesystem::path(info_.dumpfile_h5()), basedir_));
    if (exists) {
      this->load();
    } else {
      std::cerr << logger::header() << "warning: dump file not found. Restarting "
                << logger::clone(task_id_, clone_id_) << std::endl;
      is_new = true;
    }
  }
  if (is_new) worker_->init_observables(params_, measurements_);

  if (is_new && worker_->is_thermalized()) { // no thermalization steps
    BOOST_FOREACH(alps::alea::observable_set& m, measurements_) { m.reset(); }
  }

  if (is_new || worker_->progress() < 1)
    info_.start((worker_->is_thermalized()) ? "running" : "equilibrating");
  if (is_new && worker_->progress() >= 1) {
    info_.set_progress(worker_->progress());
    info_.stop();
    do_halt();
  }

  if (!is_new) timer_.reset(worker_->progress());
  loops_ = 1;
}

clone::~clone() {}

void clone::run() {
  for (clone_timer::loops_t i = 0; i < loops_; ++i) {
    bool thermalized = worker_->is_thermalized();
    double progress = worker_->progress();
    worker_->run(measurements_);
    if (!thermalized && worker_->is_thermalized()) {
      BOOST_FOREACH(alps::alea::observable_set& m, measurements_) m.reset();
      info_.stop();
      info_.start("running");
    }
    if (progress < 1 && worker_->progress() >= 1) {
      info_.set_progress(worker_->progress());
      info_.stop();
      do_halt();
      return;
    }
  }
  info_.set_progress(worker_->progress());
  loops_ = timer_.next_loops(loops_);
}

bool clone::halted() const { return !worker_; }

clone_info const& clone::info() const { return info_; }

void clone::load() {
  std::filesystem::path dump = alps::fs::absolute(std::filesystem::path(info_.dumpfile()), basedir_);
  std::filesystem::path dump_h5 =
    alps::fs::absolute(std::filesystem::path(info_.dumpfile_h5()), basedir_);
  if (exists(dump_h5)) {
    #pragma omp critical (hdf5io)
    {
      hdf5::archive h5(dump_h5.string());
      h5["/"] >> *this;
    }
  } else {
    boost::throw_exception(std::runtime_error(
      "The measurements of the clone are expected in '" + dump_h5.string() + "', which does not"
      " exist. Clones dumped in the xdr format by an earlier version of ALPS cannot be read."));
  }
  bool workerdump = (dump_policy_ == dump_policy::All) ||
    (dump_policy_ == dump_policy::RunningOnly && info_.progress() < 1);
  if (workerdump) {
    IXDRFileDump dp(dump);
    worker_->load_worker(dp);
  }
}

void clone::save() const{
  std::filesystem::path dump = alps::fs::absolute(std::filesystem::path(info_.dumpfile()), basedir_);
  std::filesystem::path dump_h5 =
    alps::fs::absolute(std::filesystem::path(info_.dumpfile_h5()), basedir_);
  #pragma omp critical (hdf5io)
  {
    hdf5::archive h5(dump_h5.string(), "a");
    h5["/"] << *this;
  }
  bool workerdump = (dump_policy_ == dump_policy::All) ||
    (dump_policy_ == dump_policy::RunningOnly && info_.progress() < 1);
  if (workerdump) {
    OXDRFileDump dp(dump);
    worker_->save_worker(dp);
  } else {
    if (exists(dump)) remove(dump);
  }
}

void clone::load(hdf5::archive & ar) {
  ar["parameters"] >> params_;
  ar["log/alps"] >> info_;
  load_observable(ar, clone_id_, measurements_);
}

void clone::save(hdf5::archive & ar) const {
  ar["parameters"] << params_;
  ar["log/alps"] << info_;
  save_observable(ar, clone_id_, measurements_);
}

void clone::checkpoint() {
  if (info_.progress() < 1) info_.stop();
  this->save();
}

void clone::suspend() {
  info_.stop();
  this->save();
  worker_.reset();
}

void clone::do_halt() {
  if (info_.progress() < 1)
    boost::throw_exception(std::logic_error("clone is not finished"));
  this->save();
  worker_.reset();
}

void clone::output() const{
  std::cout << params_;
  BOOST_FOREACH(alps::alea::observable_set const& m, measurements_) std::cout << m.results();
}

#ifdef ALPS_HAVE_MPI

//
// clone_mpi
//

clone_mpi::clone_mpi(alps::mpi::communicator const& ctrl, alps::mpi::communicator const& work,
  std::filesystem::path const& basedir, alps::parapack::option opt, clone_create_msg_t const& msg)
  : ctrl_(ctrl), work_(work), task_id_(msg.task_id), clone_id_(msg.clone_id),
    group_id_(msg.group_id), params_(msg.params), basedir_(basedir),
    dump_policy_(opt.dump_policy), timer_(opt.check_interval) {
  bool is_new = msg.is_new;
  params_["DIR_NAME"] = basedir;
  params_["BASE_NAME"] = msg.base;

  params_["TASK_ID"] = task_id_ + 1;
  params_["CLONE_ID"] = clone_id_ + 1;

  info_ = clone_info_mpi(work_, clone_id_, params_, msg.base);
  params_["WORKER_SEED"] = info_.worker_seed();
  params_["DISORDER_SEED"] = info_.disorder_seed();

  if (work_.size() > 1)
    worker_ = alps::parapack::parallel_worker_factory::make_worker(work_, params_);
  else
    worker_ = alps::parapack::worker_factory::make_worker(params_);
  if (!is_new) {
    bool exists = 
      std::filesystem::exists(alps::fs::absolute(std::filesystem::path(info_.dumpfile()), basedir_)) &&
      std::filesystem::exists(alps::fs::absolute(std::filesystem::path(info_.dumpfile_h5()), basedir_));
    int all_exist = exists;
    MPI_Allreduce(MPI_IN_PLACE, &all_exist, 1, MPI_INT, MPI_LAND, work_);
    exists = all_exist;
    if (exists) {
      this->load();
    } else {
      if (work_.rank() == 0) {
        std::cerr << logger::header() << "warning: dump file not found. Restarting "
                  << logger::clone(task_id_, clone_id_) << std::endl;
      }
      is_new = true;
    }
  }
  if (is_new) worker_->init_observables(params_, measurements_);

  if (is_new && worker_->is_thermalized()) { // no thermalization steps
    BOOST_FOREACH(alps::alea::observable_set& m, measurements_) { m.reset(); }
  }

  if (is_new || info_.progress() < 1)
    info_.start((worker_->is_thermalized()) ? "running" : "equilibrating");
  if (work_.rank() == 0 && is_new && worker_->progress() >= 1) {
    info_.set_progress(worker_->progress());
    info_.stop();
    if (group_id_ == 0) {
      do_halt();
    } else {
      send_info(mcmp_tag::clone_info);
    }
  }

  if (work_.rank() == 0 && !is_new) timer_.reset(worker_->progress());
  loops_ = 1;
}

clone_mpi::~clone_mpi() {}

void clone_mpi::run() {
  int tag = (info_.progress() < 1) ? mcmp_tag::do_step : mcmp_tag::do_nothing;
  if (work_.rank() == 0 && mpi::iprobe(ctrl_, 0, mpi::any_tag)) {
    mpi::status stat = mpi::recv(ctrl_, 0, mpi::any_tag);
    tag = stat.tag();
    if (tag == mcmp_tag::clone_suspend && info_.progress() >= 1)
      tag = mcmp_tag::clone_halt;
  }
  broadcast(work_, tag, 0);

  switch (tag) {
  case mcmp_tag::do_step :
    for (clone_timer::loops_t i = 0; i < loops_; ++i) {
      bool thermalized = worker_->is_thermalized();
      double progress = worker_->progress();
      worker_->run(measurements_);
      if (!thermalized && worker_->is_thermalized()) {
        BOOST_FOREACH(alps::alea::observable_set& m, measurements_) m.reset();
        if (work_.rank() == 0) {
          info_.stop();
          info_.start("running");
        }
      }
      if (progress < 1 && worker_->progress() >= 1) {
        if (work_.rank() == 0) {
          info_.set_progress(worker_->progress());
          info_.stop();
        }
        if (group_id_ == 0) {
          do_halt();
        } else {
          send_info(mcmp_tag::clone_info);
        }
        return;
      }
    }
    if (work_.rank() == 0) {
      info_.set_progress(worker_->progress());
      loops_ = timer_.next_loops(loops_);
    }
    broadcast(work_, loops_, 0);
    break;
  case mcmp_tag::clone_info :
    send_info(mcmp_tag::clone_info);
    break;
  case mcmp_tag::clone_checkpoint :
    do_checkpoint();
    send_info(mcmp_tag::clone_checkpoint);
    break;
  case mcmp_tag::clone_suspend :
    this->do_suspend();
    send_info(mcmp_tag::clone_suspend);
    break;
  case mcmp_tag::clone_halt :
    this->do_halt();
    send_halted();
    break;
  case mcmp_tag::do_nothing :
    break;
  default:
    if (work_.rank() == 0)
      std::cerr << "Warning: ignoring a message with an unknown tag " << tag << std::endl;
  }
}

void clone_mpi::checkpoint() {
  if (work_.rank() == 0) {
    int tag = mcmp_tag::clone_checkpoint;
    broadcast(work_, tag, 0);
    do_checkpoint();
  }
}

void clone_mpi::suspend() {
  if (work_.rank() == 0) {
    int tag = mcmp_tag::clone_suspend;
    broadcast(work_, tag, 0);
    do_suspend();
  }
}

void clone_mpi::do_checkpoint() {
  if (work_.rank() == 0 && info_.progress() < 1) info_.stop();
  this->save();
}

void clone_mpi::do_suspend() {
  if (work_.rank() == 0) info_.stop();
  this->save();
  work_.barrier();
  worker_.reset();
}

void clone_mpi::do_halt() {
  if (work_.rank() == 0 && info_.progress() < 1) {
    std::cerr << "clone is not finished\n";
    boost::throw_exception(std::logic_error("clone is not finished"));
  }
  this->save();
  work_.barrier();
  worker_.reset();
}

bool clone_mpi::halted() const { return !worker_; }

clone_info const& clone_mpi::info() const { return info_; }

void clone_mpi::load() {
  std::filesystem::path dump = alps::fs::absolute(std::filesystem::path(info_.dumpfile()), basedir_);
  std::filesystem::path dump_h5 =
    alps::fs::absolute(std::filesystem::path(info_.dumpfile_h5()), basedir_);
  if (exists(dump_h5)) {
    #pragma omp critical (hdf5io)
    {
      hdf5::archive h5(dump_h5.string());
      h5["/"] >> *this;
    }
  } else {
    boost::throw_exception(std::runtime_error(
      "The measurements of the clone are expected in '" + dump_h5.string() + "', which does not"
      " exist. Clones dumped in the xdr format by an earlier version of ALPS cannot be read."));
  }
  bool workerdump = (dump_policy_ == dump_policy::All) ||
    (dump_policy_ == dump_policy::RunningOnly && info_.progress() < 1);
  broadcast(work_, workerdump, 0);
  if (workerdump) {
    IXDRFileDump dp(alps::fs::absolute(std::filesystem::path(info_.dumpfile()), basedir_));
    worker_->load_worker(dp);
  }
}

void clone_mpi::save() const{
  std::filesystem::path dump = alps::fs::absolute(std::filesystem::path(info_.dumpfile()), basedir_);
  std::filesystem::path dump_h5 =
    alps::fs::absolute(std::filesystem::path(info_.dumpfile_h5()), basedir_);
  #pragma omp critical (hdf5io)
  {
    hdf5::archive h5(dump_h5.string(), "a");
    h5["/"] << *this;
  }
  bool workerdump = (dump_policy_ == dump_policy::All) ||
    (dump_policy_ == dump_policy::RunningOnly && info_.progress() < 1);
  broadcast(work_, workerdump, 0);
  if (workerdump) {
    OXDRFileDump dp(dump);
    worker_->save_worker(dp);
  } else {
    if (exists(dump)) remove(dump);
  }
}

void clone_mpi::load(hdf5::archive & ar) {
  ar["parameters"] >> params_;
  ar["log/alps"] >> info_;
  if (work_.size() == 1)
    load_observable(ar, clone_id_, measurements_);
  else
    load_observable(ar, clone_id_, work_.rank(), measurements_);
}

void clone_mpi::save(hdf5::archive & ar) const {
  ar["parameters"] << params_;
  ar["log/alps"] << info_;
  if (work_.size() == 1)
    save_observable(ar, clone_id_, measurements_);
  else
    save_observable(ar, clone_id_, work_.rank(), measurements_);
}

void clone_mpi::output() const{
  std::cout << params_;
  BOOST_FOREACH(alps::alea::observable_set const& m, measurements_) std::cout << m.results();
}

void clone_mpi::send_info(mcmp_tag_t tag) {
  if (work_.rank() == 0) {
    if (ctrl_.rank() == 0) std::cerr << "Error: sending to myself in clone_mpi::send_info()\n";
    mpi::send(ctrl_, 0, tag, clone_info_msg_t(task_id_, clone_id_, group_id_, info_));
  }
}

void clone_mpi::send_halted() {
  if (work_.rank() == 0) {
    if (ctrl_.rank() == 0) std::cerr << "Error: sending to myself in clone_mpi::send_halted()\n";
    mpi::send(ctrl_, 0, mcmp_tag::clone_halt, clone_halt_msg_t(task_id_, clone_id_, group_id_));
  }
}

#endif // ALPS_HAVE_MPI

} // end namespace alps
