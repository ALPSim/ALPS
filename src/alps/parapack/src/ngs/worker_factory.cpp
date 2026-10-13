/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>,
*                            Ryo Igarashi <rigarash@issp.u-tokyo.ac.jp>,
*                            Haruhiko Matsuo <halm@rist.or.jp>,
*                            Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Yuichi Motoyama <yomichi@looper.t.u-tokyo.ac.jp>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/ngs/parapack/worker_factory.h>

namespace alps {
namespace ngs_parapack {

//
// abstract_worker
//

abstract_worker::~abstract_worker() {}

void abstract_worker::load_worker(alps::hdf5::archive& ar) { this->load(ar); }

void abstract_worker::save_worker(alps::hdf5::archive& ar) const { this->save(ar); }

//
// worker_factory
//

worker_factory::worker_pointer_type worker_factory::make_worker(alps::params const& p) {
  if (!instance()->worker_creator_) {
    std::cerr << "Error: no algorithm registered\n";
    boost::throw_exception(std::runtime_error("worker_factory::make_creator()"));
  }
  return instance()->worker_creator_->create(p);
}

worker_factory* worker_factory::instance() {
  if (!instance_) instance_ = new worker_factory;
  return instance_;
}

#ifdef ALPS_HAVE_MPI

//
// parallel_worker_factory
//

parallel_worker_factory::worker_pointer_type
parallel_worker_factory::make_worker(boost::mpi::communicator const& comm,
  alps::params const& p) {
  if (!instance()->worker_creator_) {
    std::cerr << "Error: no algorithm registered\n";
    boost::throw_exception(std::runtime_error("worker_factory::make_creator()"));
  }
  return instance()->worker_creator_->create(comm, p);
}

parallel_worker_factory* parallel_worker_factory::instance() {
  if (!instance_) instance_ = new parallel_worker_factory;
  return instance_;
}

#endif // ALPS_HAVE_MPI

//
// initialization of static member pointer of factories
//

worker_factory* worker_factory::instance_ = 0;

#ifdef ALPS_HAVE_MPI

parallel_worker_factory* parallel_worker_factory::instance_ = 0;

#endif // ALPS_HAVE_MPI

} // end namespace ngs_parapack
} // end namespace alps
