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
