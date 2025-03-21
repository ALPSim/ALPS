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

#ifndef NGS_PARAPACK_CLONE_PROXY_H
#define NGS_PARAPACK_CLONE_PROXY_H

#include <alps/ngs/parapack/clone.h>

namespace alps {
namespace ngs_parapack {

class clone_proxy {
public:
  clone_proxy(clone*& clone_ptr, boost::filesystem::path const& basedir, dump_policy_t dump_policy,
    clone_timer::duration_t const& check_interval) : clone_ptr_(clone_ptr), basedir_(basedir),
    dump_policy_(dump_policy), check_interval_(check_interval) {}
  bool is_local(Process const&) const { return true; }
  void start(tid_t tid, cid_t cid, thread_group const&, alps::params const& p,
    std::string const& base, bool is_new) {
    clone_ptr_ = new clone(basedir_, dump_policy_, check_interval_, tid, cid, p, base, is_new);
  }
  clone_info const& info(Process const&) const {
    if (!clone_ptr_)
      boost::throw_exception(std::logic_error("clone_proxy::info()"));
    return clone_ptr_->info();
  }
  void checkpoint(Process const&) { if (clone_ptr_) clone_ptr_->checkpoint(); }
  void update_info(Process const&) const {}
  void suspend(Process const&) { if (clone_ptr_) clone_ptr_->suspend(); }
  void halt(Process const&) { /* if (clone_ptr_) clone_ptr_->halt(); */ }
  void destroy(Process const&) {
    if (clone_ptr_) {
      delete clone_ptr_;
      clone_ptr_ = 0;
    }
  }
private:
  clone*& clone_ptr_;
  boost::filesystem::path basedir_;
  dump_policy_t dump_policy_;
  clone_timer::duration_t check_interval_;
};

#ifdef ALPS_HAVE_MPI

class clone_proxy_mpi {
public:
  clone_proxy_mpi(clone_mpi*& clone_ptr, boost::mpi::communicator const& comm_ctrl,
    boost::mpi::communicator const& comm_work, boost::filesystem::path const& basedir,
    dump_policy_t dump_policy, clone_timer::duration_t const& check_interval)
    : clone_ptr_(clone_ptr), comm_ctrl_(comm_ctrl), comm_work_(comm_work), basedir_(basedir),
      dump_policy_(dump_policy), check_interval_(check_interval) {}
  bool is_local(Process const& proc) const { return proc == 0; }
  void start(tid_t tid, cid_t cid, process_group const& procs, alps::params const& p,
    std::string const& base, bool is_new) const {
    clone_create_msg_t msg(tid, cid, procs.group_id, p, base, is_new);
    bool worker_on_master = false;
    BOOST_FOREACH(Process p, procs.process_list) {
      if (p == 0)
        worker_on_master = true;
      else
        comm_ctrl_.send(p, mcmp_tag::clone_create, msg);
    }
    if (worker_on_master)
      clone_ptr_ =
        new clone_mpi(comm_ctrl_, comm_work_, basedir_, dump_policy_, check_interval_, msg);
  }

  clone_info const& info(Process const& proc) const {
    if (proc != 0 || !clone_ptr_) {
      std::cerr << "clone_proxy_mpi::info()\n";
      boost::throw_exception(std::logic_error("clone_proxy_mpi::info()"));
    }
    return clone_ptr_->info();
  }

  void checkpoint(Process const& proc) {
    if (proc == 0) {
      if (clone_ptr_) clone_ptr_->checkpoint();
    } else {
      comm_ctrl_.send(proc, mcmp_tag::clone_checkpoint);
    }
  }

  void update_info(Process const& proc) const {
    if (proc == 0) {
      // nothing to do
    } else {
      comm_ctrl_.send(proc, mcmp_tag::clone_info);
    }
  }

  void suspend(Process const& proc) {
    if (proc == 0) {
      if (clone_ptr_) clone_ptr_->suspend();
    } else {
      comm_ctrl_.send(proc, mcmp_tag::clone_suspend);
    }
  }

  void halt(Process const& proc) {
    if (proc == 0) {
      // nothing to do
    } else {
      comm_ctrl_.send(proc, mcmp_tag::clone_halt);
    }
  }

  void destroy(Process const& proc) {
    if (proc == 0 && clone_ptr_) {
      delete clone_ptr_;
      clone_ptr_ = 0;
    }
  }

private:
  clone_mpi*& clone_ptr_;
  boost::mpi::communicator comm_ctrl_, comm_work_;
  boost::filesystem::path basedir_;
  dump_policy_t dump_policy_;
  clone_timer::duration_t check_interval_;
};

#endif // ALPS_HAVE_MPI

} // end namespace ngs_parapack
} // end namespace alps

#endif // NGS_PARAPACK_CLONE_PROXY_H
