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

#ifndef NGS_PARAPACK_CLONE_H
#define NGS_PARAPACK_CLONE_H

#include <alps/ngs/parapack/worker_factory.h>
#include <alps/ngs/parapack/clone_info.h>
#include <alps/parapack/clone_timer.h>
#include <alps/parapack/types.h>

namespace alps {
namespace ngs_parapack {

class abstract_clone : public boost::noncopyable {
public:
  virtual ~abstract_clone() {}
  virtual void run(boost::function<bool ()> const& stop_callback,
    boost::function<void (double)> const& progress_callback) = 0;
  virtual bool halted() const = 0;
  virtual clone_info const& info() const = 0;

  virtual void load() = 0;
  virtual void save() const = 0;
  virtual void load(alps::hdf5::archive &) = 0;
  virtual void save(alps::hdf5::archive &) const = 0;

  virtual void checkpoint() = 0;
  virtual void suspend() = 0;
};

class clone : public abstract_clone {
public:
  clone(boost::filesystem::path const& basedir, dump_policy_t dump_policy, 
    clone_timer::duration_t const& check_interval, tid_t tid, cid_t cid, alps::params const& p,
    std::string const& base, bool is_new);
  virtual ~clone();

  tid_t task_id() const { return task_id_; }
  cid_t clone_id() const { return clone_id_; }

  void run(boost::function<bool ()> const& stop_callback,
    boost::function<void (double)> const& progress_callback);
  bool halted() const;
  clone_info const& info() const;

  void load();
  void save() const;
  void load(alps::hdf5::archive & ar);
  void save(alps::hdf5::archive & ar) const;

  void checkpoint();
  void suspend();

  void output() const;

protected:
  void do_halt();

private:
  tid_t task_id_;
  cid_t clone_id_;

  alps::params params_;

  boost::filesystem::path basedir_;
  clone_info info_;

  dump_policy_t dump_policy_;
  clone_timer timer_;
  clone_timer::loops_t loops_;

  boost::shared_ptr<alps::ngs_parapack::abstract_worker> worker_;
};

#ifdef ALPS_HAVE_MPI

struct clone_create_msg_t {
  clone_create_msg_t() {}
  clone_create_msg_t(tid_t tid, cid_t cid, gid_t gid, alps::params const& p, std::string const& bs,
    bool in) : task_id(tid), clone_id(cid), group_id(gid), p(p), base(bs), is_new(in) {}
  tid_t task_id;
  cid_t clone_id;
  gid_t group_id;
  alps::params p;
  std::string base;
  bool is_new;
  BOOST_SERIALIZATION_SPLIT_MEMBER()
  template<class Archive>
  void save(Archive & ar, const unsigned int) const {
    ar & task_id & clone_id & group_id & p & base & is_new;
  }
  template<class Archive>
  void load(Archive & ar, const unsigned int) {
    ar & task_id & clone_id & group_id & p & base & is_new;
  }
};

struct clone_info_msg_t {
  clone_info_msg_t() {}
  clone_info_msg_t(tid_t tid, cid_t cid, gid_t gid, clone_info const& f)
    : task_id(tid), clone_id(cid), group_id(gid), info(f) {}
  tid_t task_id;
  cid_t clone_id;
  gid_t group_id;
  clone_info info;
  template<typename Archive>
  void serialize(Archive & ar, const unsigned int) { ar & task_id & clone_id & group_id & info; }
};

struct clone_halt_msg_t {
  clone_halt_msg_t() {}
  clone_halt_msg_t(tid_t tid, cid_t cid, gid_t gid) : task_id(tid), clone_id(cid), group_id(gid) {}
  tid_t task_id;
  cid_t clone_id;
  gid_t group_id;
  template<typename Archive>
  void serialize(Archive & ar, const unsigned int) { ar & task_id & clone_id & group_id; }
};

class clone_mpi : public abstract_clone {
public:
  clone_mpi(boost::mpi::communicator const& ctrl, boost::mpi::communicator const& work,
    boost::filesystem::path const& basedir, dump_policy_t dump_policy,
    clone_timer::duration_t const& check_interval, clone_create_msg_t const& msg);
  virtual ~clone_mpi();

  tid_t task_id() const { return task_id_; }
  cid_t clone_id() const { return clone_id_; }

  void run(boost::function<bool ()> const& stop_callback,
    boost::function<void (double)> const& progress_callback);
  bool halted() const;
  clone_info const& info() const;

  void load();
  void save() const;
  void load(alps::hdf5::archive & ar);
  void save(alps::hdf5::archive & ar) const;

  void checkpoint();
  void suspend();

  void output() const;

protected:
  void do_checkpoint();
  void do_suspend();
  void do_halt();

  void send_info(mcmp_tag_t tag);
  void send_halted();

private:
  boost::mpi::communicator ctrl_, work_;

  tid_t task_id_;
  cid_t clone_id_;
  gid_t group_id_;

  alps::params params_;

  boost::filesystem::path basedir_;
  clone_info info_;

  dump_policy_t dump_policy_;
  clone_timer timer_;
  clone_timer::loops_t loops_;

  boost::shared_ptr<alps::ngs_parapack::abstract_worker> worker_;
};

#endif // ALPS_HAVE_MPI

} // end namespace ngs_parapack
} // end namespace alps

#endif // NGS_PARAPACK_CLONE_H
