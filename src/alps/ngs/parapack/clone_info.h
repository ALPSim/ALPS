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

#ifndef NGS_PARAPACK_CLONE_INFO_H
#define NGS_PARAPACK_CLONE_INFO_H

#include <alps/ngs/params.hpp>
#include <alps/parapack/process.h>
#include <alps/parapack/types.h>
#include <alps/hdf5/archive.hpp>
#include <alps/parser/xmlstream.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <vector>
#include <string>

namespace alps {
namespace ngs_parapack {

//
// clone_phase
//

class clone_phase_xml_handler;

class ALPS_DECL clone_phase {
public:
  clone_phase() {}
  clone_phase(std::vector<std::string> const& hosts, std::string const& user,
              std::string const& phase);
  void start(std::string const& phase);
  void stop();

  std::vector<std::string> const& hosts() const { return hosts_; }
  std::string const& user() const { return user_; }
  std::string const& phase() const { return phase_; }
  boost::posix_time::ptime const& start_time() const { return startt_; }
  boost::posix_time::ptime const& stop_time() const { return stopt_; }
  boost::posix_time::time_duration elapsed() const { return stopt_ - startt_; }

  void write_xml(oxstream& os) const;
  void save(alps::hdf5::archive& ar) const;
  void load(alps::hdf5::archive& ar);

  template<typename Archive>
  void serialize(Archive & ar, const unsigned int) {
    std::string start_str = to_simple_string(startt_);
    std::string stop_str = to_simple_string(stopt_);
    ar &  hosts_ & user_ & phase_ & start_str & stop_str;
    startt_ = boost::posix_time::time_from_string(start_str);
    stopt_ = boost::posix_time::time_from_string(stop_str);
  }

private:
  friend class clone_phase_xml_handler;

  std::vector<std::string> hosts_;
  std::string user_;
  std::string phase_;
  boost::posix_time::ptime startt_;
  boost::posix_time::ptime stopt_;
};

} // end namespace ngs_parapack
} // end namespace alps

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
namespace ngs_parapack {
#endif

ALPS_DECL alps::oxstream& operator<<(alps::oxstream& os, alps::ngs_parapack::clone_phase const& phase);

// ALPS_DECL alps::ODump& operator<<(alps::ODump& dump, alps::clone_phase const& phase);

// ALPS_DECL alps::IDump& operator>>(alps::IDump& dump, alps::clone_phase& phase);

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // namespace ngs_parapack
} // namespace alps
#endif

namespace alps {
namespace ngs_parapack {

//
// clone_info
//

class clone_info_xml_handler;

class ALPS_DECL clone_info {
public:
  clone_info();
  clone_info(cid_t cid);
  clone_info(cid_t cid, alps::params const& p, std::string const& dump, bool initialize = true);
  virtual ~clone_info() {}

  cid_t clone_id() const { return clone_id_; }
  std::string const& dumpfile() const { return dumpfiles_[0]; }
  std::string dumpfile_h5() const { return dumpfile() + ".h5"; }

  bool has_seed() const;
  seed_t worker_seed() const;
  seed_t disorder_seed() const;

  // the following member functions are meaningful only on master node
  void start(std::string const& phase);
  void stop();

  double progress() const { return progress_; }
  void set_progress(double p) { progress_ = p; }

  boost::posix_time::time_duration elapsed() const;

  std::vector<std::string> const& hosts() const { return phases_.back().hosts(); }
  std::string const& phase() const { return phases_.back().phase(); }
  std::vector<clone_phase> const& phases() const { return phases_; }
  std::vector<std::string> const& checkpoints() const { return dumpfiles_; }

  void write_xml(oxstream& os) const;

  // checkpoint filenames will not be saved/loaded on slave worker
  void save (alps::hdf5::archive& ar) const;
  void load (alps::hdf5::archive& ar);

  template<typename Archive>
  void serialize(Archive & ar, const unsigned int) {
    ar & clone_id_ & progress_ & phases_ & dumpfiles_ & worker_seed_ & disorder_seed_;
  }

protected:
  void init(alps::params const& p, std::string const& dump);
  virtual unsigned int num_processes() const { return 1; };
  virtual unsigned int process_id() const { return 0; }
  virtual void set_hosts(std::vector<std::string>& hosts, bool& is_master);

private:
  friend class clone_info_xml_handler;

  cid_t clone_id_;
  double progress_;
  std::vector<clone_phase> phases_;
  std::vector<std::string> dumpfiles_;

  std::vector<seed_t> worker_seed_;
  seed_t disorder_seed_;

  // list of current hostnames
  std::vector<std::string> hosts_;
  bool is_master_;
};

} // end namespace ngs_parapack
} // end namespace alps

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
namespace ngs_parapack {
#endif

ALPS_DECL alps::oxstream& operator<<(alps::oxstream& os, alps::ngs_parapack::clone_info const& info);

// ALPS_DECL alps::ODump& operator<<(alps::ODump& dp, alps::ngs_parapack::clone_info const& info);

// ALPS_DECL alps::IDump& operator>>(alps::IDump& dp, alps::ngs_parapack::clone_info& info);

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // namespace ngs_parapack
} // namespace alps
#endif

#ifdef ALPS_HAVE_MPI

namespace alps {
namespace ngs_parapack {

//
// clone_info_mpi
//

class clone_info_mpi : public clone_info {
public:
  // interprocess communication is required
  clone_info_mpi(boost::mpi::communicator const& comm, cid_t cid, alps::params const& p,
    std::string const& base);
private:
  boost::mpi::communicator comm_;
  virtual unsigned int num_processes() const;
  virtual unsigned int process_id() const;
  virtual void set_hosts(std::vector<std::string>& hosts, bool& is_master);
};

} // end namespace ngs_parapack
} // end namespace alps

#endif // ALPS_HAVE_MPI

#endif // NGS_PARAPACK_CLONE_INFO_H
