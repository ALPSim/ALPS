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

#ifndef NGS_PARAPACK_PARAPACK_H
#define NGS_PARAPACK_PARAPACK_H

#include <alps/config.h>
#include <alps/ngs/parapack/worker_factory.h>
#include <alps/ngs/parapack/job.h>
#include <alps/parapack/option.h>
#include <iostream>

namespace alps {

struct no_worker {};
    
namespace ngs_parapack {

ALPS_DECL int start_impl(int argc, char **argv);

template<typename SERIAL, typename PARALLEL = alps::no_worker>
struct start;
  
template<typename SERIAL, typename PARALLEL>
struct start {
  start(int argc, char **argv) {
    worker_factory::register_worker<SERIAL>();
#ifdef ALPS_HAVE_MPI
    parallel_worker_factory::register_worker<PARALLEL>();
#endif
    ret = start_impl(argc, argv);
  }
  operator int() { return ret; }
private:
  int ret;
};
  
template<typename PARALLEL>
struct start<alps::no_worker, PARALLEL> {
  start(int argc, char **argv) {
#ifdef ALPS_HAVE_MPI
    parallel_worker_factory::register_worker<PARALLEL>();
#endif
    ret = start_impl(argc, argv);
  }
  operator int() { return ret; }
private:
  int ret;
};

template<typename SERIAL>
struct start<SERIAL, alps::no_worker> {
  start(int argc, char **argv) {
    worker_factory::register_worker<SERIAL>();
    ret = start_impl(argc, argv);
  }
  operator int() { return ret; }
private:
  int ret;
};

ALPS_DECL int run_sequential(int argc, char **argv);

ALPS_DECL int run_sequential_mpi(int argc, char **argv);

ALPS_DECL int start_sgl(int argc, char **argv);

ALPS_DECL int start_mpi(int argc, char **argv);

ALPS_DECL void print_copyright(std::ostream& os = std::cout);

ALPS_DECL void print_license(std::ostream& os = std::cout);

ALPS_DECL std::string alps_version();

ALPS_DECL void print_taskinfo(std::ostream& os, std::vector<alps::ngs_parapack::task> const& tasks,
  task_range_t const& task_range);

// return 1 for job XML (<JOB>) file or 2 for task XML (<SIMULATION>)
ALPS_DECL int load_filename(boost::filesystem::path const& file, std::string& file_in_str,
  std::string& file_out_str);

ALPS_DECL void load_version(boost::filesystem::path const& file,
  std::vector<std::pair<std::string, std::string> >& versions);

ALPS_DECL void load_tasks(boost::filesystem::path const& file_in,
  boost::filesystem::path const& file_out, boost::filesystem::path const& basedir,
  std::string& simname, std::vector<alps::ngs_parapack::task>& tasks, bool check_parameter, bool write_xml);

ALPS_DECL void save_tasks(boost::filesystem::path const& file, std::string const& simname,
  std::string const& file_in_str, std::string const& file_out_str, std::vector<alps::ngs_parapack::task>& tasks);

} // namespace ngs_parapack
} // namespace alps

#endif // NGS_PARAPACK_PARAPACK_H
