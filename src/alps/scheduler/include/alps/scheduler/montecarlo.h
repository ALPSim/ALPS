/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_SCHEDULER_MONTECARLO_H
#define ALPS_SCHEDULER_MONTECARLO_H

#include <alps/scheduler/scheduler.h>
#include <alps/scheduler/task.h>
#include <alps/scheduler/worker.h>
#include <alps/model/model_helper.h>
#include <alps/model/sign.h>
#include <alps/lattice/graph_helper.h>
#include <alps/alea/observables.hpp>
#include <alps/config.h>

#include <alps/hdf5.hpp>

#include <memory>
#include <string>

namespace alps {
namespace scheduler {

/// One Monte Carlo run. An application derives from this class, declares its
/// observables in `measurements` and adds a sample to each in every sweep.
class ALPS_DECL MCRun : public Worker
{
public:
  static void print_copyright(std::ostream&);

  MCRun(const ProcessList&,const alps::Parameters&,int);

  void save_worker(ODump&) const;
  void load_worker(IDump&);
  
  void save(hdf5::archive &) const;
  void load(hdf5::archive &);

  virtual void save(ODump&) const;
  virtual void load(IDump&);

  void write_xml(const std::filesystem::path& name) const;
  const alps::alea::observable_set& get_measurements() const { return measurements;}

  /// What this run has measured so far
  alps::alea::result_set get_results() const;

  /// False for a run read from a file of an earlier version of ALPS: its
  /// results can be read, but it cannot go on measuring or be merged.
  bool can_continue() const { return !earlier_results_; }

  std::string work_phase();
  void run();
  virtual bool is_thermalized() const;
  virtual bool handle_message(const Process& runmaster,int32_t tag);
protected:
  alps::alea::observable_set measurements;
private:
  std::shared_ptr<const alps::alea::result_set> earlier_results_;
};


class ALPS_DECL DummyMCRun : public MCRun
{
public:
  DummyMCRun(const ProcessList& w,const alps::Parameters& p,int n);
  DummyMCRun();
  void dostep();
  double work_done() const;
// astreich, 06.23
  ResultType get_summary() const;
};


class ALPS_DECL MCSimulation : public WorkerTask
{        
public:
  MCSimulation(const ProcessList& w, const std::filesystem::path& p) 
      : WorkerTask(w,p) { 
    construct();
  }

  MCSimulation(const ProcessList& w, const Parameters& p) 
      : WorkerTask(w,p) 
  {
    construct();
  }

  static void print_copyright(std::ostream&) {}
  
  /// The results of all runs together, and those added with add_result
  alps::alea::result_set get_measurements() const;
  
  /// Add a result, such as a function of the measured ones, to those written
  /// with the simulation. It replaces a result of the same name.
  void add_result(const std::string& name, const alps::alea::observable_result& result);

//protected:
  virtual ResultType get_summary() const;
  virtual ResultType get_summary(const std::string) const;

  void save(hdf5::archive &) const;
  void load(hdf5::archive &);

private:
  std::string worker_tag() const;
  void write_xml_body(alps::oxstream&, std::filesystem::path const& fn, bool) const;
  virtual void handle_tag(std::istream&, const XMLTag&);
  // results stored with the task, and those added
  alps::alea::result_set stored_results_, added_results_;
};


/// The results as the <AVERAGES> element of the ALPS XML files. A non-negative
/// `id` is written as its attribute.
ALPS_DECL void write_xml_averages(oxstream&, const alps::alea::result_set& results, int id = -1);


template <class G=graph_helper<>::graph_type>
class LatticeMCRun : public MCRun, public graph_helper<G>
{
public:
  LatticeMCRun(const ProcessList& w,const alps::Parameters& p,int n)
   : MCRun(w,p,n), graph_helper<G>(parms)
  {}
};


template <class G=graph_helper<>::graph_type, class I=short>
class LatticeModelMCRun : public LatticeMCRun<G>, public model_helper<I>
{
public:  
  LatticeModelMCRun(const ProcessList& w, const alps::Parameters& p, int n,
                    bool issymbolic = false)
    : LatticeMCRun<G>(w, p, n),
      model_helper<I>(*this, LatticeMCRun<G>::parms, issymbolic)
  {}
  
  bool has_sign_problem() const 
  {
    return alps::has_sign_problem(model_helper<I>::model(), *this,
                                  LatticeMCRun<G>::parms);
  }
};


template <class WORKER>
class SimpleMCFactory : public BasicFactory<MCSimulation,WORKER>
{
public:
  SimpleMCFactory() {}
};

} // end namespace scheduler
} // end namespace alps

#endif
