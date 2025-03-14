/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
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

/* $Id$ */

#include <alps/scheduler/task.h>
#include <alps/scheduler/types.h>
#include <alps/scheduler/scheduler.h>
#include <alps/expression.h>
#include <alps/parser/parser.h>
#include <alps/osiris/mpdump.h>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/throw_exception.hpp>
#include <fstream>
#include <stdexcept>

#define ALPS_TRACE

namespace alps {
namespace scheduler {

WorkerTask::WorkerTask(const ProcessList& w,const boost::filesystem::path& filename)
  : Task(w,filename),
    start_time(0),
    start_work(0.)
{
}

WorkerTask::WorkerTask(const ProcessList& w,const Parameters& p)
  : Task(w,p),
    start_time(0),
    start_work(0.)
{
}

WorkerTask::~WorkerTask()
{
  for (unsigned int i=0;i<runs.size();++i)
    if(runs[i])
      delete runs[i];
}

void WorkerTask::handle_tag(std::istream& infile, const XMLTag& intag)
{
  if (intag.name!=worker_tag()) {
    Task::handle_tag(infile,intag);
    return;
  }

  XMLTag tag=intag;
  // scan for <CHECKPOINT> tag
  if (tag.type==XMLTag::SINGLE)
    boost::throw_exception(std::runtime_error("<CHECKPOINT> element missing in task file"));
  std::string worker_close ="/"+worker_tag();
  tag=parse_tag(infile,true);
  while (tag.name!="CHECKPOINT") {
    if(tag.name==worker_close)
      boost::throw_exception(std::runtime_error("<CHECKPOINT> element missing in task file"));
    skip_element(infile,tag);
    tag=parse_tag(infile,true);
  }

  // read <CHECKPOINT> tag
  CheckpointFiles files;
  while (tag.name == "CHECKPOINT") {
    if (tag.attributes["file"]=="")
      boost::throw_exception(std::runtime_error("file attribute missing in <CHECKPOINT> element in task file"));
    if (tag.attributes["format"]=="osiris")
      files.in=boost::filesystem::absolute(
      boost::filesystem::path(tag.attributes["file"]),infilename.parent_path());
    else if (tag.attributes["format"]=="hdf5")
      files.hdf5in=boost::filesystem::absolute(
      boost::filesystem::path(tag.attributes["file"]),infilename.parent_path());
    else
      boost::throw_exception(std::runtime_error("unknown format in <CHECKPOINT> element in task file"));
    skip_element(infile,tag);
    tag=parse_tag(infile,true);
  }

    // relative to XML file

  runfiles.push_back(files);
  workerstatus.push_back(RunOnDump);
  while (tag.name!=worker_close) {
    skip_element(infile,tag);
    tag=parse_tag(infile,true);
  }
}

void WorkerTask::construct() // delayed until child class is fully constructed
{
  Task::construct();
  runs.resize(workerstatus.size());
  ProcessList here(cpus());
  int j=-1; // count existing runs
  int in=0; // first available node
  for (unsigned int i=0;i<runs.size();i++) {
    j++;
    // load as many runs as possible
    if(in+cpus()<=where.size()) {// a process is available
      if(j==0&&where[in].local()) {
        // one run runs locally
#ifdef ALPS_TRACE
        std::cerr  << "Loading run 1 locally on " << where[0] << "\n";
#endif
        std::copy(where.begin()+in,where.begin()+in+cpus(),here.begin());
        runs[0]=theScheduler->make_worker(here,parms);
        runs[0]->load_from_file(runfiles[i].in,runfiles[i].hdf5in);
        theWorker = runs[0];
        workerstatus[0] = LocalRun;
        in+=cpus();
      }
      else { // load other runs onto remote nodes
#ifdef ALPS_TRACE
        std::cerr  << "Loading run " << j+1 << " remote on " << where[j] << "\n";
#endif
        std::copy(where.begin()+in,where.begin()+in+cpus(),here.begin());
        runs[j]=new RemoteWorker(here,parms);
        runs[j]->load_from_file(runfiles[i].in,runfiles[i].hdf5in);
        workerstatus[j] = RemoteRun;
        in+=cpus();
      }
    }
    else { // no node available: load information only
#ifdef ALPS_TRACE
      std::cerr  << "Loading information about run " << j+1 << " from file "
                 << runfiles[i].in.string() << "\n";
#endif
      runs[j]=theScheduler->make_worker(parms);
      runs[j]->load_from_file(runfiles[i].in,runfiles[i].hdf5in);
      workerstatus[j] = RunOnDump;
    }
  }

  if(in+cpus()<=where.size()) { // more nodes than runs dumped: create extra runs
    runs.resize(where.size()/cpus());
    workerstatus.resize(where.size()/cpus());
    runfiles.resize(where.size()/cpus());
    for(int i=j+1;in+cpus()<=where.size();i++)
    {
      std::copy(where.begin()+in,where.begin()+in+cpus(),here.begin());
      if(in==0&&here[0].local()) { // one on the local node
        runs[0]=theScheduler->make_worker(here,parms);
        theWorker = runs[0];
        parms["SEED"] = static_cast<int32_t>(parms["SEED"])+cpus();
        in +=cpus();
        workerstatus[0] = LocalRun;
#ifdef ALPS_TRACE
        std::cerr  << "Created run 1 locally\n";
#endif
      }
      else { // other runs on remote nodes
        runs[i]=new RemoteWorker(here,parms);
        parms["SEED"] = static_cast<int32_t>(parms["SEED"])+cpus();
        in +=cpus();
        workerstatus[i] = RemoteRun;
#ifdef ALPS_TRACE
        std::cerr  << "Created run " << i+1 << " remote on Host ID: "
                   << where[i]<< "\n";
#endif
      }
    }
  }
  for (unsigned int i=0;i<runs.size();++i)
    runs[i]->set_parameters(parms);
}

// start all runs which are active
void WorkerTask::start()
{
  if(!started()) {
    Task::start();
    for (unsigned int i=0; i<runs.size();i++)
      if(runs[i] && workerstatus[i] > RunNotExisting && workerstatus[i] < RunOnDump) {
        runs[i]->start_worker();
      }
  }
}


// start an extra run on a new node
void WorkerTask::add_process(const Process& p)
{
  ProcessList here(1);
  here[0]=p;

  unsigned int i;
  // look for empty entry
  for ( i=0;i<where.size() && where[i].valid();i++)
    {}
  if(i==where.size())
    where.resize(i+1);
  where[i] = p;

  unsigned int j;
  // look for run to start on this process
  for (j=0; j<runs.size() && runs[j] && workerstatus[j] != RunNotExisting
                              && workerstatus[j] != RunOnDump ; j++)
    {}

  if(i != j)
    boost::throw_exception(std::logic_error( "In Task::add_process: # running runs != # running processes"));

  if(j==runs.size() || workerstatus[j] != RunOnDump) { // start new run
    runs.resize(j+1);
    workerstatus.resize(j+1);
    runfiles.resize(j+1);
#ifdef ALPS_TRACE
    std::cerr  << "Creating additional run " << j+1 << " remote on Host: " << p << "\n";
#endif
    runs[j]=new RemoteWorker(here,parms);
    parms["SEED"] = static_cast<int32_t>(parms["SEED"])+cpus();
    workerstatus[j] = RemoteRun;
    if(started())
      runs[j]->start_worker();
  } else {// continue old run
#ifdef ALPS_TRACE
    std::cerr  << "Loading additional run " << j << " remote on Host: " << p << "\n";
#endif
    runs[j]=new RemoteWorker(here,parms);
    runs[j]->load_from_file(runfiles[j].in,runfiles[i].hdf5in);
    workerstatus[j] = RemoteRun;
  }
}



// is it finished???
bool WorkerTask::finished(double& more_time, double& percentage) const
{
  if (finished_)
    return true;

  // get work estimate
  double w = work();
  if(w<=0.)
    return true;

  percentage = 1.-w;
  if (percentage < 0.)
    percentage=0.;
  else if (percentage>1.)
    percentage=1.;
  // estinate remaining time
  if(more_time<0)
    start_time=0; // new estimate

  if(start_time==0) { // initialize timing
    start_time=time(0);
    start_work=w;
    old_work=w;
  }
  else if(start_work==old_work) {
    start_time=time(0);
    if(w!=old_work) {
      start_work=w;
      old_work=-1;
    }
  }
  else if(start_work>w) {
    // estimate remaining time
    // propose to run 1/4 of that time
    time_t now = time(0);
    more_time = 0.25*w*(now-start_time)/(start_work-w);
  }
  return false;
}

// do some work on the local run
void WorkerTask::dostep()
{
  if(theWorker)
    dynamic_cast<Worker&>(*theWorker).run();
}

// halt all active runs
void WorkerTask::halt()
{
  if(started()) {
    Task::halt();
    for(unsigned int i=0;i<runs.size();i++)
      if(runs[i] && workerstatus[i] > RunNotExisting && workerstatus[i] < RunOnDump)
        runs[i]->halt_worker();
  }
}

ResultType WorkerTask::get_summary() const
{
  ResultType res;
  res.mean=0.;
  res.error=0.;
  res.count=0.;

  ProcessList where_master;

  // add runs stored locally
  if (runs.size()) {
    for (unsigned int i=0; i<runs.size(); i++) {
      if (workerstatus[i]==RemoteRun) {
        if (!runs[i])
          boost::throw_exception(std::runtime_error("Run does not exist in Task::get_measurements"));
        where_master.push_back(dynamic_cast<RemoteWorker*>(runs[i])->process());
      }
      else if (runs[i])
        res += runs[i]->get_summary();
    }
  }

  if (where_master.size()) {
    // broadcast request to all slaves
    OMPDump send;
    send.send(where_master,MCMP_get_summary);

    // collect results
    for (unsigned int i=0; i<where_master.size(); i++) {
      // receive dump
      IMPDump receive(MCMP_summary);
      ResultType s_res;
      receive >> s_res.T >> s_res.mean >> s_res.error >> s_res.count;
      res += s_res;
    }
  }
  return res;
}

double WorkerTask::work_done()  const
{
  double w=0.;
  ProcessList where_master;

  // add runs stored locally
  if(runs.size()) {
    for (unsigned int i=0;i<runs.size();i++) {
      if(workerstatus[i]==RemoteRun) {
         if(!runs[i])
            boost::throw_exception(std::runtime_error( "run does not exist in Task::get_measurements"));
        //where_master.push_back( Process(dynamic_cast<RemoteWorker&>(*runs[i]).process()));
        where_master.push_back(dynamic_cast<RemoteWorker*>(runs[i])->process());
      }
      else if(runs[i])
        w += runs[i]->work_done();
    }
  }

  // adding measurements from remote runs:
  if(where_master.size()) {
    // broadcast request to all slaves
    OMPDump send;
    send.send(where_master,MCMP_get_run_work);

    // collect results
    for (unsigned int i=0;i<where_master.size();i++) {
      // receive dump from remote process, abort if error
      IMPDump receive(MCMP_run_work);
      w += double(receive);
    }
  }
  return w;
}

double WorkerTask::work() const
{
  if (finished_)
    return 0.;
  return (parms.defined("WORK_FACTOR") ? alps::evaluate<double>(parms["WORK_FACTOR"], parms) : 1. )
         *(1.-work_done());
}

inline boost::filesystem::path optional_complete(boost::filesystem::path const& p, boost::filesystem::path const& dir)
{
  if (dir.empty() || p.empty())
    return p;
  else
    return boost::filesystem::absolute(p,dir);
}

// checkpoint: save into a file
void WorkerTask::write_xml_body(alps::oxstream& out, const boost::filesystem::path& fn, bool) const {
  boost::filesystem::path dir=fn.parent_path();
  for (unsigned int i=0;i<runs.size();++i) {
    if(workerstatus[i] == RunNotExisting) {
      if(runs[i])
        boost::throw_exception(std::logic_error("run exists but marked as non-existing"));
    }
    else if(runs[i]==0)
      boost::throw_exception(std::logic_error("run does not exist but marked as existing"));
    else {
      if (!runfiles[i].out.empty()) {
        runfiles[i].in=optional_complete(runfiles[i].out,dir);
#ifdef ALPS_HAVE_HDF5
        if (!runfiles[i].hdf5out.empty())
          runfiles[i].hdf5in=optional_complete(runfiles[i].hdf5out,dir);
        else {
          runfiles[i].hdf5in=runfiles[i].in.parent_path()/(runfiles[i].in.filename().string()+".h5");
        }
#endif
      }
      else {
        // search file name
        int j=0;
        bool found=false;
        std::string name;
        do {
          found = false;
          name =fn.filename().string();
          name = name.substr(0, name.find_last_of('.'));
          name+= ".run" + boost::lexical_cast<std::string,int>(j+1);
          for (unsigned int k=0;k<runfiles.size();++k)
          if(runfiles[k].out.filename().string()==name)
            found=true;
          j++;
        } while (found);
        runfiles[i].out = boost::filesystem::path(name);
#ifdef ALPS_HAVE_HDF5
        runfiles[i].hdf5out = boost::filesystem::path(name+".h5");
#endif
      }
      runfiles[i].in=optional_complete(runfiles[i].out,dir);
      if(workerstatus[i] == LocalRun || workerstatus[i] == RemoteRun) {
        runs[i]->save_to_file(optional_complete(runfiles[i].out,dir),optional_complete(runfiles[i].hdf5out,dir));
      } else if (workerstatus[i] == RunOnDump) {
        if(optional_complete(runfiles[i].out,dir).string()!=runfiles[i].in.string()) {
          boost::filesystem::remove(boost::filesystem::absolute(runfiles[i].out,dir));
          boost::filesystem::copy_file(boost::filesystem::absolute(runfiles[i].in,dir),boost::filesystem::absolute(runfiles[i].out,dir));
        }
#ifdef ALPS_HAVE_HDF5
        if(optional_complete(runfiles[i].hdf5out,dir).string()!=runfiles[i].hdf5in.string()) {
          if (boost::filesystem::exists(optional_complete(runfiles[i].hdf5out,dir)))
            boost::filesystem::remove(optional_complete(runfiles[i].hdf5out,dir));
          if (boost::filesystem::exists(optional_complete(runfiles[i].hdf5in,dir)))
            boost::filesystem::copy_file(optional_complete(runfiles[i].hdf5in,dir),optional_complete(runfiles[i].hdf5out,dir));
        }
#endif
      }
      else
        boost::throw_exception(std::logic_error("incorrect status of run"));
#ifdef ALPS_HAVE_HDF5
      if (!runfiles[i].hdf5out.empty())
        runfiles[i].hdf5in=optional_complete(runfiles[i].hdf5out,dir);
      else
        runfiles[i].hdf5in=boost::filesystem::path();
#endif
    }
  }
  ProcessList remote_runs;
  for (unsigned int i = 0; i < runs.size(); ++i) {
    if(workerstatus[i] == RemoteRun)
      remote_runs.push_back( Process(dynamic_cast<const RemoteWorker*>(runs[i])->process()));
    else {
      out << alps::start_tag(worker_tag())
          << runs[i]->get_info()
          << alps::start_tag("CHECKPOINT") << alps::attribute("format","osiris")
          << alps::attribute("file",runfiles[i].out.string())
          << alps::end_tag("CHECKPOINT");
#ifdef ALPS_HAVE_HDF5
      if (boost::filesystem::exists(runfiles[i].hdf5in)) {
        out << alps::start_tag("CHECKPOINT")
            << alps::attribute("format","hdf5")
            << alps::attribute("file",runfiles[i].hdf5out.string())
            << alps::end_tag("CHECKPOINT");
      }
#endif
      out << alps::end_tag(worker_tag());
    }
  }
  if(remote_runs.size()) {
    OMPDump send;
    send.send(remote_runs, MCMP_get_run_info);
    for (std::size_t i = 0; i < remote_runs.size(); ++i) {
      IMPDump receive(MCMP_run_info);
      TaskInfo infos;
      std::string name1, name2;
      receive >> infos >> name1 >> name2;
      out << alps::start_tag(worker_tag())
          << infos
          << alps::start_tag("CHECKPOINT")
          << alps::attribute("format","osiris")
          << alps::attribute("file",name1)
          << alps::end_tag("CHECKPOINT");
#ifdef ALPS_HAVE_HDF5
      if (boost::filesystem::exists(name2))
        out << alps::start_tag("CHECKPOINT")
            << alps::attribute("format","hdf5")
            << alps::attribute("file", name2)
            << alps::end_tag("CHECKPOINT");
#endif
      out << alps::end_tag(worker_tag());
    }
  }
}

} // namespace scheduler
} // namespace alps
