/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>,
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

#include <alps/scheduler/scheduler.h>
#include <alps/scheduler/signal.hpp>
#include <alps/scheduler/types.h>
#include <alps/osiris/mpdump.h>
#include <alps/osiris/comm.h>
#include <alps/parser/parser.h>
#include <alps/parser/xmlstream.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/smart_ptr.hpp>
#include <fstream>

namespace alps {

namespace scheduler {

MasterScheduler::MasterScheduler(const NoJobfileOptions& opt,const Factory& p)
  : Scheduler(opt,p)
{
}


MasterScheduler::MasterScheduler(const Options& opt,const Factory& p)
  : Scheduler(opt,p)
{
  // the rest of the initialisation is done by the set_new_jobfile - function
  // this function can also be used to re-start the scheduler with a new
  // job file.
  set_new_jobfile(opt.jobfilename);
}

/**
 * Registers a new job file and updates everything that is related to the job
 * files - namely the task files are parsed and the tasks are created
 * accordingly.
 * This function is used for the fitting and allows to start another simulation
 * based on newly created files. Deleting and re-creating the schedulers yields
 * to problems with syncronisation and integrity of the message space.
 *
 * @params jobfilename The name of the new job file
 */
void MasterScheduler::set_new_jobfile(const boost::filesystem::path& jobfilename)
{
  // clean the 'traces' of the previous simulation
  taskfiles.clear();
  tasks.clear();
  taskstatus.clear();
  sim_results.clear();

  outfilepath = jobfilename;
  infilepath = jobfilename;
  
  infilepath=boost::filesystem::absolute(infilepath);
  outfilepath=boost::filesystem::absolute(outfilepath);

  if (!jobfilename.empty())
    parse_job_file(infilepath);
      
  tasks.resize(taskfiles.size());
  sim_results.resize(taskfiles.size());

  std::cout << "parsing task files ... \n";
  for (unsigned int i=0; i<taskfiles.size(); i++) {
#ifndef BOOST_NO_EXCEPTIONS
    try {
#endif
      tasks[i]=make_task(taskfiles[i].in);
      if (tasks[i] && 
          taskstatus[i]!= TaskFinished && tasks[i]->finished_notime())  {
        tasks[i]->start();
        std::cout << "Task " << i+1 << " is actually finished.\n";
        finish_task(i);
      }

#ifndef BOOST_NO_EXCEPTIONS
    }
    catch (const std::runtime_error& err) // file does not exist
    {
      std::cout << err.what() << "\n";
      std::cout  << "Cannot open simulation file " << taskfiles[i].in.string() 
                 << ".\n";
      tasks[i]=0;
      taskstatus[i] = TaskNotExisting;
    }
#endif
  }
  
}

void MasterScheduler::parse_job_file(const boost::filesystem::path& filename)
{
  boost::filesystem::ifstream infile(filename);
  XMLTag tag=parse_tag(infile,true);
  if (tag.name!="JOB")
    boost::throw_exception(std::runtime_error("missing <JOB> element in jobfile"));
  tag=parse_tag(infile);
  if (tag.name=="OUTPUT") {
    if(tag.attributes["file"]!="")
      outfilepath=boost::filesystem::absolute(
               boost::filesystem::path(tag.attributes["file"]),filename.parent_path());
    else
      boost::throw_exception(std::runtime_error(
               "missing 'file' attribute in <OUTPUT> element in jobfile"));
    tag=parse_tag(infile);
    if (tag.name=="/OUTPUT")
      tag=parse_tag(infile);
  }
  // make output path absolute
  while (tag.name=="TASK") {
    if (tag.attributes["status"]=="" || tag.attributes["status"]=="new")
      taskstatus.push_back(TaskNotStarted);
    else if (tag.attributes["status"]=="running")
      taskstatus.push_back(TaskHalted);
    else if (tag.attributes["status"]=="finished")
      taskstatus.push_back(TaskFinished);
    else
      boost::throw_exception(std::runtime_error(
               "illegal status attribute in <TASK> element in jobfile"));
    tag=parse_tag(infile);
    CheckpointFiles files;
    if (tag.name=="INPUT") {
      files.in=boost::filesystem::path(tag.attributes["file"]);
      if (files.in.empty())
        boost::throw_exception(std::runtime_error(
               "missing 'file' attribute in <INPUT> element in jobfile"));
      tag=parse_tag(infile);
      if (tag.name=="/INPUT")
        tag=parse_tag(infile);
    }
    else
      boost::throw_exception(std::runtime_error(
               "missing <INPUT> element in jobfile"));
    if (tag.name=="OUTPUT") {
      files.out=boost::filesystem::path(tag.attributes["file"]);
      if (files.out.empty())
        boost::throw_exception(std::runtime_error(
               "missing 'file' attribute in <OUTPUT> element in jobfile"));
      tag=parse_tag(infile);
      if (tag.name=="/OUTPUT")
        tag=parse_tag(infile);
    }
    if (files.out.empty())
      files.out=files.in;
    files.in=boost::filesystem::absolute(files.in,filename.parent_path());
    if (tag.name!="/TASK")
      boost::throw_exception(std::runtime_error(
               "missing </TASK> tag in jobfile"));
    tag = parse_tag(infile);
    taskfiles.push_back(files);
  }
  if (tag.name!="/JOB")
    boost::throw_exception(std::runtime_error("missing </JOB> tag in jobfile"));
}

// reload a simulation, this time to perform actual work
void  MasterScheduler::remake_task(ProcessList& where, const int i)
{
  if(tasks[i]==0)
    boost::throw_exception(std::logic_error(
               "cannot remake a simulation that does not exist"));
  delete tasks[i];
  tasks[i]=make_task(where,taskfiles[i].in);
}

MasterScheduler::~MasterScheduler()
{
  for (unsigned int i=0;i<tasks.size();++i)
    if(tasks[i])
      delete tasks[i];
}

void MasterScheduler::checkpoint()
{
  bool make_backup=boost::filesystem::exists(outfilepath);
  boost::filesystem::path filename=outfilepath;
  boost::filesystem::path dir=outfilepath.parent_path();
  if (make_backup)
    filename=dir/(filename.filename().string()+".bak");
  { // scope for out
    oxstream out(filename);

    out << header("UTF-8") << stylesheet(xslt_path("ALPS.xsl"));
    out << start_tag("JOB") 
        << xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
        << attribute("xsi:noNamespaceSchemaLocation",
               "http://xml.comp-phys.org/2003/8/job.xsd");
    int local_sim=-1;
    
    for (unsigned int i=0; i<tasks.size();i++) {
#ifdef ALPS_HAVE_HDF5
	    std::string task_path = taskfiles[i].out.string();
//		boost::filesystem::path task_path = taskfiles[i].out.string().substr(0, taskfiles[i].out.string().find_last_of('.')) + ".h5";
#else
		boost::filesystem::path task_path = taskfiles[i].out.string();
#endif

      if (taskstatus[i]==TaskFinished) {
        out << start_tag("TASK") << attribute("status","finished")
            << start_tag("INPUT") 
            << attribute("file",task_path)
            << end_tag() << end_tag();
        std::cout  << "Checkpointing Simulation " << i+1 << "\n";
        if (tasks[i]!=0 && boost::filesystem::absolute(taskfiles[i].out,dir).string()!=taskfiles[i].in.string()) {          
          tasks[i]->checkpoint(boost::filesystem::absolute(taskfiles[i].out,dir),write_xml);
          taskfiles[i].in=boost::filesystem::absolute(taskfiles[i].out,dir);
        }
        if (tasks[i]!=0) 
          delete tasks[i];
        tasks[i]=0;
      }
      else if(taskstatus[i]==TaskNotExisting) {
        out << start_tag("TASK") << attribute("status","finished")
            << start_tag("INPUT") << attribute("file",taskfiles[i].in.string())
            << end_tag() << end_tag();
        std::cout  << "Task# " << i+1 << " does not exist\n";
      } 
      else {
        out << start_tag("TASK") 
            << attribute("status",((taskstatus[i]==TaskNotStarted) ? "new" : "running"))
            << start_tag("INPUT") << attribute("file",task_path)
            << end_tag() << end_tag();
        if(theTask != tasks[i]) {
          std::cout  << "Checkpointing Simulation " << i+1 << "\n";
          tasks[i]->checkpoint(boost::filesystem::absolute(taskfiles[i].out,dir),write_xml);
                taskfiles[i].in=boost::filesystem::absolute(taskfiles[i].out,dir);
        }
        else
          local_sim=i;
      }
    }
    if(local_sim>=0) {
      std::cout  << "Checkpointing Simulation " << local_sim+1 << "\n";
      tasks[local_sim]->checkpoint(boost::filesystem::absolute(taskfiles[local_sim].out,dir),write_xml);
      taskfiles[local_sim].in=boost::filesystem::absolute(taskfiles[local_sim].out,dir);
    }
    out << end_tag("JOB");
  }
  if(make_backup) {
    boost::filesystem::remove(outfilepath);
    boost::filesystem::rename(filename,outfilepath);
  }
}


// store the results and delete the simulation
void MasterScheduler::finish_task(int i)
{ 
  if (tasks[i] == 0)
    return;
  tasks[i]->halt();
  taskstatus[i] = TaskHalted;
  std::cout  << "Halted Simulation " << i+1 << "\n";
  if (make_summary) {
    sim_results[i] = tasks[i]->get_summary();
  }
  tasks[i]->checkpoint(boost::filesystem::absolute(taskfiles[i].out,outfilepath.parent_path()),write_xml);
  delete tasks[i];
  tasks[i]=0;
  taskstatus[i] = TaskFinished;      
}

} // namespace scheduler
} // namespace alps
