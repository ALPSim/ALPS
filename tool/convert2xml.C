/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Simon Trebst <trebst@itp.phys.ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id$ */

#include <alps/osiris.h>
#include <alps/scheduler.h>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>
#include <fstream>
#include <stdexcept>

void convert_params(const std::string& inname, const std::string& outfilename)
{
  std::ifstream in(inname.c_str());
  alps::ParameterList list;
  in >> list;
  std::string jobname=outfilename+".in.xml";
  std::cout << "Converting parameter file " << inname << " to " <<  jobname << std::endl;
  alps::oxstream out (jobname.c_str());
  out << alps::header("UTF-8") << alps::stylesheet(alps::xslt_path("job.xsl"));
  out << alps::start_tag("JOB") << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
    << alps::attribute("xsi:noNamespaceSchemaLocation","http://xml.comp-phys.org/2003/8/job.xsd")
    << alps::start_tag("OUTPUT") << alps::attribute("file",outfilename+".out.xml") << alps::end_tag("OUTPUT");
  for (int i=0;i<list.size();++i) {
    std::string outname = boost::filesystem::path(outfilename,boost::filesystem::native).leaf();
    outname +=".task" + boost::lexical_cast<std::string,int>(i+1);
    std::string inname = outname + ".in.xml";
    std::string fullinname = outfilename + ".task" + boost::lexical_cast<std::string,int>(i+1) + ".in.xml";
    outname+=".out.xml";
    out << alps::start_tag("TASK") << alps::attribute("status","new")
      << alps::start_tag("INPUT") << alps::attribute("file",inname) << alps::end_tag("INPUT")
      << alps::start_tag("OUTPUT") << alps::attribute("file",outname) << alps::end_tag("OUTPUT")
      << alps::end_tag("TASK");
    //      out << "    <CPUS min=\"1\">\n";
    alps::oxstream task (inname.c_str());
    task << alps::header("UTF-8") << alps::stylesheet(alps::xslt_path("ALPS.xsl"));
    task << alps::start_tag("SIMULATION") << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
         << alps::attribute("xsi:noNamespaceSchemaLocation","http://xml.comp-phys.org/2002/10/QMCXML.xsd");
    task << list[i];
    task << alps::end_tag("SIMULATION");
  }
  out << alps::end_tag("JOB");
}

void convert_run(const std::string& inname, const std::string& outname)
{
  alps::IXDRFileDump dump(inname);
  std::cout << "Converting run file " << inname << " to " <<  outname+".xml" <<std::endl;
  alps::scheduler::DummyMCRun run;
  run.load_worker(dump);
  run.write_xml(outname,inname);
}

void convert_simulation(const std::string& inname, const std::string& outname)
{
  alps::IXDRFileDump dump(inname);
  if (static_cast<int>(dump)!=alps::scheduler::MCDump_task)
    boost::throw_exception(std::runtime_error("did not get a simulation on dump"));
  std::string jobname=outname+".xml";
  std::cout << "Converting simulation file " << inname << " to " <<  jobname << std::endl;
  alps::oxstream out (jobname.c_str());
  out << alps::header("UTF-8") << alps::stylesheet(alps::xslt_path("QMCXML.xsl"))
      << alps::start_tag("SIMULATION") << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
      << alps::attribute("xsi:noNamespaceSchemaLocation","http://xml.comp-phys.org/2002/10/QMCXML.xsd");
  int dummy_i;
  int version;
  int num;
  dump >> version; // version
  dump >> dummy_i;  // user version
  alps::Parameters parms;
  dump >> parms;
  out << parms;
  dump >> dummy_i; // nodes
  dump >> dummy_i; // seed
  dump >> num; // info size
  alps::scheduler::TaskInfo info;
  for (int i=0;i<num;++i)
    info.load(dump,version);
  // dump >> dummy_i; // flag if stored split
  num = static_cast<int>(dump);
  std::cout << num << " run(s)" << std::endl;
  for (int i=0;i<num;++i) {
    std::string srcname = inname+ ".run" + boost::lexical_cast<std::string,int>(i+1);
    std::string dstname = outname+ ".run" + boost::lexical_cast<std::string,int>(i+1);
    if (srcname!=dstname)
    {
      boost::filesystem::remove(dstname);
      boost::filesystem::copy_file(srcname,dstname);
    }
    out << alps::start_tag("MCRUN") << alps::start_tag("CHECKPOINT") 
        << alps::attribute("format","osiris") << alps::attribute("file=","dstname")
        << alps::end_tag("CHECKPOINT") << alps::end_tag("MCRUN");
    convert_run(srcname,dstname);
  }
  out << alps::end_tag("SIMULATION");
}
  
void convert_scheduler(const std::string& inname, const std::string& outname)
{
  std::map<int,std::string> status_text;
  status_text[alps::scheduler::MasterScheduler::TaskNotStarted]="new";
  status_text[alps::scheduler::MasterScheduler::TaskRunning]="running";
  status_text[alps::scheduler::MasterScheduler::TaskHalted]="running";
  status_text[alps::scheduler::MasterScheduler::TaskFromDump]="running";
  status_text[alps::scheduler::MasterScheduler::TaskFinished]="finished";

  alps::IXDRFileDump dump(inname);
  if (static_cast<int>(dump)!=alps::scheduler::MCDump_scheduler)
    boost::throw_exception(std::runtime_error("did not get scheduler on dump"));
  std::string jobname=outname+".xml";
  std::cout << "Converting scheduler file " << inname << " to " <<  jobname << std::endl;
  alps::oxstream out (jobname.c_str());
  out << alps::header("UTF-8") << alps::stylesheet(alps::xslt_path("job.xsl"))
    << alps::start_tag("JOB") << alps::xml_namespace("xsi","http://www.w3.org/2001/XMLSchema-instance")
    << alps::attribute("xsi:noNamespaceSchemaLocation","http://xml.comp-phys.org/2003/8/job.xsd");
  int dummy_i;
  double dummy_d;
  dump >> dummy_i; // version
  dump >> dummy_d;  // steptime
  alps::ParameterList list;
  dump >> list;
  std::vector<int> status;
  dump >> status;
  for (int i=0;i<list.size();++i) 
    if (status[i]) {
      std::string xmlname = outname;
      std::string dumpname = inname;
      xmlname += ".task" + boost::lexical_cast<std::string,int>(i+1);
      dumpname += "Sim" + boost::lexical_cast<std::string,int>(i+1);
      if(boost::filesystem::exists(dumpname)) {
        out << alps::start_tag("TASK") << alps::attribute("status",status_text[status[i]])
          << alps::start_tag("INPUT") << alps::attribute("file",xmlname+".xml")
          << alps::end_tag("INPUT") << alps::end_tag("TASK");
        convert_simulation(dumpname,xmlname);
      }
    }
   out << alps::end_tag("JOB");
}

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
    
  if (argc<2) {
    std::cerr << "Usage: " << argv[0] << " inputfile [inputfile ...]]\n";
    std::exit(-1);
  }
  for (int i=1;i<argc;++i) {    
    std::string inname=argv[i];
    if (inname.size() >= 2 && inname.substr(0, 2) == "./") inname.erase(0, 2);
    alps::IXDRFileDump dump(inname);
    int type;
    dump >> type;
    switch (type) {
    case alps::scheduler::MCDump_scheduler:
      convert_scheduler(inname,inname);
      break;
    case alps::scheduler::MCDump_task:
      convert_simulation(inname,inname);
      break;
    case alps::scheduler::MCDump_run:
      convert_run(inname,inname);
      break;
    default:
      convert_params(inname,inname);
    }
  }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  std::exit(-5);
}
#endif

}
