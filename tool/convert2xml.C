/***************************************************************************
* ALPS++ library
*
* tool/convert2xml.C   convert old scheduler files to XML
*
* $Id$
*
* Copyright (C) 2002-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
                             Simon Trebst <trebst@itp.phys.ethz.ch>
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
**************************************************************************/

#include <alps/osiris.h>
#include <alps/scheduler.h>
#include <boost/filesystem/operations.hpp>
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
  std::ofstream out (jobname.c_str());
  out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<?xml-stylesheet type=\"text/xsl\" href=\"http://xml.comp-phys.org/2002/10/job.xsl\"?>\n"
      << "<JOB xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
      << "xsi:noNamespaceSchemaLocation=\"http://xml.comp-phys.org/2002/10/job.xsd\">\n";
  out << "  <OUTPUT file=\"" << outfilename+".out.xml" << "\"/>\n";
  for (int i=0;i<list.size();++i) {
    std::string outname = outfilename;
    outname +=".task" + boost::lexical_cast<std::string,int>(i+1);
    std::string inname = outname + ".in.xml";
    outname+=".out.xml";
    out << "  <TASK status=\"new\">\n";
    out << "    <INPUT file=\"" << inname << "\"/>\n";
    out << "    <OUTPUT file=\"" << outname << "\"/>\n";
//      out << "    <CPUS min=\"1\">\n";
    out << "  </TASK>\n";
    std::ofstream task (inname.c_str());
    task << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
	 << "<?xml-stylesheet type=\"text/xsl\" href=\"http://xml.comp-phys.org/2002/10/QMCXML.xsl\"?>\n"
	 << "<SIMULATION xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
	 << "xsi:noNamespaceSchemaLocation=\"http://xml.comp-phys.org/2002/10/QMCXML.xsd\">\n";
    list[i].write_xml(task);
    task << "</SIMULATION>\n";
  }
  out << "</JOB>\n";
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
  std::ofstream out (jobname.c_str());
  out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<?xml-stylesheet type=\"text/xsl\" href=\"http://xml.comp-phys.org/2002/10/QMCXML.xsl\"?>\n"
      << "<SIMULATION xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
      << "xsi:noNamespaceSchemaLocation=\"http://xml.comp-phys.org/2002/10/QMCXML.xsd\">\n";
  int dummy_i;
  int version;
  int num;
  dump >> version; // version
  dump >> dummy_i;  // user version
  alps::Parameters parms;
  dump >> parms;
  parms.write_xml(out);
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
    out << "  <MCRUN><CHECKPOINT format=\"osiris\" file=\"" << dstname << "\"/></MCRUN>\n";
    convert_run(srcname,dstname);
  }
  out << "</SIMULATION>\n";
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
  std::ofstream out (jobname.c_str());
  out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      << "<?xml-stylesheet type=\"text/xsl\" href=\"http://xml.comp-phys.org/2002/10/job.xsl\"?>\n"
      << "<JOB xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
      << "xsi:noNamespaceSchemaLocation=\"http://xml.comp-phys.org/2002/10/job.xsd\">\n";
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
      if(boost::filesystem::exists(dumpname))
      {
        out << "  <TASK status=\"" << status_text[status[i]] << "\">\n";
        out << "    <INPUT file=\"" << xmlname+".xml" << "\"/>\n";
        out << "  </TASK>\n";
        convert_simulation(dumpname,xmlname);
      }
    }
  out << "</JOB>\n";
}

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
    
  if (argc<2 || argc>3) {
    std::cerr << "Usage: " << argv[0] << " inputfile [outputbasename]\n";
    std::exit(-1);
  }
  std::string name=argv[1];
  alps::IXDRFileDump dump(name);
  int type;
  dump >> type;
  switch (type) {
  case alps::scheduler::MCDump_scheduler:
    convert_scheduler(argv[1],argv[argc-1]);
    break;
  case alps::scheduler::MCDump_task:
    convert_simulation(argv[1],argv[argc-1]);
    break;
  case alps::scheduler::MCDump_run:
    convert_run(argv[1],argv[argc-1]);
    break;
  default:
    convert_params(argv[1],argv[argc-1]);
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
