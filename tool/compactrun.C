/***************************************************************************
* ALPS++/scheduler library
*
* tool/compactrun.C   compact a run checkpoint 
*
*
* $Id$
*
* Copyright (C) 2002 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

#include <alps/scheduler.h>
#include <boost/filesystem/operations.hpp>

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  if (argc<2 || argc>3) {
    std::cerr << "Usage: " << argv[0] << " input [output]\n";
    std::exit(-1);
  }
  std::string inname=argv[1];
  std::string outname = argv[argc-1];
  
  boost::filesystem::path inpath(inname, boost::filesystem::native);
  boost::filesystem::path outpath(outname, boost::filesystem::native);
  
  bool make_backup = boost::filesystem::exists(outpath);
  boost::filesystem::path writepath = (make_backup ? outpath.branch_path()/(outpath.leaf()+".bak") : outpath);
  
  std::cout << "Compacting run file " << inname << " to " <<  outname
	    <<std::endl;

  { // scope for life time of files
    alps::IXDRFileDump in(inpath);
    alps::OXDRFileDump out(writepath);
    alps::scheduler::DummyMCRun run;
    run.load_worker(in);
    run.save_worker(out);
  }

  if (make_backup) {
    boost::filesystem::remove(outpath);
    boost::filesystem::rename(writepath,outpath);
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
