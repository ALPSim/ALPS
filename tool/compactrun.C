/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>
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

#include <alps/scheduler/montecarlo.h>
#include <alps/osiris/xdrdump.h>
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
  
  boost::filesystem::path inpath(inname);
  boost::filesystem::path outpath(outname);
  
  bool make_backup = boost::filesystem::exists(outpath);
  boost::filesystem::path writepath = (make_backup ? outpath.parent_path()/(outpath.filename().string()+".bak") : outpath);
  
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
