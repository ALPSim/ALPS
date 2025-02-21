/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>
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
