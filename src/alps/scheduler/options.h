/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

#ifndef ALPS_SCHEDULER_OPTIONS_H
#define ALPS_SCHEDULER_OPTIONS_H

#include <alps/config.h>
#include <boost/filesystem/path.hpp>
#include <string>

namespace alps {
namespace scheduler {

//=======================================================================
// Options
//
// a class containing the options set by the user, either via command
// line switches or environment variables
//-----------------------------------------------------------------------

class ALPS_DECL NoJobfileOptions
{
public:
  std::string programname;    // name of the executable
  double min_check_time;      // minimum time between checks
  double max_check_time;      // maximum time between checks
  double checkpoint_time;     // time between two checkpoints
  int min_cpus;               // minimum number of runs per simulation
  int max_cpus;               // maximum number of runs per simulation
  double time_limit;          // time limit for the simulation
  bool use_mpi;               // should we use MPI
  bool valid;                 // shall we really run?
  bool write_xml;             // shall we write the results to XML?

  NoJobfileOptions(int argc, char** argv);
  NoJobfileOptions();
};

class Options : public NoJobfileOptions
{
public:
  boost::filesystem::path jobfilename;      // name of the jobfile

  Options(int argc, char** argv);
  Options();
};

} // end namespace
} // end namespace

#endif

