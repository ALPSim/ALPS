/***************************************************************************
* ALPS++/scheduler library
*
* scheduler/options.h   A class to store options
*
* $Id$
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*
* This software is part of the ALPS library, published under the 
* ALPS Library License; you can use, redistribute it and/or modify 
* it under the terms of the License, either version 1 or (at your option) 
* any later version.
*
* You should have received a copy of the ALPS Library License along with 
* the ALPS Library; see the file License.txt. If not, the license is also 
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
**************************************************************************/

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

class Options
{
public:
  std::string programname;    // name of the executable
  boost::filesystem::path jobfilename;      // name of the jobfile
  double min_check_time;      // minimum time between checks
  double max_check_time;      // maximum time between checks
  double checkpoint_time;     // time between two checkpoints
  int min_cpus;               // minimum number of runs per simulation
  int max_cpus;               // maximum number of runs per simulation
  double time_limit;          // time limit for the simulation

  Options(int argc=0, char** argv=0);
};

} // end namespace
} // end namespace

#endif

