/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

#ifndef OSIRIS_COMM_H
#define OSIRIS_COMM_H

#include <alps/osiris/process.h>
#include <alps/config.h>

namespace alps {


//=======================================================================
// INITIALIZATION AND CLEANUP
//
// initialize or stop the message passing library 
//-----------------------------------------------------------------------

// initialize everything

ALPS_DECL void comm_init(int& argc, char**& argv, bool=false);


// stop message passing
// the bool parameter indicates if all slave processes should be killed

ALPS_DECL void comm_exit(bool kill_slaves=false);


// do we actually run in parallel?

ALPS_DECL bool runs_parallel();

//=======================================================================
// HOST/PROCESS ENQUIRIES
//
// ask for processes, hosts, ... 
//-----------------------------------------------------------------------

namespace detail {
int local_id(); // return the id of this Process
int invalid_id(); // return an invalid id
}

ALPS_DECL bool is_master(); // is this the master Process ?

Process local_process(); // make a descriptor of the local Process
ProcessList all_processes(); // get a list of all running processes
Process master_process(); // get the master Process

} // end namespace alps

#endif // OSIRIS_COMM_H
