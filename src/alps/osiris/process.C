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

#include <alps/config.h>

#ifdef ALPS_HAVE_MPI
# undef SEEK_SET
# undef SEEK_CUR
# undef SEEK_END  
# include <mpi.h>
#endif
#include <alps/osiris/comm.h>
#include <alps/osiris/process.h>
#include <alps/osiris/dump.h>
#include <string>
#include <algorithm>
#include <functional>

namespace alps {

Process::Process(int i)
 : tid(i)
{
}

void Process::save(ODump& dump) const
{
  dump << tid;
}


void Process::load(IDump& dump)
{
  dump >> tid;
}


bool Process::local() const
{
  return (tid==detail::local_id());
}

bool Process::valid() const
{
#ifdef ALPS_HAVE_MPI

  int total;
  MPI_Comm_size(MPI_COMM_WORLD,&total);
  return ((tid>=0) && (tid < total));

#else

  return (tid==0);

#endif
}

} // end namespace alps
