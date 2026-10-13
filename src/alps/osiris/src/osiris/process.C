/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
