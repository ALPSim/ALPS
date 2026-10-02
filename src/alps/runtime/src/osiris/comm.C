/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
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

//=======================================================================
// INITIALIZATION AND CLEANUP
//
// initialize or stop the message passing library
//-----------------------------------------------------------------------

#ifdef ALPS_HAVE_MPI

namespace alps {
  int mpi_initialized=0;
}

void alps::comm_init(int& argc, char**& argv, bool usempi)
{
  if (usempi) {
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized)
      MPI_Init(&argc,&argv);
    mpi_initialized=1;
  }
  else
    mpi_initialized=0;
}

#else

void alps::comm_init(int&, char**&, bool usempi)
{
  if (usempi)
    boost::throw_exception(std::runtime_error("This program has not been compiled for use with MPI"));
}

#endif


// clean up everything
#ifdef ALPS_HAVE_MPI
void alps::comm_exit(bool kill_all)
{
  MPI_Initialized(&mpi_initialized);
  if (!mpi_initialized)
    return;
  if(kill_all)
    MPI_Abort(MPI_COMM_WORLD,-2);
  else
    MPI_Finalize();
}
#else
void alps::comm_exit(bool ) {}
#endif

//=======================================================================
// HOST/PROCESS ENQUIRIES
//
// ask for processes, hosts, ...
//-----------------------------------------------------------------------

// is this the master process ?

bool alps::is_master()
{
#ifdef ALPS_HAVE_MPI

  int num=0;
  if (mpi_initialized)
    MPI_Comm_rank(MPI_COMM_WORLD,&num);
  return (num==0);

#else
    return true; // only one CPU, always Master
#endif
}


// return an invalid host/process id

int alps::detail::invalid_id()
{
  return -1; // only one Process;
}



// return the id of the local process

int alps::detail::local_id()
{
#ifdef ALPS_HAVE_MPI
  int num=0;
  if (mpi_initialized)
    MPI_Comm_rank(MPI_COMM_WORLD,&num);
  return num;

#else
        return 0; // only one CPU, ID=0
#endif
}

// get a descriptor of this process

alps::Process alps::local_process()
{
#ifdef ALPS_HAVE_MPI
  int num=0;
  if (mpi_initialized)
    MPI_Comm_rank(MPI_COMM_WORLD,&num);
  return Process(num);
#else

  // single CPU case
  return Process(0);

#endif
}


// get a list of all processes running

alps::ProcessList alps::all_processes()
{
  ProcessList p;
#ifdef ALPS_HAVE_MPI

  int num=1;
  if (mpi_initialized)
    MPI_Comm_size(MPI_COMM_WORLD,&num);
  for (int i=0;i<num;i++)
    p.push_back(Process(i));

#else
  p.push_back(local_process());

#endif

  return p;
}


// get the parent of this process

alps::Process alps::master_process()
{
  return Process(0);
}


bool alps::runs_parallel()
{
#ifdef ALPS_HAVE_MPI
  return mpi_initialized;
#else
  return false;
#endif
}
