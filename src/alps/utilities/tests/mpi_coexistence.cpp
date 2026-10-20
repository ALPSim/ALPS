// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// ALPS's <alps/ngs/mpi.hpp> and ALPSCore's <alps/utilities/mpi.hpp> both add to
// namespace alps::mpi. Check that one program can use both, and that neither
// MPI environment wrapper finalizes MPI that it did not initialize.
#include <alps/config.h>  // alps/ngs/mpi.hpp tests ALPS_HAVE_MPI without including it
#include <alps/ngs/mpi.hpp>
#include <alps/utilities/mpi.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <gtest/gtest.h>

#include <functional>

namespace {

bool mpi_finalized()
{
    int finalized = 0;
    MPI_Finalized(&finalized);
    return finalized != 0;
}

}

TEST(UtilitiesMpiCoexistence, BothInterfacesWork)
{
    boost::mpi::communicator world;
    alps::mpi::communicator comm;
    ASSERT_GE(comm.size(), 2);
    ASSERT_EQ(world.rank(), comm.rank());
    const int rank = comm.rank();

    // ALPS's reduce, on a Boost.MPI communicator.
    int sum = 0;
    alps::mpi::reduce(world, rank + 1, sum, std::plus<int>(), 0);
    if (rank == 0)
        EXPECT_EQ(world.size() * (world.size() + 1) / 2, sum);

    // ALPSCore's broadcast and all_reduce, on its own communicator.
    int value = rank == 0 ? 42 : 0;
    alps::mpi::broadcast(comm, value, 0);
    EXPECT_EQ(42, value);
    EXPECT_EQ(comm.size(), alps::mpi::all_reduce(comm, 1, std::plus<int>()));
}

TEST(UtilitiesMpiCoexistence, EnvironmentsLeaveMpiToItsOwner)
{
    // The test's main() initialized MPI, so neither wrapper takes control of it.
    {
        alps::mpi::environment env;
    }
    EXPECT_FALSE(mpi_finalized());
    {
        boost::mpi::environment env;
    }
    EXPECT_FALSE(mpi_finalized());
    alps::mpi::communicator().barrier();
}
