// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// ALPS's <alps/ngs/mpi.hpp> and ALPSCore's <alps/utilities/mpi.hpp>, which gf
// includes, both add to namespace alps::mpi. Check that one program can use both.
#include <alps/config.h>  // alps/ngs/mpi.hpp tests ALPS_HAVE_MPI without including it
#include <alps/ngs/mpi.hpp>

#include <alps/gf/gf.hpp>

#include <boost/mpi/communicator.hpp>

#include <gtest/gtest.h>

#include <complex>
#include <functional>

TEST(GfMpiCoexistence, BothMpiInterfacesWork)
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

    // ALPSCore's broadcast, directly and through a Green's function.
    int value = rank == 0 ? 42 : 0;
    alps::mpi::broadcast(comm, value, 0);
    EXPECT_EQ(42, value);

    typedef alps::gf::matsubara_mesh<alps::gf::mesh::POSITIVE_ONLY> matsubara_mesh;
    alps::gf::omega_gf gf(matsubara_mesh(10.0, 8));
    gf.initialize();
    alps::gf::matsubara_index omega;
    omega = 3;
    if (rank == 0)
        gf(omega) = std::complex<double>(3, 4);
    gf.broadcast(comm, 0);
    EXPECT_EQ(std::complex<double>(3, 4), gf(omega));
}
