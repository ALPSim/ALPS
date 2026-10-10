// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>
#include <alps/osiris/comm.h>
#include <alps/osiris/buffer.h>

#include <iostream>
#include <stdexcept>

#ifdef ALPS_HAVE_MPI
#include <mpi.h>
#endif


TEST(OsirisCommunication, InitializationAndBufferContracts) {
        int argc = 0;
        char** argv = nullptr;
        // Serial queries must work before either Osiris or MPI initialization.
        ASSERT_TRUE((!alps::runs_parallel())) << "Osiris starts in serial mode";
        alps::comm_init(argc, argv, false);
        ASSERT_TRUE((alps::is_master())) << "Serial process is master";
        ASSERT_TRUE((alps::detail::local_id() == 0)) << "Serial process ID";
        ASSERT_TRUE((alps::detail::invalid_id() == -1)) << "Invalid process ID";
        ASSERT_TRUE((alps::local_process().local())) << "Local process descriptor";
        ASSERT_TRUE((alps::local_process() == alps::master_process())) << "Serial master descriptor";
        ASSERT_TRUE((alps::all_processes() == alps::ProcessList{alps::Process(0)})) << "Serial process enumeration";

        alps::detail::Buffer buffer;
        buffer.write(42);
        buffer.write(-1.25);
        int integer = 0;
        double real = 0;
        buffer.read(integer);
        buffer.read(real);
        ASSERT_TRUE((integer == 42 && real == -1.25)) << "Buffer round trip";

#ifdef ALPS_HAVE_MPI
        // Exercise the shared state when MPI was initialized by the embedding app.
        MPI_Init(&argc, &argv);
        ASSERT_TRUE((!alps::runs_parallel())) << "External MPI initialization does not opt Osiris in";
        alps::comm_init(argc, argv, true);
        int rank = 0, size = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        ASSERT_TRUE((alps::runs_parallel())) << "Osiris adopts initialized MPI";
        ASSERT_TRUE((alps::detail::local_id() == rank)) << "MPI process ID";
        ASSERT_TRUE((alps::is_master() == (rank == 0))) << "MPI master query";
        ASSERT_TRUE((alps::local_process().valid())) << "MPI local process validity";
        ASSERT_TRUE((alps::all_processes().size() == static_cast<std::size_t>(size))) << "MPI process enumeration";
        alps::comm_exit();
#else
        ASSERT_TRUE((alps::local_process().valid())) << "Serial process validity";
        ASSERT_TRUE((!alps::Process().valid())) << "Invalid serial process descriptor";
        bool rejected = false;
        try { alps::comm_init(argc, argv, true); }
        catch (std::runtime_error const&) { rejected = true; }
        ASSERT_TRUE((rejected)) << "Serial build rejects MPI initialization";
        alps::comm_exit();
#endif
}
