// SPDX-License-Identifier: MIT
#include <alps/osiris/comm.h>
#include <alps/osiris/buffer.h>

#include <iostream>
#include <stdexcept>

#ifdef ALPS_HAVE_MPI
#include <mpi.h>
#endif

namespace {
void require(bool condition, char const* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main(int argc, char** argv) {
    try {
        // Serial queries must work before either Osiris or MPI initialization.
        require(!alps::runs_parallel(), "Osiris starts in serial mode");
        alps::comm_init(argc, argv, false);
        require(alps::is_master(), "Serial process is master");
        require(alps::detail::local_id() == 0, "Serial process ID");
        require(alps::detail::invalid_id() == -1, "Invalid process ID");
        require(alps::local_process().local(), "Local process descriptor");
        require(alps::local_process() == alps::master_process(), "Serial master descriptor");
        require(alps::all_processes() == alps::ProcessList{alps::Process(0)},
                "Serial process enumeration");

        alps::detail::Buffer buffer;
        buffer.write(42);
        buffer.write(-1.25);
        int integer = 0;
        double real = 0;
        buffer.read(integer);
        buffer.read(real);
        require(integer == 42 && real == -1.25, "Buffer round trip");

#ifdef ALPS_HAVE_MPI
        // Exercise the shared state when MPI was initialized by the embedding app.
        MPI_Init(&argc, &argv);
        require(!alps::runs_parallel(), "External MPI initialization does not opt Osiris in");
        alps::comm_init(argc, argv, true);
        int rank = 0, size = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        require(alps::runs_parallel(), "Osiris adopts initialized MPI");
        require(alps::detail::local_id() == rank, "MPI process ID");
        require(alps::is_master() == (rank == 0), "MPI master query");
        require(alps::local_process().valid(), "MPI local process validity");
        require(alps::all_processes().size() == static_cast<std::size_t>(size),
                "MPI process enumeration");
        alps::comm_exit();
#else
        require(alps::local_process().valid(), "Serial process validity");
        require(!alps::Process().valid(), "Invalid serial process descriptor");
        bool rejected = false;
        try { alps::comm_init(argc, argv, true); }
        catch (std::runtime_error const&) { rejected = true; }
        require(rejected, "Serial build rejects MPI initialization");
        alps::comm_exit();
#endif
    } catch (std::exception const& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
