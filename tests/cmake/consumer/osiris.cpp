// SPDX-License-Identifier: MIT
#include <alps/osiris/comm.h>
#include <alps/osiris/buffer.h>
#include <stdexcept>

#ifdef ALPS_HAVE_MPI
#include <mpi.h>
#endif

int main(int argc, char** argv) {
    if (alps::runs_parallel()) throw std::runtime_error("Initial serial state");
    alps::comm_init(argc, argv, false);
    if (!alps::is_master() || alps::detail::local_id() != 0
        || alps::detail::invalid_id() != -1 || !alps::local_process().local()
        || alps::local_process() != alps::master_process()
        || alps::all_processes() != alps::ProcessList{alps::Process(0)})
        throw std::runtime_error("Serial communication exports");

    alps::detail::Buffer buffer;
    buffer.write(42);
    buffer.write(-1.25);
    int integer = 0;
    double real = 0;
    buffer.read(integer);
    buffer.read(real);
    if (integer != 42 || real != -1.25) throw std::runtime_error("Buffer round trip");

#ifdef ALPS_HAVE_MPI
    MPI_Init(&argc, &argv);
    if (alps::runs_parallel()) throw std::runtime_error("External MPI initialization");
    alps::comm_init(argc, argv, true);
    int rank = 0, size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (!alps::runs_parallel() || alps::detail::local_id() != rank
        || alps::is_master() != (rank == 0) || !alps::local_process().valid()
        || alps::all_processes().size() != static_cast<std::size_t>(size))
        throw std::runtime_error("MPI communication exports");
#else
    if (!alps::local_process().valid() || alps::Process().valid())
        throw std::runtime_error("Serial process validity");
    bool rejected = false;
    try { alps::comm_init(argc, argv, true); }
    catch (std::runtime_error const&) { rejected = true; }
    if (!rejected) throw std::runtime_error("Serial build accepted MPI initialization");
#endif
    alps::comm_exit();
}
