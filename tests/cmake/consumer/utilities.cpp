// SPDX-License-Identifier: MIT
#include <alps/utility/copyright.hpp>
#include <alps/utility/citations.hpp>
#include <alps/utility/encode.hpp>
#include <alps/utility/os.hpp>
#include <alps/utility/temporary_filename.hpp>
#include <alps/utility/vectorio.hpp>
#include <alps/utility/vmusage.hpp>
#include <alps/ngs/cast.hpp>
#include <alps/ngs/short_print.hpp>
#include <alps/ngs/sleep.hpp>
#include <boost/filesystem/operations.hpp>
#include <sstream>
#include <functional>
#include <stdexcept>
#include <vector>
#ifndef _WIN32
#include <sys/stat.h>
#endif
#ifdef ALPS_HAVE_MPI
#include <alps/ngs/boost_mpi.hpp>
#endif

int main(int argc, char** argv) {
#ifdef ALPS_HAVE_MPI
    boost::mpi::environment environment(argc, argv);
    boost::mpi::communicator world;
    std::vector<int> local{1, 2}, total;
    boost::mpi::reduce(world, local, total, std::plus<int>(), 0);
    if (world.rank() == 0 && total != std::vector<int>{world.size(), 2 * world.size()})
        throw std::runtime_error("Utility MPI reduction failed");
#endif
    std::ostringstream banner;
    alps::print_copyright(banner, "hybridization");
    if (alps::citation_text("hybridization").empty()
        || alps::citation_details("interaction").find("10.1103/PhysRevB.72.035122") == std::string::npos
        || banner.str().find("copyright (c)") == std::string::npos)
        throw std::runtime_error("Utilities citation exports failed");
    const std::string name = "a/path with spaces";
    if (alps::hdf5_name_decode(alps::hdf5_name_encode(name)) != name
        || alps::version().empty() || alps::hostname().empty()
        || alps::cast<int>(std::string("42")) != 42)
        throw std::runtime_error("Utility library contract failed");

    const auto values = alps::read_vector<std::vector<double>>("1 2 3");
    if (values != std::vector<double>{1., 2., 3.}
        || alps::write_vector(values, ",") != "1,2,3")
        throw std::runtime_error("Vector I/O requires more than the utilities library");

#ifndef _WIN32
    const auto previous_mask = umask(0);
#endif
    const std::string temporary = alps::temporary_filename("alps-utilities-contract");
#ifndef _WIN32
    umask(previous_mask);
    struct stat status;
    if (stat(temporary.c_str(), &status) != 0 || (status.st_mode & 0777) != 0600)
        throw std::runtime_error("Utility temporary file is not private to its owner");
#endif
    if (!boost::filesystem::remove(temporary))
        throw std::runtime_error("Utility temporary file was not created");
    alps::sleep(0);
    (void)alps::vmusage();

    // Call the out-of-line overloads to check DLL exports as well as templates.
    float f = 1.2345f;
    double d = 1.2345;
    long double l = 1.2345L;
    std::ostringstream output;
    output << alps::detail::short_print_proxy<float>(f, 3) << ' '
           << alps::detail::short_print_proxy<double>(d, 3) << ' '
           << alps::detail::short_print_proxy<long double>(l, 3);
    if (output.str() != "1.23 1.23 1.23")
        throw std::runtime_error("Utility formatting exports failed");
}
