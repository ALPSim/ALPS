// SPDX-License-Identifier: MIT
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/ngs/signal.hpp>
#include <boost/filesystem/operations.hpp>
#include <complex>
#include <stdexcept>
#include <vector>

int main() {
    const char* filename = "hdf5-component-contract.h5";
    const std::vector<std::complex<double>> expected{{1., 2.}, {3., -4.}};
    {
        alps::hdf5::archive writer(filename, "w");
        writer["/values"] << expected;
        writer["/values/@description"] << std::string("complex vector");
        // Multiple handles must share one archive context, including typed I/O.
        alps::hdf5::archive reader(filename, "r");
        std::vector<std::complex<double>> actual;
        std::string description;
        reader["/values"] >> actual;
        reader["/values/@description"] >> description;
        if (actual != expected || description != "complex vector")
            throw std::runtime_error("HDF5 component round trip failed");
        bool caught = false;
        try {
            int missing = 0;
            reader["/missing"] >> missing;
        } catch (const alps::hdf5::path_not_found&) {
            caught = true;
        }
        if (!caught)
            throw std::runtime_error("HDF5 component exception identity failed");
    }
    // Check exported signal state without delivering process-level signals.
    alps::ngs::signal signals;
    alps::ngs::signal::slot(17);
    if (signals.empty() || signals.top() != 17)
        throw std::runtime_error("HDF5 signal state is not shared");
    signals.pop();
    if (!signals.empty())
        throw std::runtime_error("HDF5 signal queue was not cleared");
    boost::filesystem::remove(filename);
}
