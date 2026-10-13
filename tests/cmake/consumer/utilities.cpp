// SPDX-License-Identifier: MIT
#include <alps/utility/copyright.hpp>
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
#include <stdexcept>
#include <vector>

int main() {
    const std::string name = "a/path with spaces";
    if (alps::hdf5_name_decode(alps::hdf5_name_encode(name)) != name
        || alps::version().empty() || alps::hostname().empty()
        || alps::cast<int>(std::string("42")) != 42)
        throw std::runtime_error("Utility library contract failed");

    const auto values = alps::read_vector<std::vector<double>>("1 2 3");
    if (values != std::vector<double>{1., 2., 3.}
        || alps::write_vector(values, ",") != "1,2,3")
        throw std::runtime_error("Vector I/O requires more than the utilities library");

    const std::string temporary = alps::temporary_filename("alps-utilities-contract");
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
