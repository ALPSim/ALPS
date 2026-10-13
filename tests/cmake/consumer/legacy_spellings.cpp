// SPDX-License-Identifier: MIT
// A consumer written against ALPS 3.0.0: old include spellings and alps::ngs
// names, which the installed legacy layer forwards to the canonical ones.
#include <alps/ngs/cast.hpp>
#include <alps/ngs/short_print.hpp>
#include <alps/ngs/signal.hpp>
#include <alps/ngs/stacktrace.hpp>
#include <alps/utility/cli.hpp>
#include <alps/utility/copyright.hpp>
#include <alps/utility/encode.hpp>
#include <alps/utility/numeric_cast.hpp>
#include <alps/utility/os.hpp>
#include <alps/utility/temporary_filename.hpp>
#include <alps/utility/vectorio.hpp>
#include <complex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#if defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

int main() {
    static_assert(std::is_same<alps::ngs::signal, alps::signal>::value, "alps::ngs::signal");
    (void)alps::ngs::stacktrace();
    std::ostringstream out;
    out << alps::short_print(std::vector<int>{1, 2}) << alps::cast<std::string>(7);
    std::vector<int> parsed;
    alps::read_vector_resize(std::string("1 2"), parsed);
    if (out.str() != "[1,2]7" || parsed.size() != 2
        || alps::hdf5_name_decode(alps::hdf5_name_encode("a/b")) != "a/b"
        || alps::version().empty() || alps::temporary_filename("legacy").empty()
        || alps::temp_directory_path().empty() || !alps::cli_is_master()
        || alps::numeric_cast<double>(std::complex<double>(1.5, 0.0)) != 1.5)
        throw std::runtime_error("Legacy spellings no longer reach the SDK");
}
