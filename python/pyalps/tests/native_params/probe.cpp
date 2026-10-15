// A separately compiled consumer catches ABI, cross-module and GIL errors.
#include <alps/hdf5/archive.hpp>
#include <alps/params.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <exception>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace nb = nanobind;

// A program defines "value" with its type and reads it, as solvers do.
template <typename T>
T defined(alps::params const & parameters) {
    alps::params copy(parameters);
    copy.define<T>("value", "the probed value");
    std::ostringstream missing;
    if (copy.has_missing(missing))
        throw std::invalid_argument(missing.str());
    return copy["value"].as<T>();
}

NB_MODULE(parameter_probe, module) {
    nb::module_::import_("pyalps.ngs");
    module.def("vector", &defined<std::vector<double>>);
    module.def("integer_vector", &defined<std::vector<long>>);
    module.def("text", &defined<std::string>);
    module.def("string_vector", &defined<std::vector<std::string>>);
    module.def("integer", &defined<int>);
    module.def("wide_integer", &defined<long>);
    module.def("flag", &defined<bool>);
    module.def("string_parameters", [] {
        alps::params p;
        p["value"] = std::vector<std::string>{"Sz"};
        return p;
    });
    module.def("clone", [](alps::params const & p) { return alps::params(p); });
    module.def("replace", [](alps::params & p) {
        p.define<std::vector<double>>("value", "the probed value");
        p["value"] = std::vector<double>{5., 6.};
    });
    module.def("empty_vectors", [] {
        alps::params p;
        p["integer"] = std::vector<int>();
        p["real"] = std::vector<double>();
        p["boolean"] = std::vector<bool>();
        p["text"] = std::vector<std::string>();
        return p;
    });
    module.def("load", [](alps::params & p, alps::hdf5::archive & ar) { p.load(ar); });
    module.def("save", [](alps::params const & p, alps::hdf5::archive & ar) { p.save(ar); });
    module.def("append_archive", [](std::string const & filename) {
        alps::hdf5::archive writer(filename, "a");
        writer["/from_native"] << 42;
    });
    module.def("threaded_vector", [](alps::params const & p) {
        std::vector<double> values;
        std::exception_ptr error;
        std::thread worker([&] {
            try { values = defined<std::vector<double>>(p); }
            catch (...) { error = std::current_exception(); }
        });
        worker.join();
        if (error) std::rethrow_exception(error);
        return values;
    }, nb::call_guard<nb::gil_scoped_release>());
    module.def("destroy_on_worker", [](alps::params const & p) {
        auto copy = std::make_unique<alps::params>(p);
        std::thread worker([copy = std::move(copy)]() mutable { copy.reset(); });
        worker.join();
    }, nb::call_guard<nb::gil_scoped_release>());
}
