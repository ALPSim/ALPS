// A separately compiled consumer catches ABI, cross-module and GIL errors.
#include <alps/ngs/params.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/complex.h>
#include <exception>
#include <thread>

namespace nb = nanobind;
NB_MODULE(parameter_probe, module) {
    nb::module_::import_("pyalps.ngs");
    module.def("vector", [](alps::params const & p) { return p["value"].cast<std::vector<double>>(); });
    module.def("complex_vector", [](alps::params const & p) { return p["value"].cast<std::vector<std::complex<double>>>(); });
    module.def("integer_vector", [](alps::params const & p) { return p["value"].cast<std::vector<long long>>(); });
    module.def("integer", [](alps::params const & p) { return p["value"].cast<int>(); });
    module.def("wide_integer", [](alps::params const & p) { return p["value"].cast<long long>(); });
    module.def("clone", [](alps::params const & p) { return alps::params(p); });
    module.def("replace", [](alps::params & p) { p["value"] = std::vector<double>{5., 6.}; });
    module.def("empty_vectors", [] {
        alps::params p;
        p["integer"] = std::vector<int>();
        p["real"] = std::vector<double>();
        p["complex"] = std::vector<std::complex<double>>();
        p["boolean"] = std::vector<bool>();
        p["text"] = std::vector<std::string>();
        return p;
    });
    module.def("load", [](alps::params & p, alps::hdf5::archive & ar) { p.load(ar); });
    module.def("threaded_vector", [](alps::params const & p) {
        std::vector<double> values;
        std::exception_ptr error;
        std::thread worker([&] {
            try { values = p["value"].cast<std::vector<double>>(); }
            catch (...) { error = std::current_exception(); }
        });
        worker.join();
        if (error) std::rethrow_exception(error);
        return values;
    }, nb::call_guard<nb::gil_scoped_release>());
    module.def("destroy_on_worker", [](alps::params & p) {
        auto copy = std::make_unique<alps::params>(p);
        p.erase("value");
        std::thread worker([copy = std::move(copy)]() mutable { copy.reset(); });
        worker.join();
    }, nb::call_guard<nb::gil_scoped_release>());
}
