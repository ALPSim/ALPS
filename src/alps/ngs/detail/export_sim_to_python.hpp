// Copyright (C) 2010-2012 by Lukas Gamper
//               2026      by the ALPS collaboration
// SPDX-License-Identifier: MIT
//
// Header-only nanobind support for downstream ALPS simulations.  Keeping this
// integration in an opt-in header prevents libalps itself from depending on
// Python or nanobind while preserving the historic public include path and
// ALPS_EXPORT_SIM_TO_PYTHON entry point.
#ifndef ALPS_NGS_DETAIL_EXPORT_SIM_TO_PYTHON_HPP
#define ALPS_NGS_DETAIL_EXPORT_SIM_TO_PYTHON_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/mcbase.hpp>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <string>

namespace alps {
namespace python {

namespace nb = nanobind;

template <typename Simulation>
class exported_simulation : public Simulation {
public:
    using parameters_type = typename Simulation::parameters_type;
    using result_names_type = typename Simulation::result_names_type;
    using results_type = typename Simulation::results_type;

    explicit exported_simulation(parameters_type const & parameters,
                                 std::size_t seed_offset = 0)
        : Simulation(parameters, seed_offset) {}

    ~exported_simulation() override = default;

    bool run_python(nb::object stop_callback) {
        return Simulation::run([stop_callback]() -> bool {
            nb::gil_scoped_acquire gil;
            return nb::cast<bool>(stop_callback());
        });
    }

    results_type collect_results_python(
        result_names_type const & names = result_names_type()) const {
        return names.empty() ? Simulation::collect_results()
                             : Simulation::collect_results(names);
    }

    alps::random01 & get_random() { return this->random; }
    parameters_type & get_parameters() { return this->parameters; }
    auto & get_measurements() {
        return this->measurements;
    }
};

template <typename Simulation>
void export_sim_to_python(nb::module_ & module, char const * name) {
    // nanobind's type registry is shared across extension modules. Import the
    // owning pyalps modules before declaring a derived simulation so mcbase,
    // params, archive, result and observable types are already registered.
    nb::module_::import_("pyalps.ngs");
    nb::module_::import_("pyalps.hdf5");

    using wrapper = exported_simulation<Simulation>;
    nb::class_<wrapper, alps::mcbase>(module, name)
        .def(nb::init<typename wrapper::parameters_type const &, std::size_t>(),
             nb::arg("parameters"), nb::arg("seed_offset") = 0)
        .def_prop_ro("random", &wrapper::get_random,
                     nb::rv_policy::reference_internal)
        .def_prop_ro("parameters", &wrapper::get_parameters,
                     nb::rv_policy::reference_internal)
        .def_prop_ro("measurements", &wrapper::get_measurements,
                     nb::rv_policy::reference_internal)
        .def("run", &wrapper::run_python, nb::arg("stop_callback"))
        .def("resultNames", &wrapper::result_names)
        .def("unsavedResultNames", &wrapper::unsaved_result_names)
        .def("collectResults", &wrapper::collect_results_python,
             nb::arg("names") = typename wrapper::result_names_type())
        .def("save",
             [](wrapper const & self, alps::hdf5::archive & archive) {
                 static_cast<Simulation const &>(self).save(archive);
             })
        .def("load",
             [](wrapper & self, alps::hdf5::archive & archive) {
                 static_cast<Simulation &>(self).load(archive);
             });
}

}  // namespace python
}  // namespace alps

#define ALPS_NANOBIND_EXPORT_SIM_TO_PYTHON(MODULE, NAME, CLASS) \
    ::alps::python::export_sim_to_python<CLASS>((MODULE), #NAME)

// Source-compatible spelling for old export.cpp files after changing their
// module declaration to ``NB_MODULE(module_name, m)``.
#define ALPS_EXPORT_SIM_TO_PYTHON(NAME, CLASS) \
    ALPS_NANOBIND_EXPORT_SIM_TO_PYTHON(m, NAME, CLASS)

#endif
