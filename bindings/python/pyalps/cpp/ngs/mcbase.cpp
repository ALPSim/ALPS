/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
 *               2026       by the ALPS collaboration                              *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the "Software"),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
// pyngsbase_c — nanobind port.
//
// Trampoline (PyMCBase) forwards the three pure-virtual mcbase methods
// (update / measure / fraction_completed) back into the Python subclass
// through nanobind's trampoline support. The old wrapper<mcbase>
// pattern becomes a standard trampoline-plus-alias pair.
//
// Params ingestion: the public alps::mcbase ctor wants an alps::params.
// libalps still declares a params(boost::python::dict) ctor in its
// header, but we don't want to drag boost::python through the
// nanobind bindings. Instead, we convert nb::dict → alps::params at the
// binding boundary by iterating and setitem-ing concrete C++ values
// (int/float/bool/str/list). That sidesteps the cross-registry issue
// and keeps the libalps ABI untouched.
#define PY_ARRAY_UNIQUE_SYMBOL pyngsbase_PyArrayHandle
#include <alps/mcbase.hpp>
#include <alps/hdf5/archive.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/make_iterator.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/trampoline.h>
namespace nb = nanobind;
#ifdef ALPS_HAVE_MPI
    #include <boost/mpi.hpp>
#endif
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>
namespace alps {
    namespace detail {
        // Convert a Python dict into an alps::params, extracting concrete
        // C++ values for each entry. This mirrors what the libalps
        // params(boost::python::dict) ctor does, but without routing the
        // nb::object through the boost::python::object variant alternative
        // — everything stays within the nanobind type registry.
        inline alps::params py_dict_to_params(nb::dict const & d) {
            alps::params p;
            for (auto item : d) {
                std::string k = nb::cast<std::string>(nb::str(item.first));
                nb::handle v = item.second;
                if (nb::isinstance<nb::bool_>(v))
                    p[k] = nb::cast<bool>(v);
                else if (nb::isinstance<nb::int_>(v))
                    p[k] = nb::cast<int>(v);
                else if (nb::isinstance<nb::float_>(v))
                    p[k] = nb::cast<double>(v);
                else if (nb::isinstance<nb::str>(v))
                    p[k] = nb::cast<std::string>(v);
                else if (nb::isinstance<nb::list>(v) || nb::isinstance<nb::tuple>(v))
                    p[k] = nb::cast<std::vector<double>>(v);
                else
                    throw nb::type_error((
                        "unsupported type for key '" + k +
                        "' in params dict (expected bool/int/float/str/list)").c_str());
            }
            return p;
        }
    }
}
namespace alps {
    // Trampoline: holds Python overrides for pure-virtuals. The
    // protected mcbase members (random / parameters / measurements)
    // are accessed via lambdas in the binding below, which friend-in
    // through PyMCBase (a protected member is visible to a derived
    // class's own member functions / friends).
    class PyMCBase : public mcbase {
        public:
            NB_TRAMPOLINE(mcbase, 3);
            #ifdef ALPS_HAVE_MPI
                PyMCBase(nb::dict const & arg,
                         std::size_t seed_offset = 42,
                         boost::mpi::communicator const & /*comm*/ = boost::mpi::communicator())
                    : mcbase(alps::detail::py_dict_to_params(arg), seed_offset)
                {}
            #else
                PyMCBase(nb::dict const & arg, std::size_t seed_offset = 42)
                    : mcbase(alps::detail::py_dict_to_params(arg), seed_offset)
                {}
            #endif
            void update() override {
                NB_OVERRIDE_PURE(update);
            }
            void measure() override {
                NB_OVERRIDE_PURE(measure);
            }
            double fraction_completed() const override {
                NB_OVERRIDE_PURE(fraction_completed);
            }
            // Accessors for protected mcbase members. Called from the
            // binding lambdas below (they friend-in through PyMCBase).
            alps::random01 & get_random() { return random; }
            mcbase::parameters_type & get_parameters() { return parameters; }
            alps::mcobservables & get_measurements() { return measurements; }
            // mcbase::run takes a std::function<bool()>; wrap a Python
            // callable so the stop_callback can be driven from Python.
            bool run_py(nb::object stop_callback) {
                return mcbase::run([stop_callback]() -> bool {
                    nb::gil_scoped_acquire gil;
                    return nb::cast<bool>(stop_callback());
                });
            }
    };
}
NB_MODULE(pyngsbase_c, m) {
    nb::class_<alps::mcbase>(m, "_mcbase", nb::never_destruct());
    nb::class_<alps::PyMCBase, alps::mcbase>(m, "mcbase")
        // Always expose the (dict, seed_offset) form from Python. When
        // ALPS_HAVE_MPI is on we'd *like* to offer an optional
        // communicator too, but boost::mpi::communicator is not a
        // nanobind-registered type so nb::arg(..).default_value() can't
        // materialise it. MPI simulations that actually need to hand
        // Python a communicator should do so from C++ using the
        // extended trampoline ctor directly.
        .def(nb::init<nb::dict const &, std::size_t>(),
             nb::arg("dict"),
             nb::arg("seed_offset") = 42)
        .def_prop_ro(
            "random",
            [](alps::PyMCBase & self) -> alps::random01 & { return self.get_random(); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "parameters",
            [](alps::PyMCBase & self) -> alps::mcbase::parameters_type & { return self.get_parameters(); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "measurements",
            [](alps::PyMCBase & self) -> alps::mcobservables & { return self.get_measurements(); },
            nb::rv_policy::reference_internal)
        .def("run",
             [](alps::PyMCBase & self, nb::object cb) { return self.run_py(std::move(cb)); })
        // Pure-virtual methods: bound on the base class; the trampoline's
        // The trampoline forwards the call into the Python subclass.
        .def("update",             &alps::mcbase::update)
        .def("measure",            &alps::mcbase::measure)
        .def("fraction_completed", &alps::mcbase::fraction_completed)
        .def("save", static_cast<void (alps::mcbase::*)(alps::hdf5::archive &) const>(
                         &alps::mcbase::save))
        .def("load", static_cast<void (alps::mcbase::*)(alps::hdf5::archive &)>(
                         &alps::mcbase::load));
}
