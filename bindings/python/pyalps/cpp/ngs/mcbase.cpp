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
// through nanobind's trampoline support. Binding mcbase itself with
// PyMCBase as its alias preserves the public base class of downstream
// simulations, as in the Boost.Python bindings.
//
// Params ingestion: the public alps::mcbase ctor wants an alps::params.
// We convert nb::dict → alps::params at the binding boundary through
// the shared ladder in ../dict_to_params.hpp, so mcbase, params and
// the application modules ingest parameters identically.
#define PY_ARRAY_UNIQUE_SYMBOL pyngsbase_PyArrayHandle
#include <alps/mcbase.hpp>
#include <alps/hdf5/archive.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/make_iterator.h>
#include <nanobind/stl/function.h>
#include <nanobind/trampoline.h>
namespace nb = nanobind;
#include <cstddef>
#include <functional>
#include <type_traits>
#include <utility>
#include "../dict_to_params.hpp"
namespace alps {
    static_assert(std::has_virtual_destructor<mcbase>::value,
                  "mcbase must safely destroy nanobind trampoline aliases");

    // Trampoline: holds Python overrides for pure-virtuals. The
    // protected mcbase members (random / parameters / measurements)
    // are accessed via lambdas in the binding below, which friend-in
    // through PyMCBase (a protected member is visible to a derived
    // class's own member functions / friends).
    class PyMCBase : public mcbase {
        public:
            // Slot count = the number of NB_OVERRIDE* calls below.
            // mcbase (src/alps/mcbase.hpp) declares five virtuals:
            // update / measure / fraction_completed (pure) and
            // save(archive&) / load(archive&); all five must be
            // forwarded so Python overrides are seen by C++ callers.
            NB_TRAMPOLINE(mcbase, 5);
            PyMCBase(nb::dict const & arg,
                     std::size_t seed_offset = 42,
                     nb::handle /*communicator*/ = nb::none())
                : mcbase(pyalps::params_from_dict(arg), seed_offset)
            {}
            void update() override {
                NB_OVERRIDE_PURE(update);
            }
            void measure() override {
                NB_OVERRIDE_PURE(measure);
            }
            double fraction_completed() const override {
                NB_OVERRIDE_PURE(fraction_completed);
            }
            // Non-pure: fall through to the C++ implementation when the
            // Python subclass doesn't override (NB_OVERRIDE, not _PURE).
            void save(alps::hdf5::archive & ar) const override {
                NB_OVERRIDE(save, ar);
            }
            void load(alps::hdf5::archive & ar) override {
                NB_OVERRIDE(load, ar);
            }
            // Accessors for protected mcbase members. Called from the
            // binding lambdas below (they friend-in through PyMCBase).
            alps::random01 & get_random() { return random; }
            mcbase::parameters_type & get_parameters() { return parameters; }
            alps::mcobservables & get_measurements() { return measurements; }
    };
}
NB_MODULE(pyngsbase_c, m) {
    nb::class_<alps::mcbase, alps::PyMCBase>(m, "mcbase")
        // Retain the legacy third argument without binding Boost.MPI. The
        // Boost.Python-era constructor accepted a communicator but never
        // passed it to alps::mcbase (which has no communicator constructor),
        // so accepting and ignoring it is behaviorally faithful. Python-side
        // communication is provided by pyalps.mpi's mpi4py adapter.
        .def(nb::init<nb::dict const &, std::size_t, nb::handle>(),
             nb::arg("dict"),
             nb::arg("seed_offset") = 42,
             nb::arg("communicator") = nb::none())
        .def_prop_ro(
            "random",
            [](alps::mcbase & self) -> alps::random01 & {
                return dynamic_cast<alps::PyMCBase &>(self).get_random();
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "parameters",
            [](alps::mcbase & self) -> alps::mcbase::parameters_type & {
                return dynamic_cast<alps::PyMCBase &>(self).get_parameters();
            },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "measurements",
            [](alps::mcbase & self) -> alps::mcobservables & {
                return dynamic_cast<alps::PyMCBase &>(self).get_measurements();
            },
            nb::rv_policy::reference_internal)
        .def("run",
             [](alps::mcbase & self, nb::object stop_callback) {
                 return self.run([stop_callback = std::move(stop_callback)]() -> bool {
                     nb::gil_scoped_acquire gil;
                     return nb::cast<bool>(stop_callback());
                 });
             })
        // Pure-virtual methods are bound on the base class; the trampoline
        // forwards each call into the Python subclass.
        .def("update",             &alps::mcbase::update)
        .def("measure",            &alps::mcbase::measure)
        .def("fraction_completed", &alps::mcbase::fraction_completed)
        .def("save", static_cast<void (alps::mcbase::*)(alps::hdf5::archive &) const>(
                         &alps::mcbase::save))
        .def("load", static_cast<void (alps::mcbase::*)(alps::hdf5::archive &)>(
                         &alps::mcbase::load));
}
