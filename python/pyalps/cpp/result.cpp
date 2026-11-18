// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// pyalps.ngs.result: a result of <alps/alea/observables.hpp>, held through
// pyalps::result (see accumulator_handles.hpp).

#include <nanobind/nanobind.h>
#include "accumulator_handles.hpp"
#include "archive_savable.hpp"
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <alps/hdf5/archive.hpp>

#include <cmath>
#include <cstddef>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace nb = nanobind;
namespace alea = alps::alea;
using pyalps::result;
using pyalps::wrap;

namespace {

// Allocate a heap-owned 1-D numpy array of T and copy `n` elements from `src`.
template <typename T>
nb::ndarray<nb::numpy, T> make_1d(T const * src, std::size_t n) {
    T * buf = new T[n ? n : 1];
    if (n) std::memcpy(buf, src, n * sizeof(T));
    nb::capsule owner(buf, [](void * p) noexcept {
        delete[] static_cast<T *>(p);
    });
    std::size_t shape[1] = { n };
    return nb::ndarray<nb::numpy, T>(buf, 1, shape, owner);
}

nb::object vec_to_numpy(std::vector<double> const & v) {
    return nb::cast(make_1d<double>(v.data(), v.size()));
}

// Results with one component read back as floats, others as numpy arrays.
nb::object by_size(std::vector<double> const & values) {
    if (values.size() == 1)
        return nb::float_(values[0]);
    return vec_to_numpy(values);
}

std::string print(result const & self) {
    if (self.empty())
        return "No Measurements";
    std::ostringstream out;
    out << *self.value;
    return out.str();
}

// Functions of a result propagate its error through the jackknife over its
// bins.
result apply(result const & self, double (*op)(double)) {
    return wrap(alea::apply(op, self.get()));
}

double square(double x) { return x * x; }
double cube(double x) { return x * x * x; }

} // namespace

NB_MODULE(pyngsresult_c, m) {
    nb::module_::import_("pyalps.cxx.pyngsobservable_c");

    m.def("observable2result", [](alea::accumulator const & observable) {
        return wrap(observable.result());
    });

    using R = result;
#define PYALPS_RESULT_FUNCTION(NAME, OP) \
    .def(#NAME, [](R const & self) { return apply(self, OP); })

    nb::class_<R>(m, "result")
        .def(nb::init<>())
        // An independent copy, as copy.deepcopy() gives.
        .def("__init__", [](R * self, R const & other) {
            new (self) R{other.value ? std::make_shared<alea::observable_result>(*other.value) : nullptr};
        })
        .def("__repr__", &print)
        .def("__deepcopy__", [](R const & self, nb::handle /*memo*/) {
            return R{self.value ? std::make_shared<alea::observable_result>(*self.value) : nullptr};
        })
        .def_prop_ro("mean", [](R const & self) { return by_size(self.get().mean_vector()); })
        .def_prop_ro("error", [](R const & self) { return by_size(self.get().error_vector()); })
        // The integrated autocorrelation time; NaN where it is not known.
        .def_prop_ro("tau", [](R const & self) { return by_size(self.get().tau_vector()); })
        // The bin means of a result with one component.
        .def_prop_ro("bins", [](R const & self) {
            if (self.get().size() != 1)
                throw std::runtime_error("bins are available for scalar results");
            std::vector<double> bins;
            for (std::vector<double> const & bin : self.get().bin_means())
                bins.push_back(bin[0]);
            return vec_to_numpy(bins);
        })
        .def_prop_ro("count", [](R const & self) -> std::size_t {
            return self.value ? self.value->count() : 0;
        })

        .def("__pos__", [](R const & self) { return wrap(self.get()); })
        .def("__neg__", [](R const & self) { return wrap(-self.get()); })
        .def("__abs__", [](R const & self) { return apply(self, static_cast<double (*)(double)>(&std::fabs)); })
        .def("__pow__", [](R const & self, double exponent) {
            return wrap(alea::apply([exponent](double x) { return std::pow(x, exponent); }, self.get()));
        })

        // In-place operators modify the result and return the same Python
        // object (Python's __i*__ semantics). rv_policy::none finds the
        // existing wrapper of self.
        .def("__iadd__", [](R & self, R const & o) -> R & { self.get() = self.get() + o.get(); return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__iadd__", [](R & self, double o)    -> R & { self.get() = self.get() + o; return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__isub__", [](R & self, R const & o) -> R & { self.get() = self.get() - o.get(); return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__isub__", [](R & self, double o)    -> R & { self.get() = self.get() - o; return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__imul__", [](R & self, R const & o) -> R & { self.get() = self.get() * o.get(); return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__imul__", [](R & self, double o)    -> R & { self.get() = self.get() * o; return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__itruediv__", [](R & self, R const & o) -> R & { self.get() = self.get() / o.get(); return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__itruediv__", [](R & self, double o)    -> R & { self.get() = self.get() / o; return self; }, nb::rv_policy::none, nb::is_operator())

        // Binary operators — left and right forms. nb::is_operator()
        // marks them as Python operator overloads so mixed-type
        // failures return NotImplemented rather than raising TypeError
        // (giving Python's reflected operator machinery a chance).
        .def("__add__",  [](R const & a, R const & b) { return wrap(a.get() + b.get()); }, nb::is_operator())
        .def("__radd__", [](R const & a, R const & b) { return wrap(b.get() + a.get()); }, nb::is_operator())
        .def("__add__",  [](R const & a, double b)    { return wrap(a.get() + b); }, nb::is_operator())
        .def("__radd__", [](R const & a, double b)    { return wrap(b + a.get()); }, nb::is_operator())
        .def("__sub__",  [](R const & a, R const & b) { return wrap(a.get() - b.get()); }, nb::is_operator())
        .def("__rsub__", [](R const & a, R const & b) { return wrap(b.get() - a.get()); }, nb::is_operator())
        .def("__sub__",  [](R const & a, double b)    { return wrap(a.get() - b); }, nb::is_operator())
        .def("__rsub__", [](R const & a, double b)    { return wrap(b - a.get()); }, nb::is_operator())
        .def("__mul__",  [](R const & a, R const & b) { return wrap(a.get() * b.get()); }, nb::is_operator())
        .def("__rmul__", [](R const & a, R const & b) { return wrap(b.get() * a.get()); }, nb::is_operator())
        .def("__mul__",  [](R const & a, double b)    { return wrap(a.get() * b); }, nb::is_operator())
        .def("__rmul__", [](R const & a, double b)    { return wrap(b * a.get()); }, nb::is_operator())
        .def("__truediv__",  [](R const & a, R const & b) { return wrap(a.get() / b.get()); }, nb::is_operator())
        .def("__rtruediv__", [](R const & a, R const & b) { return wrap(b.get() / a.get()); }, nb::is_operator())
        .def("__truediv__",  [](R const & a, double b)    { return wrap(a.get() / b); }, nb::is_operator())
        .def("__rtruediv__", [](R const & a, double b)    { return wrap(b / a.get()); }, nb::is_operator())

        PYALPS_RESULT_FUNCTION(sq, &square)
        PYALPS_RESULT_FUNCTION(cb, &cube)
        PYALPS_RESULT_FUNCTION(sqrt, static_cast<double (*)(double)>(&std::sqrt))
        PYALPS_RESULT_FUNCTION(cbrt, static_cast<double (*)(double)>(&std::cbrt))
        PYALPS_RESULT_FUNCTION(log, static_cast<double (*)(double)>(&std::log))
        PYALPS_RESULT_FUNCTION(exp, static_cast<double (*)(double)>(&std::exp))
        PYALPS_RESULT_FUNCTION(sin, static_cast<double (*)(double)>(&std::sin))
        PYALPS_RESULT_FUNCTION(cos, static_cast<double (*)(double)>(&std::cos))
        PYALPS_RESULT_FUNCTION(tan, static_cast<double (*)(double)>(&std::tan))
        PYALPS_RESULT_FUNCTION(sinh, static_cast<double (*)(double)>(&std::sinh))
        PYALPS_RESULT_FUNCTION(cosh, static_cast<double (*)(double)>(&std::cosh))
        PYALPS_RESULT_FUNCTION(tanh, static_cast<double (*)(double)>(&std::tanh))

        .def("save", [](R const & self, alps::hdf5::archive & ar) {
            alea::save_result(ar, ar.get_context(), self.get());
        })
        // Loads the result at the archive's context. A failed load leaves
        // the result as it was; other references keep the old value.
        .def("load", [](R & self, alps::hdf5::archive & ar) {
            self.value = std::make_shared<alea::observable_result>(alea::load_result(ar, ar.get_context()));
        });
#undef PYALPS_RESULT_FUNCTION
    pyalps::mark_archive_savable(m.attr("result"));
}
