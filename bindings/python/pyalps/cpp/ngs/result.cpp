// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <alps/hdf5/archive.hpp>
#include <alps/ngs/cast.hpp>
#include <alps/ngs/mcresult.hpp>
#include <alps/ngs/mcobservable.hpp>
#include <alps/ngs/stacktrace.hpp>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace nb = nanobind;
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
} // namespace
namespace alps {
    namespace detail {
        template <typename T> std::string short_print_python(T const & value) {
            return cast<std::string>(value);
        }
        template <typename T> std::string short_print_python(std::vector<T> const & value) {
            switch (value.size()) {
                case 0:
                    return "[]";
                case 1:
                    return "[" + short_print_python(value.front()) + "]";
                case 2:
                    return "[" + short_print_python(value.front()) + "," + short_print_python(value.back()) + "]";
                default:
                    return "[" + short_print_python(value.front()) + ",.." + short_print_python(value.size()) + "..," + short_print_python(value.back()) + "]";
            }
        }
        inline nb::object vec_to_numpy(std::vector<double> const & v) {
            return nb::cast(make_1d<double>(v.data(), v.size()));
        }
        std::string mcresult_print(alps::mcresult const & self) {
            if (self.count() == 0)
                return "No Measurements";
            else if (self.is_type<double>())
                return short_print_python(self.mean<double>()) + "(" + short_print_python(self.count()) + ") "
                    + "+/-" + short_print_python(self.error<double>()) + " "
                    + short_print_python(self.bins<double>()) + "#" + short_print_python(self.bin_size());
            else if (self.is_type<std::vector<double> >())
                return short_print_python(self.mean<std::vector<double> >()) + "(" + short_print_python(self.count()) + ") "
                    + "+/-" + short_print_python(self.error<std::vector<double> >()) + " "
                    + short_print_python(self.bins<std::vector<double> >()) + "#" + short_print_python(self.bin_size());
            else
                throw std::runtime_error("Unsupported type." + ALPS_STACKTRACE);
        }
        nb::object mcresult_mean(alps::mcresult const & self) {
            if (self.is_type<double>())
                return nb::float_(self.mean<double>());
            else if (self.is_type<std::vector<double> >())
                return vec_to_numpy(self.mean<std::vector<double> >());
            else
                throw std::runtime_error("Unsupported type." + ALPS_STACKTRACE);
        }
        nb::object mcresult_error(alps::mcresult const & self) {
            if (self.is_type<double>())
                return nb::float_(self.error<double>());
            else if (self.is_type<std::vector<double> >())
                return vec_to_numpy(self.error<std::vector<double> >());
            else
                throw std::runtime_error("Unsupported type." + ALPS_STACKTRACE);
        }
        nb::object mcresult_tau(alps::mcresult const & self) {
            if (self.is_type<double>())
                return nb::float_(self.tau<double>());
            else if (self.is_type<std::vector<double> >())
                return vec_to_numpy(self.tau<std::vector<double> >());
            else
                throw std::runtime_error("Unsupported type." + ALPS_STACKTRACE);
        }
        nb::object mcresult_variance(alps::mcresult const & self) {
            if (self.is_type<double>())
                return nb::float_(self.variance<double>());
            else if (self.is_type<std::vector<double> >())
                return vec_to_numpy(self.variance<std::vector<double> >());
            else
                throw std::runtime_error("Unsupported type." + ALPS_STACKTRACE);
        }
        nb::object mcresult_bins(alps::mcresult const & self) {
            if (self.is_type<double>())
                return vec_to_numpy(self.bins<double>());
            else
                throw std::runtime_error("Unsupported type." + ALPS_STACKTRACE);
        }
        alps::mcresult observable2result_export(alps::mcobservable const & obs) {
            return alps::mcresult(obs);
        }
    }
}
NB_MODULE(pyngsresult_c, m) {
    using namespace alps;
    using R = alps::mcresult;
    m.def("observable2result", &alps::detail::observable2result_export);
    nb::class_<R>(m, "result")
        .def(nb::init<>())
        .def(nb::init<R const &>())
        .def("__repr__", &alps::detail::mcresult_print)
        .def("__deepcopy__",
             [](R const & self, nb::handle /*memo*/) {
                 return R(self);
             })
        .def("__abs__", static_cast<R(*)(R)>(&abs))
        .def("__pow__", static_cast<R(*)(R, double)>(&pow))
        .def_prop_ro("mean",     &alps::detail::mcresult_mean)
        .def_prop_ro("error",    &alps::detail::mcresult_error)
        .def_prop_ro("tau",      &alps::detail::mcresult_tau)
        .def_prop_ro("variance", &alps::detail::mcresult_variance)
        .def_prop_ro("bins",     &alps::detail::mcresult_bins)
        .def_prop_ro("count",    &R::count)
        // mcresult's unary +/- operate on non-const self and return a
        // reference (not a new value). Wrap them in lambdas that
        // return a fresh copy, which is what Python's +obj/-obj expect.
        .def("__pos__", [](R self) { return +self; })
        .def("__neg__", [](R self) { return -self; })
        // In-place operators — return self by reference so the original
        // object is modified in place (Python's __i*__ semantics).
        .def("__iadd__", [](R & self, R const & o) -> R & { return self += o; }, nb::is_operator())
        .def("__iadd__", [](R & self, double o)    -> R & { return self += o; }, nb::is_operator())
        .def("__isub__", [](R & self, R const & o) -> R & { return self -= o; }, nb::is_operator())
        .def("__isub__", [](R & self, double o)    -> R & { return self -= o; }, nb::is_operator())
        .def("__imul__", [](R & self, R const & o) -> R & { return self *= o; }, nb::is_operator())
        .def("__imul__", [](R & self, double o)    -> R & { return self *= o; }, nb::is_operator())
        .def("__itruediv__", [](R & self, R const & o) -> R & { return self /= o; }, nb::is_operator())
        .def("__itruediv__", [](R & self, double o)    -> R & { return self /= o; }, nb::is_operator())
        // Binary operators — left and right forms. nb::is_operator()
        // marks them as Python operator overloads so mixed-type
        // failures return NotImplemented rather than raising TypeError
        // (giving Python's reflected operator machinery a chance).
        .def("__add__",  [](R const & a, R const & b) { return a + b; }, nb::is_operator())
        .def("__radd__", [](R const & a, R const & b) { return b + a; }, nb::is_operator())
        .def("__add__",  [](R const & a, double b)    { return a + b; }, nb::is_operator())
        .def("__radd__", [](R const & a, double b)    { return b + a; }, nb::is_operator())
        .def("__sub__",  [](R const & a, R const & b) { return a - b; }, nb::is_operator())
        .def("__rsub__", [](R const & a, R const & b) { return b - a; }, nb::is_operator())
        .def("__sub__",  [](R const & a, double b)    { return a - b; }, nb::is_operator())
        .def("__rsub__", [](R const & a, double b)    { return b - a; }, nb::is_operator())
        .def("__mul__",  [](R const & a, R const & b) { return a * b; }, nb::is_operator())
        .def("__rmul__", [](R const & a, R const & b) { return b * a; }, nb::is_operator())
        .def("__mul__",  [](R const & a, double b)    { return a * b; }, nb::is_operator())
        .def("__rmul__", [](R const & a, double b)    { return b * a; }, nb::is_operator())
        .def("__truediv__",  [](R const & a, R const & b) { return a / b; }, nb::is_operator())
        .def("__rtruediv__", [](R const & a, R const & b) { return b / a; }, nb::is_operator())
        .def("__truediv__",  [](R const & a, double b)    { return a / b; }, nb::is_operator())
        .def("__rtruediv__", [](R const & a, double b)    { return b / a; }, nb::is_operator())
        .def("sq",   static_cast<R(*)(R)>(&sq))
        .def("cb",   static_cast<R(*)(R)>(&cb))
        .def("sqrt", static_cast<R(*)(R)>(&sqrt))
        .def("cbrt", static_cast<R(*)(R)>(&cbrt))
        .def("exp",  static_cast<R(*)(R)>(&exp))
        .def("log",  static_cast<R(*)(R)>(&log))
        .def("sin",  static_cast<R(*)(R)>(&sin))
        .def("cos",  static_cast<R(*)(R)>(&cos))
        .def("tan",  static_cast<R(*)(R)>(&tan))
        .def("sinh", static_cast<R(*)(R)>(&sinh))
        .def("cosh", static_cast<R(*)(R)>(&cosh))
        .def("tanh", static_cast<R(*)(R)>(&tanh))
        .def("save", &R::save)
        .def("load", &R::load);
}
