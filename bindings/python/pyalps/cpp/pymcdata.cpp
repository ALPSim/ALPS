// Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
//               Lukas Gamper <gamperl@gmail.com>,
//               Matthias Troyer <troyer@itp.phys.ethz.ch>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <alps/alea/mcdata.hpp>
#include "numpy_compat.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
namespace nb = nanobind;
namespace alps {
    namespace python {
        // Build a 1-D numpy.ndarray (dtype=float64) from a
        // std::vector<double>. Allocates a fresh array via
        // numpy.empty + memcpy through the buffer protocol — no
        // numpy headers.
        inline nb::object vec_to_numpy(std::vector<double> const & v) {
            return make_numpy_array<double>(v.data(), {v.size()});
        }
        // Build a 2-D numpy.ndarray from a vector<vector<double>>.
        // Rows must be equal-length; on mismatch throw a value error
        // (mcdata doesn't produce ragged bins/jackknife tables).
        inline nb::object matrix_to_numpy(std::vector<std::vector<double>> const & m) {
            std::size_t rows = m.size();
            std::size_t cols = rows > 0 ? m.front().size() : 0;
            for (auto const & row : m) {
                if (row.size() != cols)
                    throw std::runtime_error("mcdata returned ragged 2-D table; refusing to convert to numpy");
            }
            std::vector<double> flat(rows * cols);
            double * dst = flat.data();
            for (auto const & row : m) {
                std::copy(row.begin(), row.end(), dst);
                dst += cols;
            }
            return make_numpy_array<double>(flat.data(), {rows, cols});
        }
        // __repr__ for mcdata<double> — "<mean> +/- <error>".
        template <typename T>
        std::string print_mcdata_scalar(alps::alea::mcdata<T> const & self) {
            std::ostringstream ss;
            ss << self.mean() << " +/- " << self.error();
            return ss.str();
        }
        // __repr__ for mcdata<vector<T>> — newline-joined scalar reprs.
        // The const_iterator's operator-> returns through boost::addressof
        // on a rvalue scalar mcdata<T>, which triggers the deleted
        // overload; take a local copy instead.
        template <typename T>
        std::string print_mcdata_vector(alps::alea::mcdata<std::vector<T>> const & self) {
            std::ostringstream ss;
            bool first = true;
            for (auto it = self.begin(); it != self.end(); ++it) {
                if (!first) ss << "\n";
                first = false;
                alps::alea::mcdata<T> entry = *it;
                ss << entry.mean() << " +/- " << entry.error();
            }
            return ss.str();
        }
        // __format__ for mcdata<T> — defers to builtins.format on
        // mean and error separately, then joins with " +/- ".
        template <typename T>
        std::string format_mcdata_scalar(alps::alea::mcdata<T> const & self,
                                         std::string const & spec) {
            nb::object fmt = nb::module_::import_("builtins").attr("format");
            std::string m = nb::cast<std::string>(fmt(self.mean(), spec));
            std::string e = nb::cast<std::string>(fmt(self.error(), spec));
            return m + " +/- " + e;
        }
        template <typename T>
        std::string format_mcdata_vector(alps::alea::mcdata<std::vector<T>> const & self,
                                         std::string const & spec) {
            std::ostringstream ss;
            bool first = true;
            for (auto it = self.begin(); it != self.end(); ++it) {
                if (!first) ss << "\n";
                first = false;
                alps::alea::mcdata<T> entry = *it;
                ss << format_mcdata_scalar(entry, spec);
            }
            return ss.str();
        }
        // Indexing: mcdata<vector<T>>[i] returns a scalar mcdata<T>
        // view; mcdata<vector<T>>[slice] returns a sliced
        // mcdata<vector<T>>.
        template <typename T>
        nb::object mcdata_vector_getitem(alps::alea::mcdata<std::vector<T>> & data,
                                         nb::object const & key) {
            std::size_t n = data.mean().size();
            if (nb::isinstance<nb::slice>(key)) {
                nb::slice s = nb::borrow<nb::slice>(key);
                auto [start, stop, step, slicelength] = s.compute(n);
                if (step != 1)
                    throw nb::index_error("slice step size not supported.");
                if (start > stop)
                    return nb::cast(alps::alea::mcdata<std::vector<T>>());
                return nb::cast(alps::alea::mcdata<std::vector<T>>(
                    data, start, stop));
            }
            long index = nb::cast<long>(key);
            if (index < 0)
                index += static_cast<long>(n);
            if (index < 0 || static_cast<std::size_t>(index) >= n)
                throw nb::index_error("Index out of range");
            return nb::cast(alps::alea::mcdata<T>(
                data, static_cast<std::size_t>(index)));
        }
        template <typename T>
        bool mcdata_vector_contains(alps::alea::mcdata<std::vector<T>> & data,
                                    nb::object const & key) {
            // Best-effort: accept either a scalar mcdata<T> or a value
            // that casts cleanly; anything else is "not in".
            // mcdata<std::vector<T>>::const_iterator does not fully
            // model std::ranges::input_range (the iterator traits
            // needed for borrowed_iterator_t aren't wired up), so the
            // ranges-form of std::find is ill-formed here. Keep the
            // classical iterator pair.
            try {
                auto probe = nb::cast<alps::alea::mcdata<T>>(key);
                return std::find(data.begin(), data.end(), probe) != data.end();
            } catch (nb::cast_error const &) {
            }
            return false;
        }
    }
}
NB_MODULE(pymcdata_c, m) {
    using alps::alea::mcdata;
    namespace pymod = alps::python;
    using Scalar = mcdata<double>;
    using Vector = mcdata<std::vector<double>>;
    // mcdata's transcendentals (sq, cb, sqrt, cbrt, exp, log, sin, …)
    // live in namespace alps::alea and are found via ADL when we pass
    // an mcdata<T> argument. Capture each as a lambda so the binding
    // doesn't have to cast through overloaded name lookup, which
    // fights <cmath>'s own abs/pow/sqrt etc. sitting in global scope.
    nb::class_<Scalar>(m, "MCScalarData",
        "Scalar Monte Carlo data. Supports +, -, *, /, +=, -=, *=, /=, "
        "abs, pow, sq, cb, sqrt, cbrt, exp, log, sin, cos, tan, sinh, cosh, tanh.")
        .def(nb::init<>())
        .def(nb::init<double>(), nb::arg("mean"))
        .def(nb::init<double, double>(), nb::arg("mean"), nb::arg("error"))
        .def_prop_ro("mean",     [](Scalar const & v) { return v.mean(); })
        .def_prop_ro("error",    [](Scalar const & v) { return v.error(); })
        .def_prop_ro("tau",      [](Scalar const & v) { return v.tau(); })
        .def_prop_ro("variance", [](Scalar const & v) { return v.variance(); })
        .def_prop_ro("bins",     [](Scalar const & v) {
                                      return pymod::vec_to_numpy(v.bins());
                                  })
        .def_prop_ro("jackknife",[](Scalar const & v) {
                                      return pymod::vec_to_numpy(v.jackknife());
                                  })
        .def_prop_ro("count",    &Scalar::count)
        .def("__repr__",   &pymod::print_mcdata_scalar<double>)
        .def("__format__", &pymod::format_mcdata_scalar<double>,
             nb::arg("format_spec"))
        .def("__deepcopy__",
             [](Scalar const & self, nb::handle /*memo*/) {
                 return Scalar(self);
             })
        .def("__abs__", [](Scalar x) { using alps::alea::abs; return abs(std::move(x)); })
        .def("__pow__", [](Scalar x, double e) { using alps::alea::pow; return pow(std::move(x), e); })
        // Unary - / + on mcdata produce new values; wrap manually
        // because the library's operator+()/-() signatures aren't
        // const-returning (which is what nb::self expects).
        .def("__pos__", [](Scalar self) { return +self; })
        .def("__neg__", [](Scalar self) { return -self; })
        // In-place operators — modify self in place, return reference.
        .def("__iadd__", [](Scalar & s, Scalar const & o) -> Scalar & { return s += o; }, nb::is_operator())
        .def("__iadd__", [](Scalar & s, double o)         -> Scalar & { return s += o; }, nb::is_operator())
        .def("__isub__", [](Scalar & s, Scalar const & o) -> Scalar & { return s -= o; }, nb::is_operator())
        .def("__isub__", [](Scalar & s, double o)         -> Scalar & { return s -= o; }, nb::is_operator())
        .def("__imul__", [](Scalar & s, Scalar const & o) -> Scalar & { return s *= o; }, nb::is_operator())
        .def("__imul__", [](Scalar & s, double o)         -> Scalar & { return s *= o; }, nb::is_operator())
        .def("__itruediv__", [](Scalar & s, Scalar const & o) -> Scalar & { return s /= o; }, nb::is_operator())
        .def("__itruediv__", [](Scalar & s, double o)         -> Scalar & { return s /= o; }, nb::is_operator())
        // Binary +/-/*//: forward and reflected forms. nb::is_operator()
        // marks them so mixed-type failures return NotImplemented and
        // Python's reflected operator machinery gets a turn.
        .def("__add__",  [](Scalar const & a, Scalar const & b) { return a + b; }, nb::is_operator())
        .def("__radd__", [](Scalar const & a, Scalar const & b) { return b + a; }, nb::is_operator())
        .def("__add__",  [](Scalar const & a, double b)         { return a + b; }, nb::is_operator())
        .def("__radd__", [](Scalar const & a, double b)         { return b + a; }, nb::is_operator())
        .def("__sub__",  [](Scalar const & a, Scalar const & b) { return a - b; }, nb::is_operator())
        .def("__rsub__", [](Scalar const & a, Scalar const & b) { return b - a; }, nb::is_operator())
        .def("__sub__",  [](Scalar const & a, double b)         { return a - b; }, nb::is_operator())
        .def("__rsub__", [](Scalar const & a, double b)         { return b - a; }, nb::is_operator())
        .def("__mul__",  [](Scalar const & a, Scalar const & b) { return a * b; }, nb::is_operator())
        .def("__rmul__", [](Scalar const & a, Scalar const & b) { return b * a; }, nb::is_operator())
        .def("__mul__",  [](Scalar const & a, Vector const & b) { return a * b; }, nb::is_operator())
        .def("__rmul__", [](Scalar const & a, Vector const & b) { return b * a; }, nb::is_operator())
        .def("__mul__",  [](Scalar const & a, double b)         { return a * b; }, nb::is_operator())
        .def("__rmul__", [](Scalar const & a, double b)         { return b * a; }, nb::is_operator())
        .def("__truediv__",  [](Scalar const & a, Scalar const & b) { return a / b; }, nb::is_operator())
        .def("__rtruediv__", [](Scalar const & a, Scalar const & b) { return b / a; }, nb::is_operator())
        .def("__truediv__",  [](Scalar const & a, double b)         { return a / b; }, nb::is_operator())
        .def("__rtruediv__", [](Scalar const & a, double b)         { return b / a; }, nb::is_operator())
        .def("sq",   [](Scalar x) { using alps::alea::sq;   return sq(std::move(x));   })
        .def("cb",   [](Scalar x) { using alps::alea::cb;   return cb(std::move(x));   })
        .def("sqrt", [](Scalar x) { using alps::alea::sqrt; return sqrt(std::move(x)); })
        .def("cbrt", [](Scalar x) { using alps::alea::cbrt; return cbrt(std::move(x)); })
        .def("exp",  [](Scalar x) { using alps::alea::exp;  return exp(std::move(x));  })
        .def("log",  [](Scalar x) { using alps::alea::log;  return log(std::move(x));  })
        .def("sin",  [](Scalar x) { using alps::alea::sin;  return sin(std::move(x));  })
        .def("cos",  [](Scalar x) { using alps::alea::cos;  return cos(std::move(x));  })
        .def("tan",  [](Scalar x) { using alps::alea::tan;  return tan(std::move(x));  })
        .def("sinh", [](Scalar x) { using alps::alea::sinh; return sinh(std::move(x)); })
        .def("cosh", [](Scalar x) { using alps::alea::cosh; return cosh(std::move(x)); })
        .def("tanh", [](Scalar x) { using alps::alea::tanh; return tanh(std::move(x)); })
        .def("set_bin_size",   &Scalar::set_bin_size)
        .def("set_bin_number", &Scalar::set_bin_number)
        .def("discard_bins",   &Scalar::discard_bins)
        .def("merge", static_cast<void (Scalar::*)(Scalar const &)>(&Scalar::merge))
        .def("save", static_cast<void (Scalar::*)(std::string const &, std::string const &) const>(&Scalar::save),
             nb::arg("filename"), nb::arg("observable_name"))
        .def("load", static_cast<void (Scalar::*)(std::string const &, std::string const &)>(&Scalar::load),
             nb::arg("filename"), nb::arg("observable_name"));
    nb::class_<Vector>(m, "MCVectorData",
        "Vector-valued Monte Carlo data.")
        .def(nb::init<>())
        .def(nb::init<std::vector<double>>(), nb::arg("mean"))
        .def(nb::init<std::vector<double>, std::vector<double>>(),
             nb::arg("mean"), nb::arg("error"))
        .def("__len__",
             [](Vector & v) {
                 return v.mean().size();
             })
        .def("__getitem__", &pymod::mcdata_vector_getitem<double>)
        .def("__contains__", &pymod::mcdata_vector_contains<double>)
        .def_prop_ro("mean",     [](Vector const & v) {
                                      return pymod::vec_to_numpy(v.mean());
                                  })
        .def_prop_ro("error",    [](Vector const & v) {
                                      return pymod::vec_to_numpy(v.error());
                                  })
        .def_prop_ro("tau",      [](Vector const & v) {
                                      return pymod::vec_to_numpy(v.tau());
                                  })
        .def_prop_ro("variance", [](Vector const & v) {
                                      return pymod::vec_to_numpy(v.variance());
                                  })
        .def_prop_ro("bins",     [](Vector const & v) {
                                      return pymod::matrix_to_numpy(v.bins());
                                  })
        .def_prop_ro("jackknife",[](Vector const & v) {
                                      return pymod::matrix_to_numpy(v.jackknife());
                                  })
        .def_prop_ro("count", &Vector::count)
        .def("__repr__",   &pymod::print_mcdata_vector<double>)
        .def("__format__", &pymod::format_mcdata_vector<double>,
             nb::arg("format_spec"))
        .def("__deepcopy__",
             [](Vector const & self, nb::handle /*memo*/) {
                 return Vector(self);
             })
        .def("__abs__", [](Vector x) { using alps::alea::abs; return abs(std::move(x)); })
        .def("__pow__", [](Vector x, double e) { using alps::alea::pow; return pow(std::move(x), e); })
        .def("__pos__", [](Vector self) { return +self; })
        .def("__neg__", [](Vector self) { return -self; })
        .def("__eq__",  [](Vector const & a, Vector const & b) { return a == b; }, nb::is_operator())
        // In-place operators.
        .def("__iadd__", [](Vector & s, Vector const & o)              -> Vector & { return s += o; }, nb::is_operator())
        .def("__iadd__", [](Vector & s, std::vector<double> const & o) -> Vector & { return s += o; }, nb::is_operator())
        .def("__isub__", [](Vector & s, Vector const & o)              -> Vector & { return s -= o; }, nb::is_operator())
        .def("__isub__", [](Vector & s, std::vector<double> const & o) -> Vector & { return s -= o; }, nb::is_operator())
        .def("__imul__", [](Vector & s, Vector const & o)              -> Vector & { return s *= o; }, nb::is_operator())
        .def("__imul__", [](Vector & s, std::vector<double> const & o) -> Vector & { return s *= o; }, nb::is_operator())
        .def("__itruediv__", [](Vector & s, Vector const & o)              -> Vector & { return s /= o; }, nb::is_operator())
        .def("__itruediv__", [](Vector & s, std::vector<double> const & o) -> Vector & { return s /= o; }, nb::is_operator())
        // Binary operators — Vector ↔ Vector / Scalar / vector<double> / double.
        .def("__add__",  [](Vector const & a, Vector const & b)              { return a + b; }, nb::is_operator())
        .def("__radd__", [](Vector const & a, Vector const & b)              { return b + a; }, nb::is_operator())
        .def("__add__",  [](Vector const & a, std::vector<double> const & b) { return a + b; }, nb::is_operator())
        .def("__radd__", [](Vector const & a, std::vector<double> const & b) { return b + a; }, nb::is_operator())
        .def("__sub__",  [](Vector const & a, Vector const & b)              { return a - b; }, nb::is_operator())
        .def("__rsub__", [](Vector const & a, Vector const & b)              { return b - a; }, nb::is_operator())
        .def("__sub__",  [](Vector const & a, std::vector<double> const & b) { return a - b; }, nb::is_operator())
        .def("__rsub__", [](Vector const & a, std::vector<double> const & b) { return b - a; }, nb::is_operator())
        .def("__mul__",  [](Vector const & a, Vector const & b)              { return a * b; }, nb::is_operator())
        .def("__mul__",  [](Vector const & a, Scalar const & b)              { return a * b; }, nb::is_operator())
        .def("__rmul__", [](Vector const & a, Scalar const & b)              { return b * a; }, nb::is_operator())
        .def("__rmul__", [](Vector const & a, Vector const & b)              { return b * a; }, nb::is_operator())
        .def("__mul__",  [](Vector const & a, std::vector<double> const & b) { return a * b; }, nb::is_operator())
        .def("__rmul__", [](Vector const & a, std::vector<double> const & b) { return b * a; }, nb::is_operator())
        .def("__truediv__",  [](Vector const & a, Vector const & b)              { return a / b; }, nb::is_operator())
        .def("__truediv__",  [](Vector const & a, Scalar const & b)              { return a / b; }, nb::is_operator())
        .def("__rtruediv__", [](Vector const & a, Vector const & b)              { return b / a; }, nb::is_operator())
        .def("__truediv__",  [](Vector const & a, std::vector<double> const & b) { return a / b; }, nb::is_operator())
        .def("__rtruediv__", [](Vector const & a, std::vector<double> const & b) { return b / a; }, nb::is_operator())
        // Vector ↔ double.
        .def("__add__",  [](Vector const & a, double b) { return a + b; }, nb::is_operator())
        .def("__radd__", [](Vector const & a, double b) { return b + a; }, nb::is_operator())
        .def("__sub__",  [](Vector const & a, double b) { return a - b; }, nb::is_operator())
        .def("__rsub__", [](Vector const & a, double b) { return b - a; }, nb::is_operator())
        .def("__mul__",  [](Vector const & a, double b) { return a * b; }, nb::is_operator())
        .def("__rmul__", [](Vector const & a, double b) { return b * a; }, nb::is_operator())
        .def("__truediv__",  [](Vector const & a, double b) { return a / b; }, nb::is_operator())
        .def("__rtruediv__", [](Vector const & a, double b) { return b / a; }, nb::is_operator())
        .def("sq",   [](Vector x) { using alps::alea::sq;   return sq(std::move(x));   })
        .def("cb",   [](Vector x) { using alps::alea::cb;   return cb(std::move(x));   })
        .def("sqrt", [](Vector x) { using alps::alea::sqrt; return sqrt(std::move(x)); })
        .def("cbrt", [](Vector x) { using alps::alea::cbrt; return cbrt(std::move(x)); })
        .def("exp",  [](Vector x) { using alps::alea::exp;  return exp(std::move(x));  })
        .def("log",  [](Vector x) { using alps::alea::log;  return log(std::move(x));  })
        .def("sin",  [](Vector x) { using alps::alea::sin;  return sin(std::move(x));  })
        .def("cos",  [](Vector x) { using alps::alea::cos;  return cos(std::move(x));  })
        .def("tan",  [](Vector x) { using alps::alea::tan;  return tan(std::move(x));  })
        .def("sinh", [](Vector x) { using alps::alea::sinh; return sinh(std::move(x)); })
        .def("cosh", [](Vector x) { using alps::alea::cosh; return cosh(std::move(x)); })
        .def("tanh", [](Vector x) { using alps::alea::tanh; return tanh(std::move(x)); })
        .def("set_bin_size",   &Vector::set_bin_size)
        .def("set_bin_number", &Vector::set_bin_number)
        .def("discard_bins",   &Vector::discard_bins)
        .def("merge", static_cast<void (Vector::*)(Vector const &)>(&Vector::merge))
        .def("save", static_cast<void (Vector::*)(std::string const &, std::string const &) const>(&Vector::save),
             nb::arg("filename"), nb::arg("observable_name"))
        .def("load", static_cast<void (Vector::*)(std::string const &, std::string const &)>(&Vector::load),
             nb::arg("filename"), nb::arg("observable_name"));
}
