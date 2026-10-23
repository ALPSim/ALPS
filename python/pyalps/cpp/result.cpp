// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// pyalps.ngs.result: a result of ALPSCore's accumulators, held through
// pyalps::result (see accumulator_handles.hpp).
#include <nanobind/nanobind.h>
#include "accumulator_handles.hpp"
#include "archive_savable.hpp"
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <alps/hdf5/archive.hpp>
#include <boost/function.hpp>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace nb = nanobind;
namespace acc = alps::accumulators;
using pyalps::result;
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
// Scalar results read back as floats, vector results as numpy arrays.
template <typename Scalar, typename Vector>
nb::object by_type(result const & self, Scalar scalar, Vector vector) {
    acc::result_wrapper & r = self.get();
    if (pyalps::holds<double>(r))
        return nb::float_(scalar(r));
    if (pyalps::holds<std::vector<double>>(r))
        return vec_to_numpy(vector(r));
    throw std::runtime_error("pyalps handles results of double or vector of double values");
}
result wrap(acc::result_wrapper && value) {
    return result{std::make_shared<acc::result_wrapper>(std::move(value))};
}
std::string print(result const & self) {
    if (self.empty())
        return "No Measurements";
    std::ostringstream out;
    self.value->print(out, true);
    return out.str();
}
// ALPSCore's results have no exp() or pow(); apply them as transformations,
// which propagate errors through the jackknife where the result has bins.
result apply(result const & self, double (*op)(double)) {
    acc::result_wrapper & r = self.get();
    if (pyalps::holds<double>(r))
        return wrap(r.transform<double>(boost::function<double(double)>(op)));
    return wrap(r.transform<std::vector<double>>(boost::function<std::vector<double>(std::vector<double>)>(
        [op](std::vector<double> values) {
            for (double & value : values) value = op(value);
            return values;
        })));
}
result power(result const & self, double exponent) {
    acc::result_wrapper & r = self.get();
    if (pyalps::holds<double>(r))
        return wrap(r.transform<double>(boost::function<double(double)>(
            [exponent](double value) { return std::pow(value, exponent); })));
    return wrap(r.transform<std::vector<double>>(boost::function<std::vector<double>(std::vector<double>)>(
        [exponent](std::vector<double> values) {
            for (double & value : values) value = std::pow(value, exponent);
            return values;
        })));
}
double exponential(double value) { return std::exp(value); }
} // namespace
NB_MODULE(pyngsresult_c, m) {
    nb::module_::import_("pyalps.cxx.pyngsobservable_c");
    m.def("observable2result", [](acc::accumulator_wrapper const & observable) {
        return result{observable.result()};
    });
    using R = result;
#define PYALPS_RESULT_FUNCTION(NAME) \
    .def(#NAME, [](R const & self) { return wrap(acc::NAME(self.get())); })
    nb::class_<R>(m, "result")
        .def(nb::init<>())
        // An independent copy, as copy.deepcopy() gives.
        .def("__init__", [](R * self, R const & other) {
            new (self) R{other.value ? std::make_shared<acc::result_wrapper>(*other.value) : nullptr};
        })
        .def("__repr__", &print)
        .def("__deepcopy__", [](R const & self, nb::handle /*memo*/) {
            return R{self.value ? std::make_shared<acc::result_wrapper>(*self.value) : nullptr};
        })
        .def_prop_ro("mean", [](R const & self) {
            return by_type(self, [](acc::result_wrapper & r) { return r.mean<double>(); },
                           [](acc::result_wrapper & r) { return r.mean<std::vector<double>>(); });
        })
        .def_prop_ro("error", [](R const & self) {
            return by_type(self, [](acc::result_wrapper & r) { return r.error<double>(); },
                           [](acc::result_wrapper & r) { return r.error<std::vector<double>>(); });
        })
        // The integrated autocorrelation time, for results with binning.
        .def_prop_ro("tau", [](R const & self) {
            return by_type(self, [](acc::result_wrapper & r) { return r.autocorrelation<double>(); },
                           [](acc::result_wrapper & r) { return r.autocorrelation<std::vector<double>>(); });
        })
        // The bin means of a scalar result with full binning.
        .def_prop_ro("bins", [](R const & self) {
            acc::result_wrapper & r = self.get();
            if (!pyalps::holds<double>(r))
                throw std::runtime_error("bins are available for scalar results");
            return vec_to_numpy(
                r.extract<acc::FullBinningAccumulator<double>::result_type>().max_num_binning().bins());
        })
        .def_prop_ro("count", [](R const & self) -> std::size_t {
            return self.value ? self.value->count() : 0;
        })
        .def("__pos__", [](R const & self) { return wrap(+self.get()); })
        .def("__neg__", [](R const & self) { return wrap(-self.get()); })
        .def("__abs__", [](R const & self) { return wrap(acc::abs(self.get())); })
        .def("__pow__", &power)
        // In-place operators modify the result, which every reference to it
        // sees, and return the same Python object (Python's __i*__
        // semantics). rv_policy::none finds the existing wrapper of self.
        .def("__iadd__", [](R & self, R const & o) -> R & { self.get() += o.get(); return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__iadd__", [](R & self, double o)    -> R & { self.get() += o; return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__isub__", [](R & self, R const & o) -> R & { self.get() -= o.get(); return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__isub__", [](R & self, double o)    -> R & { self.get() -= o; return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__imul__", [](R & self, R const & o) -> R & { self.get() *= o.get(); return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__imul__", [](R & self, double o)    -> R & { self.get() *= o; return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__itruediv__", [](R & self, R const & o) -> R & { self.get() /= o.get(); return self; }, nb::rv_policy::none, nb::is_operator())
        .def("__itruediv__", [](R & self, double o)    -> R & { self.get() /= o; return self; }, nb::rv_policy::none, nb::is_operator())
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
        PYALPS_RESULT_FUNCTION(sq)
        PYALPS_RESULT_FUNCTION(cb)
        PYALPS_RESULT_FUNCTION(sqrt)
        PYALPS_RESULT_FUNCTION(cbrt)
        PYALPS_RESULT_FUNCTION(log)
        PYALPS_RESULT_FUNCTION(sin)
        PYALPS_RESULT_FUNCTION(cos)
        PYALPS_RESULT_FUNCTION(tan)
        PYALPS_RESULT_FUNCTION(sinh)
        PYALPS_RESULT_FUNCTION(cosh)
        PYALPS_RESULT_FUNCTION(tanh)
        .def("exp", [](R const & self) { return apply(self, &exponential); })
        .def("save", [](R const & self, alps::hdf5::archive & ar) { self.get().save(ar); })
        // Loads the result at the archive's context. A failed load leaves
        // the result as it was; other references keep the old value.
        .def("load", [](R & self, alps::hdf5::archive & ar) { self.value = pyalps::load_result(ar); });
#undef PYALPS_RESULT_FUNCTION
    pyalps::mark_archive_savable(m.attr("result"));
}
