// Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
//               Lukas Gamper <gamperl@gmail.com>,
//               Matthias Troyer <troyer@itp.phys.ethz.ch>,
//               Maximilian Poprawe <poprawem@ethz.ch>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <alps/alea/detailedbinning.h>
#include <alps/alea/mcanalyze.hpp>
#include <alps/alea/mcdata.hpp>
#include <alps/alea/value_with_error.hpp>
#include <alps/hdf5.hpp>
#include <alps/numeric/vector_functions.hpp>
#include "numpy_compat.hpp"
#include "save_observable_to_hdf5.hpp"
#include <alps/random.h>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <valarray>
#include <vector>
namespace nb = nanobind;
namespace alps {
    namespace alea {
        // Wraps a scalar-valarray observable (RealVectorObservable /
        // RealVectorTimeSeriesObservable) so Python sees numpy arrays
        // in / out instead of std::valarray<double>.
        template<typename T>
        class WrappedValarrayObservable {
            using element_type = typename T::value_type::value_type;
          public:
            WrappedValarrayObservable(std::string const & name, int s = 0)
                : obs(name, s)
            {}
            // Copy the ndarray into a valarray and feed it to the observable.
            void push(nb::handle arr) {
                auto view = alps::python::as_contiguous<double>(arr);
                if (view.ndim() != 1)
                    throw std::invalid_argument(
                        "RealVectorObservable.push: expected 1-D array");
                std::size_t n = static_cast<std::size_t>(view.shape(0));
                std::valarray<element_type> v(n);
                double const * data = view.data();
                for (std::size_t i = 0; i < n; ++i)
                    v[i] = static_cast<element_type>(data[i]);
                obs << v;
            }
            std::string representation() const { return obs.representation(); }
            // Turn an alps::numeric std::valarray-like view into a
            // 1-D numpy.ndarray (dtype=float64) by copying through
            // numpy.empty + buffer protocol.
            template <typename U>
            static nb::object _to_numpy(U const & v) {
                std::size_t n = static_cast<std::size_t>(v.size());
                std::vector<double> tmp(n);
                for (std::size_t i = 0; i < n; ++i)
                    tmp[i] = static_cast<double>(v[i]);
                return alps::python::make_numpy_array<double>(tmp.data(), {n});
            }
            nb::object mean()     const { return _to_numpy(obs.mean());     }
            nb::object error()    const { return _to_numpy(obs.error());    }
            nb::object tau()      const { return _to_numpy(obs.tau());      }
            nb::object variance() const { return _to_numpy(obs.variance()); }
            void save(std::string const & filename) const {
                alps::hdf5::archive ar(filename, "a");
                ar["/simulation/results/" + obs.representation()] << obs;
            }
            typename T::count_type       count()            const { return obs.count(); }
            typename T::convergence_type converged_errors() const { return obs.converged_errors(); }
          private:
            T obs;
        };
    } // namespace alea
} // namespace alps
namespace {
// Build a 1-D numpy.ndarray (dtype=float64) from any sequence-like
// alps container (std::vector<double>, std::valarray<double>).
template <typename Container>
nb::object seq_to_numpy(Container const & v) {
    std::size_t n = static_cast<std::size_t>(v.size());
    std::vector<double> tmp(n);
    for (std::size_t i = 0; i < n; ++i)
        tmp[i] = static_cast<double>(v[i]);
    return alps::python::make_numpy_array<double>(tmp.data(), {n});
}
// Copy a numpy array into a std::vector<ValueType>. Used when the
// caller still instantiates mctimeseries<ValueType> from a
// Python array.
template <typename ValueType>
std::vector<ValueType>
numpy_to_vector(nb::handle arr) {
    auto view = alps::python::as_contiguous<double>(arr);
    if (view.ndim() != 1)
        throw std::invalid_argument(
            "mctimeseries ctor: expected 1-D array");
    std::size_t n = static_cast<std::size_t>(view.shape(0));
    std::vector<ValueType> out(n);
    double const * data = view.data();
    for (std::size_t i = 0; i < n; ++i)
        out[i] = static_cast<ValueType>(data[i]);
    return out;
}
// __repr__ helper for any ALPS type with an ostream operator.
template <typename T>
std::string stream_repr(T const & x) {
    std::ostringstream ss;
    ss << x;
    return ss.str();
}
template <typename T>
std::string value_with_error_repr(alps::alea::value_with_error<T> const & v) {
    std::ostringstream ss;
    ss << v.mean() << " +/- " << v.error();
    return ss.str();
}
// Numpy-returning wrappers for vector-valued alps::alea free
// functions (mean, variance, uncorrelated_error, binning_error) —
// the scalar overloads are bound directly and nanobind casts
// their `double` return to a Python float automatically.
template <typename T>
nb::object mean_vector(T const & x) {
    return seq_to_numpy(alps::alea::mean(x));
}
template <typename T>
nb::object variance_vector(T const & x) {
    return seq_to_numpy(alps::alea::variance(x));
}
// mctimeseries.timeseries() returns std::vector<ValueType>; hand
// back to Python as numpy. For scalar ValueType we pack 1-D; for
// vector<double> ValueType we pack 2-D. mctimeseries_view has the
// same surface.
template <typename ValueType, typename TS>
nb::object ts_to_numpy_scalar(TS const & ts) {
    auto const & v = ts.timeseries();
    return seq_to_numpy(v);
}
template <typename TS>
nb::object ts_to_numpy_vector_rows(TS const & ts) {
    auto const & rows = ts.timeseries();
    if (rows.empty())
        return alps::python::make_numpy_array<double>(
            static_cast<double const*>(nullptr), {std::size_t{0}, std::size_t{0}});
    std::size_t nrows = rows.size();
    std::size_t ncols = rows.front().size();
    for (auto const & row : rows)
        if (row.size() != ncols)
            throw std::runtime_error("mctimeseries has ragged rows; cannot shape as numpy 2-D");
    std::vector<double> flat(nrows * ncols);
    double * dst = flat.data();
    for (auto const & row : rows) {
        for (std::size_t j = 0; j < ncols; ++j)
            *dst++ = static_cast<double>(row[j]);
    }
    return alps::python::make_numpy_array<double>(flat.data(), {nrows, ncols});
}
} // namespace
NB_MODULE(pyalea_c, m) {
    m.doc() = "ALPS alea bindings (nanobind)";
    // ─── scalar-valarray observables ─────────────────────────────────
    using RealVecObs = alps::alea::WrappedValarrayObservable<alps::RealVectorObservable>;
    using RealVecTsObs = alps::alea::WrappedValarrayObservable<alps::RealVectorTimeSeriesObservable>;
    #define ALPS_PY_EXPORT_VECTOROBSERVABLE(Wrapper, PyName)                              \
        nb::class_<Wrapper>(m, PyName)                                                    \
            .def("__init__",                                                              \
                 [](Wrapper * self, std::string name, int bins) {                         \
                     new (self) Wrapper(name, bins);                                      \
                 },                                                                       \
                 nb::arg("name"), nb::arg("bins") = 0)                                    \
            .def("__repr__",     &Wrapper::representation)                                \
            .def("__deepcopy__",                                                          \
                 [](Wrapper const & self, nb::handle /*memo*/) {                          \
                     return Wrapper(self);                                                \
                 })                                                                       \
            .def("__lshift__",   &Wrapper::push, nb::arg("array"))                        \
            .def("save",         &Wrapper::save, nb::arg("filename"))                     \
            .def_prop_ro("mean",     &Wrapper::mean)                                      \
            .def_prop_ro("error",    &Wrapper::error)                                     \
            .def_prop_ro("tau",      &Wrapper::tau)                                       \
            .def_prop_ro("variance", &Wrapper::variance)                                  \
            .def_prop_ro("count",    &Wrapper::count)                                     \
            .def_prop_ro("converged_errors", &Wrapper::converged_errors)
    ALPS_PY_EXPORT_VECTOROBSERVABLE(RealVecObs,   "RealVectorObservable");
    ALPS_PY_EXPORT_VECTOROBSERVABLE(RealVecTsObs, "RealVectorTimeSeriesObservable");
    #undef ALPS_PY_EXPORT_VECTOROBSERVABLE
    // ─── scalar simple observables ───────────────────────────────────
    #define ALPS_PY_EXPORT_SIMPLEOBSERVABLE(AlpsClass, PyName)                            \
        nb::class_<alps::AlpsClass>(m, PyName)                                            \
            .def("__init__",                                                              \
                 [](alps::AlpsClass * self, std::string name, int bins) {                 \
                     new (self) alps::AlpsClass(name, bins);                              \
                 },                                                                       \
                 nb::arg("name"), nb::arg("bins") = 0)                                    \
            .def("__deepcopy__",                                                          \
                 [](alps::AlpsClass const & self, nb::handle /*memo*/) {                  \
                     return alps::AlpsClass(self);                                        \
                 })                                                                       \
            .def("__repr__",   &alps::AlpsClass::representation)                          \
            .def("__lshift__", &alps::AlpsClass::operator<<)                              \
            .def("save",       &alps::python::save_observable_to_hdf5<alps::AlpsClass>,   \
                 nb::arg("filename"))                                                     \
            .def_prop_ro("mean", &alps::AlpsClass::mean)                                  \
            .def_prop_ro("error",                                                         \
                static_cast<alps::AlpsClass::result_type (alps::AlpsClass::*)() const>(   \
                    &alps::AlpsClass::error))                                             \
            .def_prop_ro("tau",      &alps::AlpsClass::tau)                               \
            .def_prop_ro("variance", &alps::AlpsClass::variance)                          \
            .def_prop_ro("count",    &alps::AlpsClass::count)                             \
            .def_prop_ro("converged_errors", &alps::AlpsClass::converged_errors)
    ALPS_PY_EXPORT_SIMPLEOBSERVABLE(RealObservable,           "RealObservable");
    ALPS_PY_EXPORT_SIMPLEOBSERVABLE(RealTimeSeriesObservable, "RealTimeSeriesObservable");
    #undef ALPS_PY_EXPORT_SIMPLEOBSERVABLE
    // ─── value_with_error ────────────────────────────────────────────
    nb::class_<alps::alea::value_with_error<double>>(m, "ValueWithError")
        .def(nb::init<double, double>(),
             nb::arg("mean") = 0.0, nb::arg("error") = 0.0)
        .def_prop_ro("mean",  &alps::alea::value_with_error<double>::mean)
        .def_prop_ro("error", &alps::alea::value_with_error<double>::error)
        .def("__repr__", &value_with_error_repr<double>);
    // ─── StdPairDouble ─────────────────────────────────────────────
    // nanobind's STL caster already registered std::pair<double, double> as a
    // Python tuple converter via the stl.h header; nanobind takes the
    // same path via <nanobind/stl/pair.h>. Binding it again as a
    // class_ would fight that, so use a thin attribute-access
    // wrapper instead, and give integrated_autocorrelation_time a
    // Python-side signature that accepts either StdPairDouble or a
    // plain (float, float) tuple.
    struct StdPairDouble {
        double first{0.0};
        double second{0.0};
        StdPairDouble() = default;
        StdPairDouble(double f, double s) : first(f), second(s) {}
        operator std::pair<double, double>() const { return {first, second}; }
    };
    nb::class_<StdPairDouble>(m, "StdPairDouble",
        "Pair of (fit slope, fit intercept) returned by the autocorrelation fit helpers.")
        .def(nb::init<>())
        .def(nb::init<double, double>(), nb::arg("first"), nb::arg("second"))
        .def_rw("first",  &StdPairDouble::first)
        .def_rw("second", &StdPairDouble::second)
        .def("__repr__", [](StdPairDouble const & p) {
            std::ostringstream ss;
            ss << "StdPairDouble(" << p.first << ", " << p.second << ")";
            return ss.str();
        });
    // ─── mctimeseries<T> / mctimeseries_view<T> bindings ────────────
    //
    // Numpy-ctor: take any handle and copy through numpy_to_vector.
    // The timeseries() method returns numpy arrays directly; the
    // 2-D overload for vector<double> goes through the
    // ts_to_numpy_vector_rows helper.
    #define ALPS_PY_EXPORT_MCTIMESERIES_SCALAR(Value, PyName)                             \
        nb::class_<alps::alea::mctimeseries<Value>>(m, PyName)                            \
            .def(nb::init<>())                                                            \
            .def("__init__",                                                              \
                 [](alps::alea::mctimeseries<Value> * self, nb::handle a) {               \
                     new (self) alps::alea::mctimeseries<Value>(                          \
                         numpy_to_vector<Value>(a));                                      \
                 })                                                                       \
            .def(nb::init<alps::alea::mcdata<Value>>())                                   \
            .def("timeseries", [](alps::alea::mctimeseries<Value> const & self) {         \
                    return ts_to_numpy_scalar<Value>(self);                               \
                })                                                                        \
            .def_prop_ro("size", &alps::alea::mctimeseries<Value>::size)                  \
            .def("__repr__", &stream_repr<alps::alea::mctimeseries<Value>>);              \
        nb::class_<alps::alea::mctimeseries_view<Value>>(m, PyName "View")                \
            .def(nb::init<alps::alea::mctimeseries<Value>>())                             \
            .def(nb::init<alps::alea::mctimeseries_view<Value>>())                        \
            .def("timeseries", [](alps::alea::mctimeseries_view<Value> const & self) {    \
                    return ts_to_numpy_scalar<Value>(self);                               \
                })                                                                        \
            .def_prop_ro("size", &alps::alea::mctimeseries_view<Value>::size)             \
            .def("__repr__", &stream_repr<alps::alea::mctimeseries_view<Value>>)
    ALPS_PY_EXPORT_MCTIMESERIES_SCALAR(double, "MCScalarTimeseries");
    #undef ALPS_PY_EXPORT_MCTIMESERIES_SCALAR
    // Vector-valued mctimeseries: ctor from a 2-D numpy array, rows
    // are time samples. timeseries() returns 2-D.
    using VecTs    = alps::alea::mctimeseries<std::vector<double>>;
    using VecTsV   = alps::alea::mctimeseries_view<std::vector<double>>;
    using VecMcD   = alps::alea::mcdata<std::vector<double>>;
    nb::class_<VecTs>(m, "MCVectorTimeseries")
        .def(nb::init<>())
        .def("__init__", [](VecTs * self, nb::handle a) {
                auto view = alps::python::as_contiguous<double>(a);
                if (view.ndim() != 2)
                    throw std::invalid_argument(
                        "MCVectorTimeseries ctor: expected 2-D array");
                std::size_t nrows = static_cast<std::size_t>(view.shape(0));
                std::size_t ncols = static_cast<std::size_t>(view.shape(1));
                std::vector<std::vector<double>> rows(nrows);
                double const * data = view.data();
                for (std::size_t i = 0; i < nrows; ++i) {
                    rows[i].resize(ncols);
                    for (std::size_t j = 0; j < ncols; ++j)
                        rows[i][j] = data[i * ncols + j];
                }
                new (self) VecTs(rows);
            })
        .def(nb::init<VecMcD>())
        .def("timeseries", [](VecTs const & self) { return ts_to_numpy_vector_rows(self); })
        .def_prop_ro("size", &VecTs::size)
        .def("__repr__", &stream_repr<VecTs>);
    nb::class_<VecTsV>(m, "MCVectorTimeseriesView")
        .def(nb::init<VecTs>())
        .def(nb::init<VecTsV>())
        .def("timeseries", [](VecTsV const & self) { return ts_to_numpy_vector_rows(self); })
        .def_prop_ro("size", &VecTsV::size)
        .def("__repr__", &stream_repr<VecTsV>);
    // ─── alps::alea free functions over mcdata / mctimeseries ────────
    #define DEF_ALL(name, fn)                                                             \
        /* scalar-valued */                                                               \
        m.def(name, static_cast<double (*)(alps::alea::mcdata<double> const &)>(&fn));    \
        m.def(name, static_cast<double (*)(alps::alea::mctimeseries<double> const &)>(&fn)); \
        m.def(name, static_cast<double (*)(alps::alea::mctimeseries_view<double> const &)>(&fn));
    // size — works for both scalar and vector value types.
    m.def("size", static_cast<std::size_t (*)(alps::alea::mcdata<double> const &)>(&alps::size));
    m.def("size", static_cast<std::size_t (*)(alps::alea::mctimeseries<double> const &)>(&alps::size));
    m.def("size", static_cast<std::size_t (*)(alps::alea::mctimeseries_view<double> const &)>(&alps::size));
    m.def("size", static_cast<std::size_t (*)(alps::alea::mcdata<std::vector<double>> const &)>(&alps::size));
    m.def("size", static_cast<std::size_t (*)(alps::alea::mctimeseries<std::vector<double>> const &)>(&alps::size));
    m.def("size", static_cast<std::size_t (*)(alps::alea::mctimeseries_view<std::vector<double>> const &)>(&alps::size));
    // mean — scalar overloads return double, vector overloads return numpy.
    DEF_ALL("mean", alps::alea::mean)
    m.def("mean", &mean_vector<alps::alea::mcdata<std::vector<double>>>);
    m.def("mean", &mean_vector<alps::alea::mctimeseries<std::vector<double>>>);
    m.def("mean", &mean_vector<alps::alea::mctimeseries_view<std::vector<double>>>);
    // variance — same pattern.
    DEF_ALL("variance", alps::alea::variance)
    m.def("variance", &variance_vector<alps::alea::mcdata<std::vector<double>>>);
    m.def("variance", &variance_vector<alps::alea::mctimeseries<std::vector<double>>>);
    m.def("variance", &variance_vector<alps::alea::mctimeseries_view<std::vector<double>>>);
    // integrated_autocorrelation_time — scalar only. The C++ signature
    // takes the (slope, intercept) pair by const-ref.
    m.def("integrated_autocorrelation_time",
          static_cast<double (*)(alps::alea::mctimeseries<double> const &,
                                 std::pair<double, double> const &)>(
              &alps::alea::integrated_autocorrelation_time));
    m.def("integrated_autocorrelation_time",
          static_cast<double (*)(alps::alea::mctimeseries_view<double> const &,
                                 std::pair<double, double> const &)>(
              &alps::alea::integrated_autocorrelation_time));
    // running_mean / reverse_running_mean — scalar only, mctimeseries-valued.
    m.def("running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mcdata<double> const &)>(
              &alps::alea::running_mean));
    m.def("running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mctimeseries<double> const &)>(
              &alps::alea::running_mean));
    m.def("running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mctimeseries_view<double> const &)>(
              &alps::alea::running_mean));
    m.def("reverse_running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mcdata<double> const &)>(
              &alps::alea::reverse_running_mean));
    m.def("reverse_running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mctimeseries<double> const &)>(
              &alps::alea::reverse_running_mean));
    m.def("reverse_running_mean",
          static_cast<alps::alea::mctimeseries<double> (*)(alps::alea::mctimeseries_view<double> const &)>(
              &alps::alea::reverse_running_mean));
    #undef DEF_ALL
}
