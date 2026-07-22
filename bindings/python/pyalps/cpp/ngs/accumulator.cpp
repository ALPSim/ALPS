// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <alps/ngs/accumulator/feature.hpp>
#include <alps/ngs/accumulator/feature/count.hpp>
#include <alps/ngs/accumulator/feature/mean.hpp>
#include <alps/ngs/accumulator/feature/error.hpp>
#include <alps/ngs/accumulator/feature/binning_analysis.hpp>
#include <alps/ngs/accumulator/feature/max_num_binning.hpp>
#include <alps/hdf5/vector.hpp>
#include <nanobind/nanobind.h>
#include <sstream>
#include <string>
namespace nb = nanobind;
namespace {
template <typename T>
std::string print_value(T const & value) {
    std::stringstream stream;
    value.print(stream);
    return stream.str();
}
template <typename Accumulator>
typename Accumulator::result_type make_result(Accumulator const & value) {
    return typename Accumulator::result_type(value);
}
template <typename Result>
void bind_result_operators(nb::class_<Result> & cls) {
    cls
        .def("__neg__", [](Result value) { value.negate(); return value; })
        .def("__iadd__", [](Result & self, Result const & other) -> Result & { self += other; return self; }, nb::rv_policy::reference_internal)
        .def("__iadd__", [](Result & self, double value) -> Result & { self += value; return self; }, nb::rv_policy::reference_internal)
        .def("__isub__", [](Result & self, Result const & other) -> Result & { self -= other; return self; }, nb::rv_policy::reference_internal)
        .def("__isub__", [](Result & self, double value) -> Result & { self -= value; return self; }, nb::rv_policy::reference_internal)
        .def("__imul__", [](Result & self, Result const & other) -> Result & { self *= other; return self; }, nb::rv_policy::reference_internal)
        .def("__imul__", [](Result & self, double value) -> Result & { self *= value; return self; }, nb::rv_policy::reference_internal)
        .def("__itruediv__", [](Result & self, Result const & other) -> Result & { self /= other; return self; }, nb::rv_policy::reference_internal)
        .def("__itruediv__", [](Result & self, double value) -> Result & { self /= value; return self; }, nb::rv_policy::reference_internal)
        .def("__add__", [](Result value, Result const & other) { value += other; return value; }, nb::is_operator())
        .def("__add__", [](Result value, double other) { value += other; return value; }, nb::is_operator())
        .def("__radd__", [](Result value, double other) { value += other; return value; }, nb::is_operator())
        .def("__sub__", [](Result value, Result const & other) { value -= other; return value; }, nb::is_operator())
        .def("__sub__", [](Result value, double other) { value -= other; return value; }, nb::is_operator())
        .def("__rsub__", [](Result value, double other) { value.negate(); value += other; return value; }, nb::is_operator())
        .def("__mul__", [](Result value, Result const & other) { value *= other; return value; }, nb::is_operator())
        .def("__mul__", [](Result value, double other) { value *= other; return value; }, nb::is_operator())
        .def("__rmul__", [](Result value, double other) { value *= other; return value; }, nb::is_operator())
        .def("__truediv__", [](Result value, Result const & other) { value /= other; return value; }, nb::is_operator())
        .def("__truediv__", [](Result value, double other) { value /= other; return value; }, nb::is_operator())
        .def("__rtruediv__", [](Result value, double other) { value.inverse(); value *= other; return value; }, nb::is_operator())
        .def("sin",  [](Result value) { value.sin();  return value; })
        .def("cos",  [](Result value) { value.cos();  return value; })
        .def("tan",  [](Result value) { value.tan();  return value; })
        .def("sinh", [](Result value) { value.sinh(); return value; })
        .def("cosh", [](Result value) { value.cosh(); return value; })
        .def("tanh", [](Result value) { value.tanh(); return value; })
        .def("abs",  [](Result value) { value.abs();  return value; })
        .def("sqrt", [](Result value) { value.sqrt(); return value; })
        .def("log",  [](Result value) { value.log();  return value; });
}
template <typename T>
void bind_serializable(nb::class_<T> & cls) {
    cls.def("__str__", &print_value<T>)
       .def("save", &T::save)
       .def("load", &T::load)
       .def("reset", &T::reset);
}
} // namespace
NB_MODULE(pyngsaccumulator_c, m) {
    using namespace alps::accumulator::impl;
    using count_accumulator = Accumulator<double, alps::accumulator::count_tag, AccumulatorBase<double>>;
    using count_result = count_accumulator::result_type;
    nb::class_<count_accumulator> count_acc(m, "count_accumulator");
    count_acc.def(nb::init<>()).def("__call__", [](count_accumulator & self, double value) { self(value); })
             .def("result", &make_result<count_accumulator>).def("count", &count_accumulator::count);
    bind_serializable(count_acc);
    nb::class_<count_result> count_res(m, "count_result");
    count_res.def(nb::init<>()).def("count", &count_result::count);
    bind_serializable(count_res); bind_result_operators(count_res);
    using mean_accumulator = Accumulator<double, alps::accumulator::mean_tag, count_accumulator>;
    using mean_result = mean_accumulator::result_type;
    nb::class_<mean_accumulator> mean_acc(m, "mean_accumulator");
    mean_acc.def(nb::init<>()).def("__call__", [](mean_accumulator & self, double value) { self(value); })
            .def("result", &make_result<mean_accumulator>).def("count", &mean_accumulator::count)
            .def("mean", &mean_accumulator::mean);
    bind_serializable(mean_acc);
    nb::class_<mean_result> mean_res(m, "mean_result");
    mean_res.def(nb::init<>()).def("count", &mean_result::count).def("mean", &mean_result::mean);
    bind_serializable(mean_res); bind_result_operators(mean_res);
    using error_accumulator = Accumulator<double, alps::accumulator::error_tag, mean_accumulator>;
    using error_result = error_accumulator::result_type;
    nb::class_<error_accumulator> error_acc(m, "error_accumulator");
    error_acc.def(nb::init<>()).def("__call__", [](error_accumulator & self, double value) { self(value); })
             .def("result", &make_result<error_accumulator>).def("count", &error_accumulator::count)
             .def("mean", &error_accumulator::mean).def("error", &error_accumulator::error);
    bind_serializable(error_acc);
    nb::class_<error_result> error_res(m, "error_result");
    error_res.def(nb::init<>()).def("count", &error_result::count).def("mean", &error_result::mean)
             .def("error", &error_result::error);
    bind_serializable(error_res); bind_result_operators(error_res);
    using binning_accumulator = Accumulator<double, alps::accumulator::binning_analysis_tag, error_accumulator>;
    using binning_result = binning_accumulator::result_type;
    nb::class_<binning_accumulator> binning_acc(m, "binning_analysis_accumulator");
    binning_acc.def(nb::init<>()).def("__call__", [](binning_accumulator & self, double value) { self(value); })
               .def("result", &make_result<binning_accumulator>).def("count", &binning_accumulator::count)
               .def("mean", &binning_accumulator::mean).def("error", &binning_accumulator::error);
    bind_serializable(binning_acc);
    nb::class_<binning_result> binning_res(m, "binning_analysis_result");
    binning_res.def(nb::init<>()).def("count", &binning_result::count).def("mean", &binning_result::mean)
               .def("error", &binning_result::error);
    bind_serializable(binning_res); bind_result_operators(binning_res);
    using maxbin_accumulator = Accumulator<double, alps::accumulator::max_num_binning_tag, error_accumulator>;
    using maxbin_result = maxbin_accumulator::result_type;
    nb::class_<maxbin_accumulator> maxbin_acc(m, "max_num_binning_accumulator");
    maxbin_acc.def(nb::init<>()).def("__call__", [](maxbin_accumulator & self, double value) { self(value); })
              .def("result", &make_result<maxbin_accumulator>).def("count", &maxbin_accumulator::count)
              .def("mean", &maxbin_accumulator::mean).def("error", &maxbin_accumulator::error);
    bind_serializable(maxbin_acc);
    nb::class_<maxbin_result> maxbin_res(m, "max_num_binning_result");
    maxbin_res.def(nb::init<>()).def("count", &maxbin_result::count).def("mean", &maxbin_result::mean)
              .def("error", &maxbin_result::error);
    bind_serializable(maxbin_res); bind_result_operators(maxbin_res);
}
