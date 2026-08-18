// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Numpy interop without numpy headers. Construct numpy.ndarray
// instances from C++ buffers and consume incoming numpy arrays through
// nb::ndarray's DLPack/buffer view. The numpy package itself is loaded
// at runtime via nb::module_::import_("numpy"); pyalps already requires
// numpy as a runtime dependency.
#ifndef ALPS_PYTHON_NUMPY_COMPAT_HPP
#define ALPS_PYTHON_NUMPY_COMPAT_HPP
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <atomic>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>
namespace alps {
    namespace python {
        namespace nb_ = nanobind;
        // numpy dtype strings, indexed by the corresponding C++ type.
        // Used by make_numpy_array() / as_contiguous() to drive the
        // numpy.empty(dtype=…) / numpy.ascontiguousarray(dtype=…) calls.
        template <typename T> struct numpy_dtype;
        template <> struct numpy_dtype<bool>                 { static constexpr char const* name = "bool"; };
        template <> struct numpy_dtype<std::int8_t>          { static constexpr char const* name = "int8"; };
        template <> struct numpy_dtype<std::int16_t>         { static constexpr char const* name = "int16"; };
        template <> struct numpy_dtype<std::int32_t>         { static constexpr char const* name = "int32"; };
        template <> struct numpy_dtype<std::int64_t>         { static constexpr char const* name = "int64"; };
        template <> struct numpy_dtype<std::uint8_t>         { static constexpr char const* name = "uint8"; };
        template <> struct numpy_dtype<std::uint16_t>        { static constexpr char const* name = "uint16"; };
        template <> struct numpy_dtype<std::uint32_t>        { static constexpr char const* name = "uint32"; };
        template <> struct numpy_dtype<std::uint64_t>        { static constexpr char const* name = "uint64"; };
        template <> struct numpy_dtype<float>                { static constexpr char const* name = "float32"; };
        template <> struct numpy_dtype<double>               { static constexpr char const* name = "float64"; };
        template <> struct numpy_dtype<std::complex<float>>  { static constexpr char const* name = "complex64"; };
        template <> struct numpy_dtype<std::complex<double>> { static constexpr char const* name = "complex128"; };
        // Cached numpy module. Importing per call was a sys.modules
        // lookup + import-lock acquisition on every array conversion.
        // Not a function-local static: the winning thread's import can
        // release the GIL, so blocking a GIL-holding second thread on
        // the C++ static-init guard would deadlock. With the atomic
        // double-check, racing first callers both import (idempotent
        // under the import lock) and the loser drops its reference.
        // The winning reference is deliberately leaked so it stays
        // valid until interpreter shutdown regardless of static
        // destruction order.
        inline nb_::handle numpy_module() {
            static std::atomic<PyObject *> cached{nullptr};
            PyObject * mod = cached.load(std::memory_order_acquire);
            if (!mod) {
                mod = nb_::module_::import_("numpy").release().ptr();
                PyObject * expected = nullptr;
                if (!cached.compare_exchange_strong(expected, mod,
                                                    std::memory_order_acq_rel)) {
                    Py_DECREF(mod);
                    mod = expected;
                }
            }
            return mod;
        }
        // Allocates numpy.empty(shape, dtype=numpy_dtype<T>::name) and
        // memcpy's `data` (length = product(shape)) into it. Returns
        // a writable numpy.ndarray.
        template <typename T>
        inline nb_::object make_numpy_array(T const* data,
                                            std::vector<std::size_t> const& shape) {
            nb_::handle np = numpy_module();
            nb_::tuple shape_tuple = nb_::steal<nb_::tuple>(PyTuple_New(static_cast<Py_ssize_t>(shape.size())));
            if (!shape_tuple.is_valid())
                throw nb_::python_error();
            // PyTuple_SetItem (not the SET_ITEM macro): the macro pokes
            // tuple internals directly and is unavailable under the
            // limited API, which is otherwise within reach for these
            // bindings. SetItem steals the reference to dim.
            for (std::size_t i = 0; i < shape.size(); ++i) {
                PyObject * dim = PyLong_FromUnsignedLongLong(shape[i]);
                if (!dim)
                    throw nb_::python_error();
                PyTuple_SetItem(shape_tuple.ptr(), static_cast<Py_ssize_t>(i), dim);
            }
            nb_::object arr = np.attr("empty")(
                shape_tuple, nb_::arg("dtype") = numpy_dtype<T>::name);
            // Bridge the freshly-allocated numpy buffer through nb::ndarray
            // to get a writable raw pointer.
            auto nd = nb_::cast<nb_::ndarray<T, nb_::c_contig>>(arr);
            std::size_t total = 1;
            for (auto s : shape) total *= s;
            if (total > 0)
                std::memcpy(nd.data(), data, total * sizeof(T));
            return arr;
        }
        template <typename T>
        inline nb_::object make_numpy_array(std::vector<T> const& v) {
            return make_numpy_array<T>(v.data(), {v.size()});
        }
        // Strong-ref'd C-contiguous view onto a numpy array of dtype T.
        // The owner handle keeps the array alive for the lifetime of
        // the view; data() / shape() / ndim() forward to the ndarray.
        template <typename T>
        struct contiguous_view {
            nb_::object owner;
            nb_::ndarray<T, nb_::c_contig> nd;
            T const* data() const { return nd.data(); }
            std::size_t ndim() const { return nd.ndim(); }
            std::size_t shape(int i) const { return nd.shape(i); }
        };
        // Coerces `obj` to a C-contiguous numpy.ndarray of dtype T via
        // numpy.ascontiguousarray. Always produces a contiguous +
        // correctly-typed buffer (numpy copies if the input doesn't
        // already match). Equivalent in spirit to nanobind's
        // py::array_t<T, py::array::c_style | py::array::forcecast>
        // parameter form, just routed through numpy at runtime instead
        // of through the numpy C headers at compile time.
        template <typename T>
        inline contiguous_view<T> as_contiguous(nb_::handle obj) {
            nb_::handle np = numpy_module();
            nb_::object arr = np.attr("ascontiguousarray")(
                obj, nb_::arg("dtype") = numpy_dtype<T>::name);
            auto nd = nb_::cast<nb_::ndarray<T, nb_::c_contig>>(arr);
            return contiguous_view<T>{std::move(arr), std::move(nd)};
        }
    } // namespace python
} // namespace alps
#endif // ALPS_PYTHON_NUMPY_COMPAT_HPP
