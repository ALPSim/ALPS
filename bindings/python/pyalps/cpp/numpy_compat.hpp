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
        // Allocates numpy.empty(shape, dtype=numpy_dtype<T>::name) and
        // memcpy's `data` (length = product(shape)) into it. Returns
        // a writable numpy.ndarray.
        template <typename T>
        inline nb_::object make_numpy_array(T const* data,
                                            std::vector<std::size_t> const& shape) {
            nb_::object np = nb_::module_::import_("numpy");
            nb_::tuple shape_tuple = nb_::steal<nb_::tuple>(PyTuple_New(static_cast<Py_ssize_t>(shape.size())));
            for (std::size_t i = 0; i < shape.size(); ++i)
                PyTuple_SET_ITEM(shape_tuple.ptr(), static_cast<Py_ssize_t>(i),
                                 PyLong_FromUnsignedLongLong(shape[i]));
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
            nb_::object np = nb_::module_::import_("numpy");
            nb_::object arr = np.attr("ascontiguousarray")(
                obj, nb_::arg("dtype") = numpy_dtype<T>::name);
            auto nd = nb_::cast<nb_::ndarray<T, nb_::c_contig>>(arr);
            return contiguous_view<T>{std::move(arr), std::move(nd)};
        }
    } // namespace python
} // namespace alps
#endif // ALPS_PYTHON_NUMPY_COMPAT_HPP
