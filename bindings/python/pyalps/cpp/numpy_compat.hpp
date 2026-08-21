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
#include <type_traits>
#include <utility>
#include <vector>
namespace alps {
    namespace python {
        namespace nb_ = nanobind;
        // numpy dtype strings, indexed by the corresponding C++ type. Used by
        // as_contiguous() and by the integer arm of make_numpy_array, where
        // the dtype has to be named exactly (see the comment there).
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
        // Hand a heap-owned std::vector to NumPy through nb::ndarray. The
        // capsule keeps the vector alive for as long as the array (or any
        // view of it) exists, so the data is not copied on the way out.
        //
        // This used to allocate numpy.empty(shape, dtype=...) through the
        // cached numpy module and memcpy into it: two copies per conversion
        // (callers build a vector, which was then copied again) plus a Python
        // call. nb::ndarray hands NumPy the buffer C++ already owns.
        //
        // The returned array is writable: NumPy takes no ownership, only a
        // reference to the capsule, matching what numpy.empty + memcpy
        // produced. An empty shape yields a 0-d array, as numpy.empty(())
        // did.
        template <typename T>
        inline nb_::object make_numpy_array(std::vector<T> values,
                                            std::vector<std::size_t> const& shape) {
            if constexpr (std::is_integral_v<T>) {
                // Integers keep the numpy.empty + memcpy route. nb::ndarray
                // describes the buffer through DLPack, which encodes only
                // "signed 64-bit" and so cannot distinguish `long` from
                // `long long`. NumPy can: they are different dtype objects
                // with different .char ('l' vs 'q') and different reprs, and
                // routing int64 through DLPack changed
                // archive["/ints"] from int64 to longlong. Same bits, but a
                // visible change in returned dtype, which this migration has
                // no business making. Passing the dtype name keeps it exact.
                nb_::handle np = numpy_module();
                nb_::tuple shape_tuple = nb_::steal<nb_::tuple>(
                    PyTuple_New(static_cast<Py_ssize_t>(shape.size())));
                if (!shape_tuple.is_valid())
                    throw nb_::python_error();
                for (std::size_t i = 0; i < shape.size(); ++i) {
                    PyObject* dim = PyLong_FromUnsignedLongLong(shape[i]);
                    if (!dim)
                        throw nb_::python_error();
                    PyTuple_SetItem(shape_tuple.ptr(),
                                    static_cast<Py_ssize_t>(i), dim);
                }
                nb_::object array = np.attr("empty")(
                    shape_tuple, nb_::arg("dtype") = numpy_dtype<T>::name);
                if (!values.empty()) {
                    auto view = nb_::cast<nb_::ndarray<T, nb_::c_contig>>(array);
                    std::memcpy(view.data(), values.data(),
                                values.size() * sizeof(T));
                }
                return array;
            } else {
                // Floating point and complex have one unambiguous NumPy dtype
                // each, so the zero-copy route is exact.
                auto* owned = new std::vector<T>(std::move(values));
                nb_::capsule owner(owned, [](void* pointer) noexcept {
                    delete static_cast<std::vector<T>*>(pointer);
                });
                return nb_::cast(nb_::ndarray<nb_::numpy, T>(
                    owned->data(), shape.size(), shape.data(), owner));
            }
        }
        template <typename T>
        inline nb_::object make_numpy_array(std::vector<T> values) {
            std::size_t const size = values.size();
            return make_numpy_array<T>(std::move(values), {size});
        }
        // Copying entry point, for callers whose source is not a vector they
        // can give up (a single stack value, or a null pointer standing for an
        // empty extent). `data` may be null when the shape has no elements.
        template <typename T>
        inline nb_::object make_numpy_array(T const* data,
                                            std::vector<std::size_t> const& shape) {
            std::size_t total = 1;
            for (auto extent : shape) total *= extent;
            std::vector<T> values(total);
            if (total > 0 && data != nullptr)
                std::memcpy(values.data(), data, total * sizeof(T));
            return make_numpy_array<T>(std::move(values), shape);
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
