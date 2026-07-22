// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
/// Header-only template: dispatches a nb::object to a visitor based
/// on the object's Python type name. For numpy arrays the buffer
/// protocol (PEP 3118) provides the raw data + shape + format string;
/// for numpy scalars nb::cast handles the conversion through the
/// scalar's __int__/__float__/__complex__ methods. No dependence on
/// <numpy/arrayobject.h>.
#ifndef PYALPS_NGS_EXTRACT_FROM_PYOBJECT_HPP
#define PYALPS_NGS_EXTRACT_FROM_PYOBJECT_HPP
    #include <alps/ngs/cast.hpp>
    #include <alps/ngs/config.hpp>
    #include <nanobind/nanobind.h>
    #include <nanobind/ndarray.h>
    #include <nanobind/stl/complex.h>
    #include <nanobind/stl/string.h>
    #include <complex>
    #include <cstdint>
    #include <stdexcept>
    #include <string>
    #include <vector>
    namespace alps {
        namespace detail {
            namespace nb_ = nanobind;
            // True iff the ndarray view is C-contiguous: strides match
            // the canonical row-major layout (rightmost stride =
            // itemsize, each leftward stride = previous * shape[i+1]).
            // Zero-rank scalars are trivially contiguous.
            template <typename Arr>
            inline bool ndarray_is_c_contiguous(Arr const & arr) {
                if (arr.ndim() == 0) return true;
                int64_t expected = 1;
                for (int64_t i = arr.ndim() - 1; i >= 0; --i) {
                    if (arr.shape(i) == 0) return true; // empty array
                    if (arr.stride(i) != expected) return false;
                    expected *= arr.shape(i);
                }
                return true;
            }
            /// Dispatches `data` to `visitor` based on Python's type
            /// name. Visitor must be callable with bool / int / long /
            /// double / std::complex<double> / std::string /
            /// nb::list / nb::dict, plus the two-arg numpy form
            /// `visitor(T const*, std::vector<std::size_t>)` for each
            /// supported native numpy element type.
            template<typename T> void extract_from_pyobject_py11(T & visitor, nb_::handle data) {
                std::string dtype = data.ptr()->ob_type->tp_name;
                if (dtype == "bool") visitor(nb_::cast<bool>(data));
                else if (dtype == "int") visitor(nb_::cast<int>(data));
                else if (dtype == "long") visitor(nb_::cast<long>(data));
                else if (dtype == "float") visitor(nb_::cast<double>(data));
                else if (dtype == "complex") visitor(nb_::cast<std::complex<double>>(data));
                else if (dtype == "str") visitor(nb_::cast<std::string>(data));
                else if (dtype == "list") visitor(nb_::borrow<nb_::list>(data));
                else if (dtype == "tuple") {
                    // materialise the tuple as a list so the visitor only
                    // needs one sequence overload.
                    nb_::list as_list = nb_::steal<nb_::list>(
                        PySequence_List(data.ptr()));
                    visitor(as_list);
                }
                else if (dtype == "dict") visitor(nb_::borrow<nb_::dict>(data));
                // numpy scalars: extract through the scalar's own
                // __int__/__float__/__complex__ — no numpy C macros.
                else if (dtype == "numpy.str_" || dtype == "numpy.str")
                    visitor(std::string(nb_::cast<std::string>(nb_::str(data.attr("__str__")()))));
                else if (dtype == "numpy.bool_" || dtype == "numpy.bool")
                    visitor(nb_::cast<bool>(data));
                else if (dtype == "numpy.int8")  visitor(static_cast<std::int8_t>(nb_::cast<long>(data)));
                else if (dtype == "numpy.int16") visitor(static_cast<std::int16_t>(nb_::cast<long>(data)));
                else if (dtype == "numpy.int32") visitor(static_cast<std::int32_t>(nb_::cast<long>(data)));
                else if (dtype == "numpy.int64") visitor(static_cast<std::int64_t>(nb_::cast<long long>(data)));
                else if (dtype == "numpy.uint8")  visitor(static_cast<std::uint8_t>(nb_::cast<unsigned long>(data)));
                else if (dtype == "numpy.uint16") visitor(static_cast<std::uint16_t>(nb_::cast<unsigned long>(data)));
                else if (dtype == "numpy.uint32") visitor(static_cast<std::uint32_t>(nb_::cast<unsigned long>(data)));
                else if (dtype == "numpy.uint64") visitor(static_cast<std::uint64_t>(nb_::cast<unsigned long long>(data)));
                else if (dtype == "numpy.float32") visitor(static_cast<float>(nb_::cast<double>(data)));
                else if (dtype == "numpy.float64") visitor(nb_::cast<double>(data));
                else if (dtype == "numpy.complex64")
                    visitor(std::complex<float>(
                          nb_::cast<double>(data.attr("real").attr("__float__")())
                        , nb_::cast<double>(data.attr("imag").attr("__float__")())
                    ));
                else if (dtype == "numpy.complex128")
                    visitor(nb_::cast<std::complex<double>>(data));
                else if (dtype == "numpy.ndarray") {
                    // Raw buffer access via nb::ndarray, with a strict
                    // dtype match — nb::cast<nb::ndarray<T>>(arr) of a
                    // mismatched-dtype array silently coerces (e.g.
                    // int → bool yields all-true), so we inspect
                    // .dtype() ourselves and pick the matching arm.
                    // We require C-contiguity; the typical save path
                    // is bulk contiguous data and silently copying
                    // behind the user's back was the old
                    // PyArray_GETCONTIGUOUS behaviour we don't want
                    // to inherit.
                    auto arr_any = nb_::cast<nb_::ndarray<nb_::ro, nb_::c_contig>>(data);
                    std::vector<std::size_t> sizes;
                    sizes.reserve(arr_any.ndim());
                    for (std::size_t i = 0; i < arr_any.ndim(); ++i)
                        sizes.push_back(static_cast<std::size_t>(arr_any.shape(i)));
                    auto dt = arr_any.dtype();
                    #define DISPATCH_DTYPE(T)                                                           \
                        if (dt == nb_::dtype<T>())                                                      \
                            return visitor(static_cast<T const *>(arr_any.data()), sizes);
                    DISPATCH_DTYPE(bool)
                    DISPATCH_DTYPE(signed char)
                    DISPATCH_DTYPE(unsigned char)
                    DISPATCH_DTYPE(short)
                    DISPATCH_DTYPE(unsigned short)
                    DISPATCH_DTYPE(int)
                    DISPATCH_DTYPE(unsigned)
                    DISPATCH_DTYPE(long)
                    DISPATCH_DTYPE(unsigned long)
                    DISPATCH_DTYPE(long long)
                    DISPATCH_DTYPE(unsigned long long)
                    DISPATCH_DTYPE(float)
                    DISPATCH_DTYPE(double)
                    DISPATCH_DTYPE(std::complex<float>)
                    DISPATCH_DTYPE(std::complex<double>)
                    #undef DISPATCH_DTYPE
                    throw std::runtime_error(
                        "Unknown numpy element dtype at save site" + ALPS_STACKTRACE);
                } else
                    throw std::runtime_error("Unsupported type: " + dtype + ALPS_STACKTRACE);
            }
        } // namespace detail
    } // namespace alps
#endif  // PYALPS_NGS_EXTRACT_FROM_PYOBJECT_HPP
