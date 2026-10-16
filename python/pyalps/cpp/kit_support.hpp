// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Conversions shared by the models_c, ed_c and classic_dmrg_c extensions:
// dense local matrices and column-major eigenvector blocks as NumPy arrays,
// the kit's enumerations accepted by name as well as by value, and the
// row-major integer state tables of the exact-diagonalization bases.
#ifndef PYALPS_KIT_SUPPORT_HPP
#define PYALPS_KIT_SUPPORT_HPP

#include "numpy_compat.hpp"

#include <alps/models/lattice.hpp>
#include <alps/models/physical_basis.hpp>
#include <alps/models/qubit_encoding.hpp>
#include <alps/models/site.hpp>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pyalps::kit {

namespace nb = nanobind;
using Complex = std::complex<double>;

// -- NumPy interchange -------------------------------------------------------

// A row-major dim x dim matrix as a (dim, dim) complex128 array (a copy).
inline nb::object square_matrix(std::vector<Complex> const& values, std::size_t dim) {
    if (values.size() != dim * dim)
        throw std::logic_error("local matrix size does not match its dimension");
    return alps::python::make_numpy_array<Complex>(std::vector<Complex>(values), {dim, dim});
}

// Any square real or complex array (or nested sequence) as a row-major
// complex vector; `dim` receives the dimension.
inline std::vector<Complex> square_matrix_from(nb::handle array, std::size_t& dim, char const* what) {
    auto const view = alps::python::as_contiguous<Complex>(array);
    if (view.ndim() != 2 || view.shape(0) != view.shape(1))
        throw nb::type_error((std::string(what) + " must be a square two-dimensional array").c_str());
    dim = view.shape(0);
    return std::vector<Complex>(view.data(), view.data() + dim * dim);
}

// One-dimensional array (or sequence) of T.
template <typename T>
inline std::vector<T> vector_from(nb::handle array, char const* what) {
    auto const view = alps::python::as_contiguous<T>(array);
    if (view.ndim() != 1)
        throw nb::type_error((std::string(what) + " must be a one-dimensional array").c_str());
    return std::vector<T>(view.data(), view.data() + view.shape(0));
}

template <typename T>
inline nb::object vector_to_numpy(std::vector<T> values) {
    return alps::python::make_numpy_array<T>(std::move(values));
}

template <typename T>
inline nb::object matrix_to_numpy(std::vector<T> values, std::size_t rows, std::size_t columns) {
    if (values.size() != rows * columns)
        throw std::logic_error("matrix buffer does not match its shape");
    return alps::python::make_numpy_array<T>(std::move(values), {rows, columns});
}

// A column-major rows x columns block (eigenvectors[column * rows + row]) as a
// Fortran-ordered (rows, columns) array without a copy: column j of the array
// is eigenvector j.
template <typename T>
inline nb::object column_major(std::vector<T> values, std::size_t rows, std::size_t columns) {
    if (values.size() != rows * columns)
        throw std::logic_error("column-major buffer does not match its shape");
    if (values.empty())
        return alps::python::make_numpy_array<T>(std::move(values), {rows, columns});
    auto* owned = new std::vector<T>(std::move(values));
    nb::capsule owner(owned, [](void* pointer) noexcept { delete static_cast<std::vector<T>*>(pointer); });
    std::size_t const shape[2] = {rows, columns};
    std::int64_t const strides[2] = {1, static_cast<std::int64_t>(rows)};
    return nb::cast(nb::ndarray<nb::numpy, T>(owned->data(), 2, shape, owner, strides));
}

// A row-major (count, num_sites) integer table of local-state indices.
struct StateTable {
    std::vector<std::int64_t> states;
    std::size_t count = 0;
};

inline StateTable states_from(nb::handle array, std::size_t num_sites, char const* what) {
    auto const view = alps::python::as_contiguous<std::int64_t>(array);
    if (view.ndim() != 2)
        throw nb::type_error((std::string(what) + " must be a two-dimensional integer array").c_str());
    if (view.shape(1) != num_sites)
        throw nb::value_error((std::string(what) + " must have one column per site").c_str());
    StateTable table;
    table.count = view.shape(0);
    table.states.assign(view.data(), view.data() + table.count * num_sites);
    return table;
}

// -- enumerations by value or by name ------------------------------------------

template <typename E>
E enum_or_name(nb::handle value, std::initializer_list<std::pair<char const*, E>> names, char const* what) {
    if (nb::isinstance<E>(value)) return nb::cast<E>(value);
    if (nb::isinstance<nb::str>(value)) {
        std::string const text = nb::cast<std::string>(value);
        for (auto const& [name, item] : names)
            if (text == name) return item;
    }
    std::string accepted;
    for (auto const& [name, item] : names) accepted += (accepted.empty() ? "'" : ", '") + std::string(name) + "'";
    throw nb::type_error((std::string(what) + " must be one of " + accepted).c_str());
}

inline alps::lattice::Boundary boundary_of(nb::handle value) {
    using alps::lattice::Boundary;
    return enum_or_name<Boundary>(value, {{"periodic", Boundary::periodic}, {"open", Boundary::open}}, "boundary");
}

// A single boundary (uniform over the dimensions) or one per dimension.
inline std::vector<alps::lattice::Boundary> boundaries_of(nb::handle value, int dimension) {
    using alps::lattice::Boundary;
    if (nb::isinstance<nb::str>(value) || nb::isinstance<Boundary>(value))
        return std::vector<Boundary>(static_cast<std::size_t>(dimension), boundary_of(value));
    std::vector<Boundary> result;
    for (nb::handle item : nb::borrow<nb::iterable>(value)) result.push_back(boundary_of(item));
    if (static_cast<int>(result.size()) != dimension)
        throw nb::value_error("one boundary per dimension is needed");
    return result;
}

inline alps::models::FermionConvention convention_of(nb::handle value) {
    using alps::models::FermionConvention;
    return enum_or_name<FermionConvention>(value, {{"strict", FermionConvention::strict},
                                                   {"lattice", FermionConvention::lattice}}, "fermion convention");
}

inline alps::models::Statistics statistics_of(nb::handle value) {
    using alps::models::Statistics;
    return enum_or_name<Statistics>(value, {{"bosonic", Statistics::bosonic},
                                            {"fermionic", Statistics::fermionic}}, "statistics");
}

inline alps::models::Parity parity_of(nb::handle value) {
    using alps::models::Parity;
    return enum_or_name<Parity>(value, {{"even", Parity::even}, {"odd", Parity::odd}}, "parity");
}

inline alps::models::QubitEncodingScheme scheme_of(nb::handle value) {
    using alps::models::QubitEncodingScheme;
    return enum_or_name<QubitEncodingScheme>(value, {{"binary", QubitEncodingScheme::binary},
                                                     {"gray", QubitEncodingScheme::gray}}, "qubit encoding scheme");
}

inline std::string format_complex(Complex value) {
    char buffer[64];
    if (value.imag() == 0.) std::snprintf(buffer, sizeof(buffer), "%g", value.real());
    else std::snprintf(buffer, sizeof(buffer), "(%g%+gj)", value.real(), value.imag());
    return buffer;
}

} // namespace pyalps::kit

#endif // PYALPS_KIT_SUPPORT_HPP
