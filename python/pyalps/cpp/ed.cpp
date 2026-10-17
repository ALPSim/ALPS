// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// pyalps.ed: exact diagonalization of compiled operators (ALPS::ed). Product
// and projected problems, sparse assembly, the dense and Lanczos solvers,
// charge/momentum/point-group resolved spectra and thermodynamics. Matrices
// cross as NumPy arrays (eigenvectors Fortran-ordered, column j = state j)
// and canonical CSR values with a to_scipy() conversion; every solve releases
// the GIL.
//
// The surface follows src/alps/ed; the skilledwolf/alpsx bindings
// (python/alps/cpp/diag_bindings.cpp, projected_ed_bindings.hpp at
// d37e51227) exchanged dictionaries and were consulted for the GIL, ownership
// and column-major conventions only.
#include "kit_support.hpp"

#include <alps/ed.hpp>
#include <alps/models.hpp>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <complex>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace nb = nanobind;
namespace ed = alps::ed;
namespace mdl = alps::models;
using namespace pyalps::kit;

namespace {

// -- controls -----------------------------------------------------------------------

ed::EigenstateOptions eigen_options(int count, std::string const& method, std::string const& which,
                                    std::optional<double> tolerance, std::optional<int> max_iterations) {
    ed::EigenstateOptions options;
    options.count = count;
    if (method == "auto" || method == "automatic") options.method = ed::EigenstateMethod::automatic;
    else if (method == "dense") options.method = ed::EigenstateMethod::dense;
    else if (method == "lanczos") options.method = ed::EigenstateMethod::lanczos;
    else throw nb::value_error("method must be 'auto', 'dense' or 'lanczos'");
    options.which = which;
    options.tolerance = tolerance;
    options.max_iterations = max_iterations;
    return options;
}

char const* method_name(ed::EigenstateMethod method) {
    switch (method) {
    case ed::EigenstateMethod::automatic: return "auto";
    case ed::EigenstateMethod::dense: return "dense";
    case ed::EigenstateMethod::lanczos: return "lanczos";
    }
    return "auto";
}

#define PYALPS_EIGEN_ARGS                                                                                    \
    nb::arg("count") = 6, nb::arg("method") = "auto", nb::arg("which") = "SA", nb::arg("tolerance") = nb::none(), \
    nb::arg("max_iterations") = nb::none()

#define PYALPS_SPECTRUM_ARGS                                                                                 \
    nb::arg("count") = nb::none(), nb::arg("which") = "SA", nb::arg("tolerance") = nb::none(),               \
    nb::arg("max_iterations") = nb::none()

// -- sparse matrices ----------------------------------------------------------------

template <class Scalar>
nb::object scipy_csr(std::vector<Scalar> const& values, std::vector<std::int64_t> const& indices,
                     std::vector<std::int64_t> const& offsets, std::size_t rows, std::size_t columns) {
    nb::object sparse = nb::module_::import_("scipy.sparse");
    return sparse.attr("csr_matrix")(
        nb::make_tuple(vector_to_numpy(std::vector<Scalar>(values)), vector_to_numpy(std::vector<std::int64_t>(indices)),
                       vector_to_numpy(std::vector<std::int64_t>(offsets))),
        nb::arg("shape") = nb::make_tuple(rows, columns));
}

template <class Scalar>
void bind_sparse(nb::module_& m, char const* name, char const* doc) {
    using Matrix = ed::SparseMatrix<Scalar>;
    nb::class_<Matrix>(m, name, doc)
        .def(nb::init<>())
        .def_prop_ro("shape", [](Matrix const& s) { return std::make_pair(s.rows, s.columns); })
        .def_prop_ro("rows", [](Matrix const& s) { return s.rows; })
        .def_prop_ro("columns", [](Matrix const& s) { return s.columns; })
        .def_prop_ro("offsets", [](Matrix const& s) { return vector_to_numpy(std::vector<std::int64_t>(s.offsets)); },
                     "Row offsets (indptr), length rows + 1.")
        .def_prop_ro("indices", [](Matrix const& s) { return vector_to_numpy(std::vector<std::int64_t>(s.indices)); },
                     "Column indices, sorted within every row.")
        .def_prop_ro("values", [](Matrix const& s) { return vector_to_numpy(std::vector<Scalar>(s.values)); })
        .def_prop_ro("nnz", [](Matrix const& s) { return s.values.size(); })
        .def("to_scipy", [](Matrix const& s) { return scipy_csr(s.values, s.indices, s.offsets, s.rows, s.columns); },
             "The matrix as a scipy.sparse.csr_matrix (row = output state, column = input state).")
        .def("toarray", [](Matrix const& s) {
            return scipy_csr(s.values, s.indices, s.offsets, s.rows, s.columns).attr("toarray")();
        }, "The dense NumPy matrix.")
        .def("__repr__", [name](Matrix const& s) {
            return std::string(name) + "(" + std::to_string(s.rows) + "x" + std::to_string(s.columns) + ", nnz="
                 + std::to_string(s.values.size()) + ")";
        });
}

// COO contributions with the shape they belong to.
struct CooMatrix {
    ed::CooFromLocal coo;
    std::size_t rows = 0, columns = 0;
};

nb::object csr_of(CooMatrix const& coo, bool real) {
    if (real) return nb::cast(ed::canonical_csr<double>(coo.coo, coo.rows, coo.columns));
    return nb::cast(ed::canonical_csr<Complex>(coo.coo, coo.rows, coo.columns));
}

void bind_sparse_assembly(nb::module_& m) {
    bind_sparse<Complex>(m, "SparseMatrix", "Owned canonical complex CSR: sorted columns, summed duplicates, no exact zeros.");
    bind_sparse<double>(m, "RealSparseMatrix", "Owned canonical real CSR: sorted columns, summed duplicates, no exact zeros.");

    nb::class_<CooMatrix>(m, "CooMatrix", "Uncombined (row, column, value) contributions of an assembled operator.")
        .def_prop_ro("rows", [](CooMatrix const& c) { return vector_to_numpy(std::vector<std::int64_t>(c.coo.rows)); })
        .def_prop_ro("cols", [](CooMatrix const& c) { return vector_to_numpy(std::vector<std::int64_t>(c.coo.cols)); })
        .def_prop_ro("data", [](CooMatrix const& c) { return vector_to_numpy(std::vector<Complex>(c.coo.data)); })
        .def_prop_ro("is_complex", [](CooMatrix const& c) { return c.coo.is_complex; })
        .def_prop_ro("shape", [](CooMatrix const& c) { return std::make_pair(c.rows, c.columns); })
        .def_prop_ro("nnz", [](CooMatrix const& c) { return c.coo.data.size(); })
        .def("to_csr", [](CooMatrix const& c, bool real) { return csr_of(c, real); }, nb::arg("real") = false,
             "Canonical CSR; `real` requires every imaginary part to be exactly zero.")
        .def("to_scipy", [](CooMatrix const& c) {
            nb::object sparse = nb::module_::import_("scipy.sparse");
            return sparse.attr("coo_matrix")(
                nb::make_tuple(vector_to_numpy(std::vector<Complex>(c.coo.data)),
                               nb::make_tuple(vector_to_numpy(std::vector<std::int64_t>(c.coo.rows)),
                                              vector_to_numpy(std::vector<std::int64_t>(c.coo.cols)))),
                nb::arg("shape") = nb::make_tuple(c.rows, c.columns));
        }, "The contributions as a scipy.sparse.coo_matrix (duplicates sum on conversion).")
        .def("__repr__", [](CooMatrix const& c) {
            return "CooMatrix(" + std::to_string(c.rows) + "x" + std::to_string(c.columns) + ", nnz="
                 + std::to_string(c.coo.data.size()) + ")";
        });

    nb::class_<ed::ProductBasis>(m, "ProductBasis", "Enumerated product states, site 0 the most significant digit.")
        .def_prop_ro("dimension", [](ed::ProductBasis const& b) { return b.dimension; })
        .def_prop_ro("num_sites", [](ed::ProductBasis const& b) {
            return b.dimension ? b.states.size() / b.dimension : std::size_t{0};
        })
        .def_prop_ro("states", [](ed::ProductBasis const& b) {
            std::size_t const width = b.dimension ? b.states.size() / b.dimension : 0;
            return matrix_to_numpy(std::vector<std::int32_t>(b.states), b.dimension, width);
        }, "Row-major (dimension, num_sites) local-state indices in lexicographic order.")
        .def("__len__", [](ed::ProductBasis const& b) { return b.dimension; })
        .def("__repr__", [](ed::ProductBasis const& b) { return "ProductBasis(dimension=" + std::to_string(b.dimension) + ")"; });
    m.def("enumerate_product_basis", [](std::vector<int> const& site_dims,
                                        std::vector<std::pair<std::vector<std::vector<double>>, double>> const& constraints,
                                        std::size_t max_full_dimension) {
        std::vector<ed::BasisConstraint> native;
        for (auto const& [values, target] : constraints) native.push_back({values, target});
        nb::gil_scoped_release release;
        return ed::enumerate_product_basis(site_dims, native, max_full_dimension);
    }, nb::arg("site_dims"), nb::arg("constraints") = std::vector<std::pair<std::vector<std::vector<double>>, double>>{},
       nb::arg("max_full_dimension") = std::size_t{10'000'000},
       "The full product basis or one charge sector: `constraints` is a list of (values[site][state], target).");

    m.def("build_coo_from_local", [](mdl::CompiledOperator const& op, nb::handle states,
                                     std::vector<std::vector<double>> const& parity) {
        auto const table = states_from(states, op.site_dims.size(), "states");
        nb::gil_scoped_release release;
        return CooMatrix{ed::build_coo_from_local(op, table.states, table.count, parity), table.count, table.count};
    }, nb::arg("op"), nb::arg("states"), nb::arg("parity") = std::vector<std::vector<double>>{},
       "Assemble any operator on an explicit (count, num_sites) basis; `parity` gives one +/-1 diagonal per site "
       "for the Jordan-Wigner strings (all +1 when empty). Transitions leaving the basis are dropped.");
    m.def("build_csr_from_local", [](mdl::CompiledOperator const& op, nb::handle states,
                                     std::vector<std::vector<double>> const& parity, bool real) {
        auto const table = states_from(states, op.site_dims.size(), "states");
        CooMatrix coo;
        {
            nb::gil_scoped_release release;
            coo = CooMatrix{ed::build_coo_from_local(op, table.states, table.count, parity), table.count, table.count};
        }
        return csr_of(coo, real);
    }, nb::arg("op"), nb::arg("states"), nb::arg("parity") = std::vector<std::vector<double>>{}, nb::arg("real") = false,
       "build_coo_from_local reduced to canonical CSR.");

    nb::class_<ed::SparseFromLocal>(m, "SparseFromLocal",
        "Real CSR of a CompiledModel on its enumerated (least-significant-first) basis, with the state codes.")
        .def_prop_ro("dim", [](ed::SparseFromLocal const& s) { return s.dim; })
        .def_prop_ro("shape", [](ed::SparseFromLocal const& s) { return std::make_pair(s.dim, s.dim); })
        .def_prop_ro("offsets", [](ed::SparseFromLocal const& s) { return vector_to_numpy(std::vector<std::int64_t>(s.indptr)); })
        .def_prop_ro("indices", [](ed::SparseFromLocal const& s) { return vector_to_numpy(std::vector<std::int64_t>(s.indices)); })
        .def_prop_ro("values", [](ed::SparseFromLocal const& s) { return vector_to_numpy(std::vector<double>(s.data)); })
        .def_prop_ro("states", [](ed::SparseFromLocal const& s) { return vector_to_numpy(std::vector<std::uint64_t>(s.states)); },
                     "Mixed-radix codes of the basis states, site 0 least significant, ascending.")
        .def("to_scipy", [](ed::SparseFromLocal const& s) {
            return scipy_csr(s.data, s.indices, s.indptr, static_cast<std::size_t>(s.dim), static_cast<std::size_t>(s.dim));
        })
        .def("__repr__", [](ed::SparseFromLocal const& s) {
            return "SparseFromLocal(dim=" + std::to_string(s.dim) + ", nnz=" + std::to_string(s.data.size()) + ")";
        });
    m.def("build_sparse_from_local", [](mdl::CompiledModel const& model, std::map<int, std::vector<double>> const& parity_by_type,
                                        std::vector<std::pair<std::string, double>> const& constraints,
                                        std::uint64_t max_full_dim) {
        std::vector<ed::QnConstraint> native;
        for (auto const& [name, target] : constraints) native.push_back({name, target});
        nb::gil_scoped_release release;
        return ed::build_sparse_from_local(model, parity_by_type, native, max_full_dim);
    }, nb::arg("model"), nb::arg("parity_by_type") = std::map<int, std::vector<double>>{},
       nb::arg("constraints") = std::vector<std::pair<std::string, double>>{},
       nb::arg("max_full_dim") = std::uint64_t{1} << 24,
       "Real CSR of a compiled model, optionally in the sector where the named diag observables sum to targets.");
}

// -- eigensolvers --------------------------------------------------------------------

struct DenseSolution {
    std::size_t dimension = 0;
    std::vector<double> eigenvalues;
    std::variant<std::vector<double>, std::vector<Complex>> eigenvectors;
    bool has_vectors = false;
};

struct LanczosSolution {
    std::size_t dimension = 0;
    std::vector<double> eigenvalues;
    std::variant<std::vector<double>, std::vector<Complex>> eigenvectors;
    bool has_vectors = false;
    std::vector<double> residual_norms;
    int nconv = 0;
    double tolerance = 0.;
    int max_iterations = 0, iterations = 0, runs = 0;
};

nb::object vectors_of(std::variant<std::vector<double>, std::vector<Complex>> const& buffer, std::size_t rows,
                      std::size_t columns) {
    return std::visit([&](auto const& values) { return column_major(values, rows, columns); }, buffer);
}

bool is_complex_object(nb::handle array) {
    return nb::cast<bool>(alps::python::numpy_module().attr("iscomplexobj")(array));
}

void bind_eigensolvers(nb::module_& m) {
    nb::class_<DenseSolution>(m, "DenseEigenResult", "Ascending eigenvalues and, when requested, eigenvectors as columns.")
        .def_prop_ro("eigenvalues", [](DenseSolution const& s) { return vector_to_numpy(std::vector<double>(s.eigenvalues)); })
        .def_prop_ro("eigenvectors", [](DenseSolution const& s) -> nb::object {
            if (!s.has_vectors) return nb::none();
            return vectors_of(s.eigenvectors, s.dimension, s.eigenvalues.size());
        }, "(dimension, count) array, column j the eigenvector of eigenvalue j; None when not requested.")
        .def("__len__", [](DenseSolution const& s) { return s.eigenvalues.size(); })
        .def("__repr__", [](DenseSolution const& s) {
            return "DenseEigenResult(" + std::to_string(s.eigenvalues.size()) + " of " + std::to_string(s.dimension) + ")";
        });
    m.def("dense_eigensystem", [](nb::handle matrix, bool eigenvectors, std::optional<std::size_t> count,
                                  std::string const& which) {
        nb::object array = alps::python::numpy_module().attr("asarray")(matrix);
        bool const complex_input = is_complex_object(array);
        DenseSolution solution;
        auto solve = [&](auto tag) {
            using Scalar = decltype(tag);
            auto const view = alps::python::as_contiguous<Scalar>(array);
            if (view.ndim() != 2 || view.shape(0) != view.shape(1))
                throw nb::type_error("the matrix must be a square two-dimensional array");
            std::size_t const n = view.shape(0);
            std::vector<Scalar> packed(n * n);
            for (std::size_t row = 0; row < n; ++row)
                for (std::size_t column = 0; column < n; ++column) packed[column * n + row] = view.data()[row * n + column];
            nb::gil_scoped_release release;
            auto result = ed::dense_eigensystem(std::move(packed), n, eigenvectors, count, which);
            solution.dimension = n;
            solution.eigenvalues = std::move(result.eigenvalues);
            solution.eigenvectors = std::move(result.eigenvectors_column_major);
            solution.has_vectors = eigenvectors;
        };
        if (complex_input) solve(Complex{});
        else solve(double{});
        return solution;
    }, nb::arg("matrix"), nb::arg("eigenvectors") = true, nb::arg("count") = nb::none(), nb::arg("which") = "SA",
       "LAPACK diagonalization of the Hermitian operator given by the matrix's lower triangle. "
       "`count` selects that many states by `which` (SA, LA, SM, LM); omitted returns the full spectrum.");
    m.def("select_eigenvalues", &ed::select_eigenvalues, nb::arg("values"), nb::arg("count"), nb::arg("which"),
          "Indices of the requested SA/LA/SM/LM subset, in ascending energy order.");
    m.def("validate_eigensolver_controls", &ed::validate_eigensolver_controls, nb::arg("which"),
          nb::arg("tolerance") = nb::none(), nb::arg("max_iterations") = nb::none());

    nb::class_<LanczosSolution>(m, "LanczosResult", "A deflated Lanczos solve: eigenpairs with multiplicity and the controls in effect.")
        .def_prop_ro("eigenvalues", [](LanczosSolution const& s) { return vector_to_numpy(std::vector<double>(s.eigenvalues)); })
        .def_prop_ro("eigenvectors", [](LanczosSolution const& s) -> nb::object {
            if (!s.has_vectors) return nb::none();
            return vectors_of(s.eigenvectors, s.dimension, s.eigenvalues.size());
        })
        .def_prop_ro("residual_norms", [](LanczosSolution const& s) { return vector_to_numpy(std::vector<double>(s.residual_norms)); })
        .def_ro("nconv", &LanczosSolution::nconv)
        .def_ro("tolerance", &LanczosSolution::tolerance)
        .def_ro("max_iterations", &LanczosSolution::max_iterations)
        .def_ro("iterations", &LanczosSolution::iterations)
        .def_ro("runs", &LanczosSolution::runs)
        .def("__len__", [](LanczosSolution const& s) { return s.eigenvalues.size(); })
        .def("__repr__", [](LanczosSolution const& s) {
            return "LanczosResult(" + std::to_string(s.eigenvalues.size()) + " eigenvalues, " + std::to_string(s.iterations)
                 + " steps in " + std::to_string(s.runs) + " runs)";
        });
    m.def("lanczos_csr", [](nb::handle offsets, nb::handle indices, nb::handle values, std::int64_t dimension, int count,
                            std::string const& which, double tolerance, int max_iterations, bool want_vectors) {
        auto const native_offsets = vector_from<std::int64_t>(offsets, "offsets");
        auto const native_indices = vector_from<std::int64_t>(indices, "indices");
        LanczosSolution solution;
        auto solve = [&](auto tag) {
            using Scalar = decltype(tag);
            auto const native_values = vector_from<Scalar>(values, "values");
            nb::gil_scoped_release release;
            auto result = ed::lanczos_csr(native_offsets, native_indices, native_values, dimension, count, which, tolerance,
                                          max_iterations, want_vectors);
            solution.dimension = static_cast<std::size_t>(dimension);
            solution.eigenvalues = std::move(result.eigenvalues);
            solution.eigenvectors = std::move(result.eigenvectors);
            solution.has_vectors = want_vectors;
            solution.residual_norms = std::move(result.residual_norms);
            solution.nconv = result.nconv;
            solution.tolerance = result.tolerance;
            solution.max_iterations = result.max_iterations;
            solution.iterations = result.iterations;
            solution.runs = result.runs;
        };
        if (is_complex_object(values)) solve(Complex{});
        else solve(double{});
        return solution;
    }, nb::arg("offsets"), nb::arg("indices"), nb::arg("values"), nb::arg("dimension"), nb::arg("count"),
       nb::arg("which") = "SA", nb::arg("tolerance") = 0., nb::arg("max_iterations") = 0, nb::arg("want_vectors") = true,
       "Deflated IETL Lanczos on a real symmetric or complex Hermitian CSR matrix (SA, LA or LM; "
       "tolerance <= 0 and max_iterations <= 0 select the driver defaults). Degenerate levels come back with multiplicity.");
}

// -- product problems ------------------------------------------------------------------

nb::object product_hamiltonian(ed::ProductProblem const& p) {
    return std::visit([&](auto const& values) {
        using Scalar = typename std::decay_t<decltype(values)>::value_type;
        ed::SparseMatrix<Scalar> matrix{p.dimension(), p.dimension(), p.row_offsets(), p.column_indices(), values};
        return nb::cast(std::move(matrix));
    }, p.matrix_values());
}

void bind_product(nb::module_& m) {
    nb::class_<ed::ProductProblem>(m, "ProductProblem",
        "A prepared Hermitian problem on a product basis (full space or one charge sector): the assembled "
        "canonical CSR, solved densely or by Lanczos.")
        .def_prop_ro("dimension", &ed::ProductProblem::dimension)
        .def_prop_ro("site_dims", &ed::ProductProblem::site_dims)
        .def_prop_ro("physical_basis", &ed::ProductProblem::physical_basis)
        .def_prop_ro("sector", &ed::ProductProblem::sector)
        .def_prop_ro("states", [](ed::ProductProblem const& p) {
            return matrix_to_numpy(std::vector<std::int32_t>(p.states()), p.dimension(), p.site_dims().size());
        }, "Row-major (dimension, num_sites) product states, site 0 most significant.")
        .def_prop_ro("hamiltonian", &product_hamiltonian, "The assembled matrix (real or complex canonical CSR).")
        .def("eigenstates", [](ed::ProductProblem const& p, int count, std::string const& method, std::string const& which,
                               std::optional<double> tolerance, std::optional<int> max_iterations) {
            auto const options = eigen_options(count, method, which, tolerance, max_iterations);
            nb::gil_scoped_release release;
            return p.eigenstates(options);
        }, PYALPS_EIGEN_ARGS,
           "The lowest (`which`) `count` eigenpairs; 'auto' uses Lanczos above dimension 1000.")
        .def("operator_matrix", [](ed::ProductProblem const& p, mdl::CompiledOperator const& op) {
            nb::gil_scoped_release release;
            return p.operator_matrix(op);
        }, nb::arg("op"), "Any operator on the problem's local spaces, assembled on its basis (projected to the sector).")
        .def("__repr__", [](ed::ProductProblem const& p) { return "ProductProblem(dimension=" + std::to_string(p.dimension()) + ")"; });

    nb::class_<ed::ProductEigenResult>(m, "ProductEigenResult", "Eigenpairs of a ProductProblem.")
        .def_ro("problem", &ed::ProductEigenResult::problem)
        .def_prop_ro("eigenvalues", [](ed::ProductEigenResult const& r) { return vector_to_numpy(std::vector<double>(r.eigenvalues)); })
        .def_prop_ro("eigenvectors", [](ed::ProductEigenResult const& r) {
            return vectors_of(r.eigenvectors, r.problem.dimension(), r.eigenvalues.size());
        }, "(dimension, count) array in the problem basis; column j belongs to eigenvalue j.")
        .def_prop_ro("method", [](ed::ProductEigenResult const& r) { return method_name(r.method); })
        .def_ro("sparse_requested", &ed::ProductEigenResult::sparse_requested)
        .def_ro("tolerance", &ed::ProductEigenResult::tolerance)
        .def_ro("max_iterations", &ed::ProductEigenResult::max_iterations)
        .def_ro("iterations", &ed::ProductEigenResult::iterations)
        .def_prop_ro("residual_norms", [](ed::ProductEigenResult const& r) { return vector_to_numpy(std::vector<double>(r.residual_norms)); })
        .def("expectation_values", [](ed::ProductEigenResult const& r, mdl::CompiledOperator const& op) {
            std::vector<Complex> values;
            {
                nb::gil_scoped_release release;
                values = ed::expectation_values(r, op);
            }
            return vector_to_numpy(std::move(values));
        }, nb::arg("op"), "<psi_i|O|psi_i> for every eigenvector, O assembled on the problem's basis.")
        .def("__len__", [](ed::ProductEigenResult const& r) { return r.eigenvalues.size(); })
        .def("__repr__", [](ed::ProductEigenResult const& r) {
            return "ProductEigenResult(" + std::to_string(r.eigenvalues.size()) + " states, " + method_name(r.method) + ")";
        });
    m.def("expectation_values", [](ed::ProductEigenResult const& r, mdl::CompiledOperator const& op) {
        std::vector<Complex> values;
        {
            nb::gil_scoped_release release;
            values = ed::expectation_values(r, op);
        }
        return vector_to_numpy(std::move(values));
    }, nb::arg("result"), nb::arg("op"), "<psi_i|O|psi_i> for every eigenvector of a product result.");

    m.def("prepare_product", [](mdl::Hamiltonian const& model, std::map<std::string, double> sector, std::size_t max_full_dimension) {
        ed::ProductOptions options;
        options.sector = std::move(sector);
        options.max_full_dimension = max_full_dimension;
        nb::gil_scoped_release release;
        return ed::prepare_product(model, options);
    }, nb::arg("model"), nb::arg("sector") = std::map<std::string, double>{}, nb::arg("max_full_dimension") = std::size_t{10'000'000},
       "Validate Hermiticity, compile the model and assemble it on the full product basis or the charge sector "
       "`sector` (name -> total of a declared charge).");
    m.def("prepare_product", [](mdl::TermSum const& terms, std::vector<mdl::SiteBasisPtr> const& sites,
                                std::vector<std::string> const& conserved, nb::handle convention,
                                std::map<std::string, double> sector, std::size_t max_full_dimension) {
        auto const conv = convention_of(convention);
        ed::ProductOptions options;
        options.sector = std::move(sector);
        options.max_full_dimension = max_full_dimension;
        nb::gil_scoped_release release;
        return ed::prepare_product(terms, sites, conserved, conv, options);
    }, nb::arg("terms"), nb::arg("sites"), nb::arg("conserved_quantum_numbers") = std::vector<std::string>{},
       nb::arg("fermion_convention") = "strict", nb::arg("sector") = std::map<std::string, double>{},
       nb::arg("max_full_dimension") = std::size_t{10'000'000},
       "The same for located terms on ordered sites.");
    m.def("prepare_product", [](mdl::CompiledOperator const& op, mdl::PhysicalBasis const& physical,
                                std::map<std::string, double> sector, std::size_t max_full_dimension) {
        ed::ProductOptions options;
        options.sector = std::move(sector);
        options.max_full_dimension = max_full_dimension;
        nb::gil_scoped_release release;
        return ed::prepare_product(op, physical, options);
    }, nb::arg("op"), nb::arg("physical_basis"), nb::arg("sector") = std::map<std::string, double>{},
       nb::arg("max_full_dimension") = std::size_t{10'000'000},
       "The same for a compiled operator (Jordan-Wigner strings included) on the local spaces of a PhysicalBasis.");
}

// -- symmetry spectra, finite groups, projected bases ------------------------------------

std::string block_repr(ed::SpectrumBlock const& b) {
    std::string out = "SpectrumBlock(dimension=" + std::to_string(b.block_dimension) + ", eigenvalues="
                    + std::to_string(b.eigenvalues.size());
    if (!b.momentum.empty()) {
        out += ", momentum=[";
        for (std::size_t i = 0; i < b.momentum.size(); ++i) {
            char buffer[32];
            std::snprintf(buffer, sizeof(buffer), "%g", b.momentum[i]);
            out += (i ? ", " : "") + std::string(buffer);
        }
        out += "]";
    }
    if (!b.quantum_numbers.empty()) {
        out += ", charges=[";
        for (std::size_t i = 0; i < b.quantum_numbers.size(); ++i) {
            char buffer[32];
            std::snprintf(buffer, sizeof(buffer), "%g", b.quantum_numbers[i]);
            out += (i ? ", " : "") + std::string(buffer);
        }
        out += "]";
    }
    if (!b.irrep.empty()) out += ", irrep='" + b.irrep + "'";
    return out + ", " + b.method + ")";
}

void bind_symmetry(nb::module_& m) {
    nb::class_<ed::SpectrumLocalSpace>(m, "SpectrumLocalSpace",
        "One site's parity signs and charge rows (one row per conserved charge, in the caller's order).")
        .def(nb::init<>())
        .def("__init__", [](ed::SpectrumLocalSpace* self, std::vector<double> parity, std::vector<std::vector<double>> quantum_numbers) {
            new (self) ed::SpectrumLocalSpace{std::move(parity), std::move(quantum_numbers)};
        }, nb::arg("parity"), nb::arg("quantum_numbers") = std::vector<std::vector<double>>{})
        .def_rw("parity", &ed::SpectrumLocalSpace::parity)
        .def_rw("quantum_numbers", &ed::SpectrumLocalSpace::quantum_numbers)
        .def_prop_ro("dim", [](ed::SpectrumLocalSpace const& s) { return s.parity.size(); });
    nb::class_<ed::SpectrumTranslation>(m, "SpectrumTranslation",
        "A primitive translation mapping source site -> translated site, with its period.")
        .def("__init__", [](ed::SpectrumTranslation* self, std::vector<std::size_t> sites, std::size_t period) {
            new (self) ed::SpectrumTranslation{std::move(sites), period};
        }, nb::arg("sites"), nb::arg("period"))
        .def_rw("sites", &ed::SpectrumTranslation::sites)
        .def_rw("period", &ed::SpectrumTranslation::period);
    nb::class_<ed::SpectrumBlock>(m, "SpectrumBlock",
        "One joint charge/momentum (or irrep) block: its dimension, labels, ascending eigenvalues and solver report.")
        .def_ro("block_dimension", &ed::SpectrumBlock::block_dimension)
        .def_ro("momentum", &ed::SpectrumBlock::momentum)
        .def_ro("quantum_numbers", &ed::SpectrumBlock::quantum_numbers)
        .def_prop_ro("eigenvalues", [](ed::SpectrumBlock const& b) { return vector_to_numpy(std::vector<double>(b.eigenvalues)); })
        .def_ro("method", &ed::SpectrumBlock::method)
        .def_ro("tolerance", &ed::SpectrumBlock::tolerance)
        .def_ro("max_iterations", &ed::SpectrumBlock::max_iterations)
        .def_ro("iterations", &ed::SpectrumBlock::iterations)
        .def_ro("nconv", &ed::SpectrumBlock::nconv)
        .def_ro("irrep", &ed::SpectrumBlock::irrep)
        .def_ro("irrep_dimension", &ed::SpectrumBlock::irrep_dimension)
        .def_ro("group_order", &ed::SpectrumBlock::group_order)
        .def("__len__", [](ed::SpectrumBlock const& b) { return b.eigenvalues.size(); })
        .def("__repr__", &block_repr);
    m.def("symmetry_spectrum", [](nb::handle offsets, nb::handle indices, nb::handle values, std::size_t dimension,
                                  nb::handle product_states, std::vector<ed::SpectrumLocalSpace> const& local_spaces,
                                  std::vector<ed::SpectrumTranslation> const& translations, int count,
                                  std::vector<std::size_t> const& momentum_indices, std::string const& which,
                                  std::optional<double> tolerance, std::optional<int> max_iterations) {
        auto const native_offsets = vector_from<std::int64_t>(offsets, "offsets");
        auto const native_indices = vector_from<std::int64_t>(indices, "indices");
        auto const native_values = vector_from<Complex>(values, "values");
        auto const table = states_from(product_states, local_spaces.size(), "product_states");
        nb::gil_scoped_release release;
        return ed::symmetry_spectrum(native_offsets, native_indices, native_values, dimension, table.states, local_spaces,
                                     translations, count, momentum_indices, which, tolerance, max_iterations);
    }, nb::arg("offsets"), nb::arg("indices"), nb::arg("values"), nb::arg("dimension"), nb::arg("product_states"),
       nb::arg("local_spaces"), nb::arg("translations") = std::vector<ed::SpectrumTranslation>{}, nb::arg("count") = 6,
       nb::arg("momentum_indices") = std::vector<std::size_t>{}, nb::arg("which") = "SA", nb::arg("tolerance") = nb::none(),
       nb::arg("max_iterations") = nb::none(),
       "Resolve a Hermitian product-basis CSR matrix into joint charge and momentum blocks.");

    nb::class_<ed::FiniteGroupIrrep>(m, "FiniteGroupIrrep",
        "A finite site-permutation group closed from its generators, with one unitary irrep given by the "
        "generators' row-major matrices (checked for the relations, unitarity and irreducibility).")
        .def("__init__", [](ed::FiniteGroupIrrep* self, std::vector<std::vector<std::size_t>> generators,
                            std::vector<nb::object> const& generator_matrices, std::size_t dimension, std::string label,
                            std::size_t max_order) {
            std::vector<std::vector<Complex>> matrices;
            for (auto const& matrix : generator_matrices) {
                std::size_t dim = 0;
                matrices.push_back(square_matrix_from(matrix, dim, "generator matrix"));
                if (dim != dimension) throw nb::value_error("generator matrices must have the irrep dimension");
            }
            nb::gil_scoped_release release;
            new (self) ed::FiniteGroupIrrep(std::move(generators), std::move(matrices), dimension, std::move(label), max_order);
        }, nb::arg("generators"), nb::arg("generator_matrices"), nb::arg("dimension"), nb::arg("label"),
           nb::arg("max_order") = std::size_t{4096})
        .def_prop_ro("order", &ed::FiniteGroupIrrep::order)
        .def_prop_ro("num_sites", &ed::FiniteGroupIrrep::num_sites)
        .def_prop_ro("dimension", &ed::FiniteGroupIrrep::dimension)
        .def_prop_ro("label", &ed::FiniteGroupIrrep::label)
        .def_prop_ro("generators", &ed::FiniteGroupIrrep::generators)
        .def_prop_ro("permutations", &ed::FiniteGroupIrrep::permutations, "Every element as a site permutation, identity first.")
        .def_prop_ro("matrices", [](ed::FiniteGroupIrrep const& g) {
            std::vector<nb::object> out;
            for (auto const& matrix : g.matrices()) out.push_back(square_matrix(matrix, g.dimension()));
            return out;
        }, "The irrep matrix of every element, in permutation order.")
        .def_prop_ro("characters", [](ed::FiniteGroupIrrep const& g) { return vector_to_numpy(std::vector<Complex>(g.characters())); })
        .def("__repr__", [](ed::FiniteGroupIrrep const& g) {
            return "FiniteGroupIrrep('" + g.label() + "', order=" + std::to_string(g.order()) + ", dimension="
                 + std::to_string(g.dimension()) + ")";
        });
    m.def("dihedral_irrep", &ed::dihedral_irrep, nb::arg("n"), nb::arg("angular_index") = 0, nb::arg("reflection") = 1,
          "D_n on a ring of n sites: rotation i -> i+1, reflection i -> -i; k = angular_index in [0, n/2].");
    m.def("validate_irrep_family", &ed::validate_irrep_family, nb::arg("irreps"),
          "Require a complete family of inequivalent irreps of one group.");

    nb::class_<ed::ProjectedBasis>(m, "ProjectedBasis",
        "An orthonormal symmetry basis: orbit representatives of a charge sector under commuting translations "
        "at one momentum, or the isotypic component of a finite-group irrep.")
        .def_prop_ro("dimension", &ed::ProjectedBasis::dimension)
        .def_prop_ro("num_sites", &ed::ProjectedBasis::num_sites)
        .def_prop_ro("parent_dimension", &ed::ProjectedBasis::parent_dimension)
        .def_prop_ro("momentum_index", &ed::ProjectedBasis::momentum_index)
        .def_prop_ro("site_dims", &ed::ProjectedBasis::site_dims)
        .def_prop_ro("local_spaces", &ed::ProjectedBasis::local_spaces)
        .def_prop_ro("translations", &ed::ProjectedBasis::translations)
        .def_prop_ro("charge_targets", &ed::ProjectedBasis::charge_targets)
        .def_prop_ro("momentum", &ed::ProjectedBasis::momentum)
        .def_prop_ro("representatives", [](ed::ProjectedBasis const& b) {
            return matrix_to_numpy(std::vector<std::int32_t>(b.representatives()), b.dimension(), b.num_sites());
        }, "Row-major (dimension, num_sites) orbit representatives.")
        .def_prop_ro("orbit_sizes", &ed::ProjectedBasis::orbit_sizes)
        .def_prop_ro("finite_group", [](ed::ProjectedBasis const& b) -> nb::object {
            if (!b.finite_group()) return nb::none();
            return nb::cast(*b.finite_group());
        }, "The irrep of a finite-group basis, None for a translation basis.")
        .def("__len__", &ed::ProjectedBasis::dimension)
        .def("__repr__", [](ed::ProjectedBasis const& b) {
            return "ProjectedBasis(dimension=" + std::to_string(b.dimension()) + ", parent=" + std::to_string(b.parent_dimension()) + ")";
        });
    m.def("make_projected_basis", [](std::vector<ed::SpectrumLocalSpace> const& local_spaces,
                                     std::vector<ed::SpectrumTranslation> const& translations,
                                     std::vector<std::optional<double>> const& charge_targets, std::size_t momentum_index) {
        nb::gil_scoped_release release;
        return ed::make_projected_basis(local_spaces, translations, charge_targets, momentum_index);
    }, nb::arg("local_spaces"), nb::arg("translations") = std::vector<ed::SpectrumTranslation>{},
       nb::arg("charge_targets") = std::vector<std::optional<double>>{}, nb::arg("momentum_index") = 0,
       "Empty targets leave every charge unrestricted; otherwise one target or None per charge row.");
    m.def("make_finite_group_basis", [](std::vector<ed::SpectrumLocalSpace> const& local_spaces, ed::FiniteGroupIrrep const& irrep,
                                        std::vector<std::optional<double>> const& charge_targets) {
        nb::gil_scoped_release release;
        return ed::make_finite_group_basis(local_spaces, irrep, charge_targets);
    }, nb::arg("local_spaces"), nb::arg("irrep"), nb::arg("charge_targets") = std::vector<std::optional<double>>{});
    m.def("projected_operator", [](mdl::CompiledOperator const& op, ed::ProjectedBasis const& source, ed::ProjectedBasis const& target) {
        nb::gil_scoped_release release;
        return CooMatrix{ed::projected_operator(op, source, target), target.dimension(), source.dimension()};
    }, nb::arg("op"), nb::arg("source"), nb::arg("target"), "P_target^dagger O P_source as COO contributions.");
    m.def("projected_sparse_operator", [](mdl::CompiledOperator const& op, ed::ProjectedBasis const& source, ed::ProjectedBasis const& target) {
        nb::gil_scoped_release release;
        return ed::projected_sparse_operator<Complex>(op, source, target);
    }, nb::arg("op"), nb::arg("source"), nb::arg("target"), "P_target^dagger O P_source as canonical CSR.");
    m.def("apply_projected_operator", [](mdl::CompiledOperator const& op, ed::ProjectedBasis const& source,
                                         ed::ProjectedBasis const& target, nb::handle state) {
        auto const vector = vector_from<Complex>(state, "state");
        std::vector<Complex> result;
        {
            nb::gil_scoped_release release;
            result = ed::apply_projected_operator(op, source, target, vector);
        }
        return vector_to_numpy(std::move(result));
    }, nb::arg("op"), nb::arg("source"), nb::arg("target"), nb::arg("state"),
       "O applied to a state of the source basis, returned in the target basis.");
    m.def("expect_projected_operator", [](mdl::CompiledOperator const& op, ed::ProjectedBasis const& basis, nb::handle state) {
        auto const vector = vector_from<Complex>(state, "state");
        nb::gil_scoped_release release;
        return ed::expect_projected_operator(op, basis, vector);
    }, nb::arg("op"), nb::arg("basis"), nb::arg("state"), "<state|O|state> in a projected basis.");
    nb::class_<ed::ProjectedExpansion>(m, "ProjectedExpansion", "A projected state expanded into its parent product sector.")
        .def_ro("basis", &ed::ProjectedExpansion::basis)
        .def_prop_ro("amplitudes", [](ed::ProjectedExpansion const& e) { return vector_to_numpy(std::vector<Complex>(e.amplitudes)); });
    m.def("expand_projected_state", [](ed::ProjectedBasis const& basis, nb::handle state) {
        auto const vector = vector_from<Complex>(state, "state");
        nb::gil_scoped_release release;
        return ed::expand_projected_state(basis, vector);
    }, nb::arg("basis"), nb::arg("state"));
    m.def("symmetry_spectrum_from_local", [](mdl::CompiledOperator const& op, std::vector<ed::SpectrumLocalSpace> const& local_spaces,
                                             std::vector<ed::SpectrumTranslation> const& translations,
                                             std::vector<std::optional<double>> const& charge_targets, int count,
                                             std::vector<std::size_t> const& momentum_indices, std::string const& which,
                                             std::optional<double> tolerance, std::optional<int> max_iterations) {
        nb::gil_scoped_release release;
        return ed::symmetry_spectrum_from_local(op, local_spaces, translations, charge_targets, count, momentum_indices, which,
                                                tolerance, max_iterations);
    }, nb::arg("op"), nb::arg("local_spaces"), nb::arg("translations") = std::vector<ed::SpectrumTranslation>{},
       nb::arg("charge_targets") = std::vector<std::optional<double>>{}, nb::arg("count") = 6,
       nb::arg("momentum_indices") = std::vector<std::size_t>{}, nb::arg("which") = "SA", nb::arg("tolerance") = nb::none(),
       nb::arg("max_iterations") = nb::none(), "Joint charge/momentum spectra without a parent matrix.");
    m.def("finite_spectrum_from_local", [](mdl::CompiledOperator const& op, std::vector<ed::SpectrumLocalSpace> const& local_spaces,
                                           std::vector<ed::FiniteGroupIrrep> const& irreps,
                                           std::vector<std::optional<double>> const& charge_targets, int count,
                                           std::string const& which, std::optional<double> tolerance, std::optional<int> max_iterations) {
        nb::gil_scoped_release release;
        return ed::finite_spectrum_from_local(op, local_spaces, irreps, charge_targets, count, which, tolerance, max_iterations);
    }, nb::arg("op"), nb::arg("local_spaces"), nb::arg("irreps"), nb::arg("charge_targets") = std::vector<std::optional<double>>{},
       nb::arg("count") = 6, nb::arg("which") = "SA", nb::arg("tolerance") = nb::none(), nb::arg("max_iterations") = nb::none());
}

// -- projected problems and model spectra ---------------------------------------------------

ed::ModelSpectrumOptions model_options(std::map<std::string, double> sector, std::vector<double> momentum, bool all_momenta,
                                       std::optional<int> count, std::string const& which, std::optional<double> tolerance,
                                       std::optional<int> max_iterations) {
    ed::ModelSpectrumOptions options;
    options.sector = std::move(sector);
    options.momentum = std::move(momentum);
    options.all_momenta = all_momenta;
    options.count = count;
    options.which = which;
    options.tolerance = tolerance;
    options.max_iterations = max_iterations;
    return options;
}

ed::FiniteSpectrumOptions finite_options(std::map<std::string, double> sector, std::optional<int> count, std::string const& which,
                                         std::optional<double> tolerance, std::optional<int> max_iterations) {
    ed::FiniteSpectrumOptions options;
    options.sector = std::move(sector);
    options.count = count;
    options.which = which;
    options.tolerance = tolerance;
    options.max_iterations = max_iterations;
    return options;
}

std::vector<double> all_levels(ed::ModelSpectrum const& spectrum) {
    std::vector<double> levels;
    for (auto const& block : spectrum.blocks) levels.insert(levels.end(), block.eigenvalues.begin(), block.eigenvalues.end());
    std::sort(levels.begin(), levels.end());
    return levels;
}

#define PYALPS_MODEL_SPECTRUM_ARGS                                                                               \
    nb::arg("sector") = std::map<std::string, double>{}, nb::arg("momentum") = std::vector<double>{},            \
    nb::arg("all_momenta") = false, PYALPS_SPECTRUM_ARGS

void bind_projected(nb::module_& m) {
    nb::class_<ed::MomentumSelection>(m, "MomentumSelection", "A row of the momentum grid and its canonical components.")
        .def_ro("index", &ed::MomentumSelection::index)
        .def_ro("momentum", &ed::MomentumSelection::momentum);
    m.def("select_momentum", &ed::select_momentum, nb::arg("periods"), nb::arg("wavevector"),
          "Snap a wavevector (radians per translation) onto the grid 2*pi*j/period.");

    nb::class_<ed::ProjectedSpaces>(m, "ProjectedSpaces",
        "Sorted conserved charge names and one SpectrumLocalSpace per site with a row per charge.")
        .def(nb::init<>())
        .def("__init__", [](ed::ProjectedSpaces* self, std::vector<std::string> charges, std::vector<ed::SpectrumLocalSpace> local_spaces) {
            new (self) ed::ProjectedSpaces{std::move(charges), std::move(local_spaces)};
        }, nb::arg("charges"), nb::arg("local_spaces"))
        .def_rw("charges", &ed::ProjectedSpaces::charges)
        .def_rw("local_spaces", &ed::ProjectedSpaces::local_spaces)
        .def_prop_ro("num_sites", [](ed::ProjectedSpaces const& s) { return s.local_spaces.size(); });
    m.def("projected_spaces", &ed::projected_spaces, nb::arg("physical_basis"),
          "The parity and charge rows of a complete physical basis (parity +1 under the lattice convention).");

    nb::class_<ed::ProjectedProblem>(m, "ProjectedProblem",
        "A prepared Hermitian problem in a projected basis, its reduced matrix assembled once.")
        .def_prop_ro("basis", &ed::ProjectedProblem::basis)
        .def_prop_ro("dimension", &ed::ProjectedProblem::dimension)
        .def_prop_ro("charges", &ed::ProjectedProblem::charges)
        .def_prop_ro("hamiltonian", &ed::ProjectedProblem::hamiltonian, "The reduced Hamiltonian as canonical complex CSR.")
        .def("eigenstates", [](ed::ProjectedProblem const& p, int count, std::string const& method, std::string const& which,
                               std::optional<double> tolerance, std::optional<int> max_iterations) {
            auto const options = eigen_options(count, method, which, tolerance, max_iterations);
            nb::gil_scoped_release release;
            return p.eigenstates(options);
        }, PYALPS_EIGEN_ARGS, "Eigenpairs in the projected basis (vectors are complex for every momentum).")
        .def("__repr__", [](ed::ProjectedProblem const& p) { return "ProjectedProblem(dimension=" + std::to_string(p.dimension()) + ")"; });
    nb::class_<ed::ProjectedEigenResult>(m, "ProjectedEigenResult", "Eigenpairs of a ProjectedProblem.")
        .def_ro("problem", &ed::ProjectedEigenResult::problem)
        .def_ro("spectrum", &ed::ProjectedEigenResult::spectrum, "The SpectrumBlock: eigenvalues, method, momentum, controls.")
        .def_prop_ro("eigenvalues", [](ed::ProjectedEigenResult const& r) {
            return vector_to_numpy(std::vector<double>(r.spectrum.eigenvalues));
        })
        .def_prop_ro("eigenvectors", [](ed::ProjectedEigenResult const& r) {
            return column_major(std::vector<Complex>(r.eigenvectors_column_major), r.problem.dimension(), r.spectrum.eigenvalues.size());
        }, "(dimension, count) complex array in the projected basis; column j belongs to eigenvalue j.")
        .def_prop_ro("method", [](ed::ProjectedEigenResult const& r) { return r.spectrum.method; })
        .def_prop_ro("residual_norms", [](ed::ProjectedEigenResult const& r) { return vector_to_numpy(std::vector<double>(r.residual_norms)); })
        .def_ro("sparse_requested", &ed::ProjectedEigenResult::sparse_requested)
        .def("expectation_values", [](ed::ProjectedEigenResult const& r, mdl::CompiledOperator const& op) {
            std::size_t const dimension = r.problem.dimension();
            std::size_t const count = r.spectrum.eigenvalues.size();
            std::vector<Complex> values;
            {
                nb::gil_scoped_release release;
                for (std::size_t column = 0; column < count; ++column) {
                    std::vector<Complex> const state(r.eigenvectors_column_major.begin() + column * dimension,
                                                     r.eigenvectors_column_major.begin() + (column + 1) * dimension);
                    values.push_back(ed::expect_projected_operator(op, r.problem.basis(), state));
                }
            }
            return vector_to_numpy(std::move(values));
        }, nb::arg("op"), "<psi_i|O|psi_i> for every eigenvector; O must preserve the basis's symmetry sector.")
        .def("__len__", [](ed::ProjectedEigenResult const& r) { return r.spectrum.eigenvalues.size(); })
        .def("__repr__", [](ed::ProjectedEigenResult const& r) {
            return "ProjectedEigenResult(" + std::to_string(r.spectrum.eigenvalues.size()) + " states, " + r.spectrum.method + ")";
        });

    m.def("prepare_projected", [](mdl::CompiledOperator const& op, ed::ProjectedBasis const& basis) {
        nb::gil_scoped_release release;
        return ed::prepare_projected(op, basis);
    }, nb::arg("op"), nb::arg("basis"), "Validate and assemble an operator on an explicitly constructed basis.");
    m.def("prepare_projected", [](mdl::CompiledOperator const& op, ed::ProjectedSpaces const& spaces,
                                  std::vector<ed::SpectrumTranslation> const& translations, std::map<std::string, double> sector,
                                  std::vector<double> momentum) {
        ed::ProjectedOptions options;
        options.sector = std::move(sector);
        options.momentum = std::move(momentum);
        nb::gil_scoped_release release;
        return ed::prepare_projected(op, spaces, translations, options);
    }, nb::arg("op"), nb::arg("spaces"), nb::arg("translations") = std::vector<ed::SpectrumTranslation>{},
       nb::arg("sector") = std::map<std::string, double>{}, nb::arg("momentum") = std::vector<double>{},
       "A compiled operator on the charge sector `sector` at `momentum` (one component per translation).");
    m.def("prepare_projected", [](mdl::TermSum const& terms, std::vector<mdl::SiteBasisPtr> const& sites,
                                  std::vector<std::string> const& conserved, nb::handle convention,
                                  std::vector<ed::SpectrumTranslation> const& translations, std::map<std::string, double> sector,
                                  std::vector<double> momentum) {
        auto const conv = convention_of(convention);
        ed::ProjectedOptions options;
        options.sector = std::move(sector);
        options.momentum = std::move(momentum);
        nb::gil_scoped_release release;
        return ed::prepare_projected(terms, sites, conserved, conv, translations, options);
    }, nb::arg("terms"), nb::arg("sites"), nb::arg("conserved_quantum_numbers") = std::vector<std::string>{},
       nb::arg("fermion_convention") = "strict", nb::arg("translations") = std::vector<ed::SpectrumTranslation>{},
       nb::arg("sector") = std::map<std::string, double>{}, nb::arg("momentum") = std::vector<double>{});
    m.def("prepare_projected", [](mdl::Hamiltonian const& model, std::map<std::string, double> sector, std::vector<double> momentum) {
        ed::ProjectedOptions options;
        options.sector = std::move(sector);
        options.momentum = std::move(momentum);
        nb::gil_scoped_release release;
        return ed::prepare_projected(model, options);
    }, nb::arg("model"), nb::arg("sector") = std::map<std::string, double>{}, nb::arg("momentum") = std::vector<double>{},
       "A model on the charge sector `sector`; `momentum` (one component per periodic lattice direction) selects "
       "the translation-projected basis, empty a charge-only basis.");

    m.def("prepare_finite", [](mdl::CompiledOperator const& op, ed::ProjectedSpaces const& spaces, ed::FiniteGroupIrrep const& irrep,
                               std::map<std::string, double> const& sector) {
        nb::gil_scoped_release release;
        return ed::prepare_finite(op, spaces, irrep, sector);
    }, nb::arg("op"), nb::arg("spaces"), nb::arg("irrep"), nb::arg("sector") = std::map<std::string, double>{});
    m.def("prepare_finite", [](mdl::TermSum const& terms, std::vector<mdl::SiteBasisPtr> const& sites, ed::FiniteGroupIrrep const& irrep,
                               std::vector<std::string> const& conserved, nb::handle convention,
                               std::map<std::string, double> const& sector) {
        auto const conv = convention_of(convention);
        nb::gil_scoped_release release;
        return ed::prepare_finite(terms, sites, conserved, conv, irrep, sector);
    }, nb::arg("terms"), nb::arg("sites"), nb::arg("irrep"), nb::arg("conserved_quantum_numbers") = std::vector<std::string>{},
       nb::arg("fermion_convention") = "strict", nb::arg("sector") = std::map<std::string, double>{});
    m.def("prepare_finite", [](mdl::Hamiltonian const& model, ed::FiniteGroupIrrep const& irrep, std::map<std::string, double> const& sector) {
        nb::gil_scoped_release release;
        return ed::prepare_finite(model, irrep, sector);
    }, nb::arg("model"), nb::arg("irrep"), nb::arg("sector") = std::map<std::string, double>{},
       "A model on the isotypic component of `irrep` within the charge sector.");

    nb::class_<ed::ModelSpectrum>(m, "ModelSpectrum", "Blocks of a model's spectrum resolved by charges and momenta or irreps.")
        .def_ro("charges", &ed::ModelSpectrum::charges, "Sorted charge names: the order of SpectrumBlock.quantum_numbers.")
        .def_ro("momentum", &ed::ModelSpectrum::momentum)
        .def_ro("blocks", &ed::ModelSpectrum::blocks)
        .def_ro("sector", &ed::ModelSpectrum::sector)
        .def_ro("directions", &ed::ModelSpectrum::directions)
        .def_prop_ro("eigenvalues", [](ed::ModelSpectrum const& s) { return vector_to_numpy(all_levels(s)); },
                     "Every block's eigenvalues merged in ascending order.")
        .def("__len__", [](ed::ModelSpectrum const& s) { return s.blocks.size(); })
        .def("__repr__", [](ed::ModelSpectrum const& s) {
            return "ModelSpectrum(" + std::to_string(s.blocks.size()) + " blocks, " + std::to_string(all_levels(s).size()) + " levels)";
        });
    m.def("model_spectrum", [](mdl::Hamiltonian const& model, std::map<std::string, double> sector, std::vector<double> momentum,
                               bool all_momenta, std::optional<int> count, std::string const& which, std::optional<double> tolerance,
                               std::optional<int> max_iterations) {
        auto const options = model_options(std::move(sector), std::move(momentum), all_momenta, count, which, tolerance, max_iterations);
        nb::gil_scoped_release release;
        return ed::model_spectrum(model, options);
    }, nb::arg("model"), PYALPS_MODEL_SPECTRUM_ARGS,
       "The spectrum resolved by every declared charge not fixed in `sector` and, with `momentum` or `all_momenta`, "
       "by lattice momentum; `count` limits the levels per block.");
    m.def("model_spectrum", [](mdl::TermSum const& terms, std::vector<mdl::SiteBasisPtr> const& sites,
                               std::vector<std::string> const& conserved, nb::handle convention,
                               std::vector<ed::SpectrumTranslation> const& translations, std::map<std::string, double> sector,
                               std::vector<double> momentum, bool all_momenta, std::optional<int> count, std::string const& which,
                               std::optional<double> tolerance, std::optional<int> max_iterations) {
        auto const conv = convention_of(convention);
        auto const options = model_options(std::move(sector), std::move(momentum), all_momenta, count, which, tolerance, max_iterations);
        nb::gil_scoped_release release;
        return ed::model_spectrum(terms, sites, conserved, conv, translations, options);
    }, nb::arg("terms"), nb::arg("sites"), nb::arg("conserved_quantum_numbers") = std::vector<std::string>{},
       nb::arg("fermion_convention") = "strict", nb::arg("translations") = std::vector<ed::SpectrumTranslation>{},
       PYALPS_MODEL_SPECTRUM_ARGS);
    m.def("model_spectrum", [](mdl::CompiledOperator const& op, ed::ProjectedSpaces const& spaces,
                               std::vector<ed::SpectrumTranslation> const& translations, std::map<std::string, double> sector,
                               std::vector<double> momentum, bool all_momenta, std::optional<int> count, std::string const& which,
                               std::optional<double> tolerance, std::optional<int> max_iterations) {
        auto const options = model_options(std::move(sector), std::move(momentum), all_momenta, count, which, tolerance, max_iterations);
        nb::gil_scoped_release release;
        return ed::model_spectrum(op, spaces, translations, options);
    }, nb::arg("op"), nb::arg("spaces"), nb::arg("translations") = std::vector<ed::SpectrumTranslation>{}, PYALPS_MODEL_SPECTRUM_ARGS);

    m.def("finite_spectrum", [](mdl::Hamiltonian const& model, std::vector<ed::FiniteGroupIrrep> const& irreps,
                                std::map<std::string, double> sector, std::optional<int> count, std::string const& which,
                                std::optional<double> tolerance, std::optional<int> max_iterations) {
        auto const options = finite_options(std::move(sector), count, which, tolerance, max_iterations);
        nb::gil_scoped_release release;
        return ed::finite_spectrum(model, irreps, options);
    }, nb::arg("model"), nb::arg("irreps"), nb::arg("sector") = std::map<std::string, double>{}, PYALPS_SPECTRUM_ARGS,
       "The spectrum resolved by every irrep of a complete family and every unfixed declared charge.");
    m.def("finite_spectrum", [](mdl::TermSum const& terms, std::vector<mdl::SiteBasisPtr> const& sites,
                                std::vector<ed::FiniteGroupIrrep> const& irreps, std::vector<std::string> const& conserved,
                                nb::handle convention, std::map<std::string, double> sector,
                                std::optional<int> count, std::string const& which, std::optional<double> tolerance,
                                std::optional<int> max_iterations) {
        auto const conv = convention_of(convention);
        auto const options = finite_options(std::move(sector), count, which, tolerance, max_iterations);
        nb::gil_scoped_release release;
        return ed::finite_spectrum(terms, sites, conserved, conv, irreps, options);
    }, nb::arg("terms"), nb::arg("sites"), nb::arg("irreps"), nb::arg("conserved_quantum_numbers") = std::vector<std::string>{},
       nb::arg("fermion_convention") = "strict", nb::arg("sector") = std::map<std::string, double>{},
       PYALPS_SPECTRUM_ARGS);
    m.def("finite_spectrum", [](mdl::CompiledOperator const& op, ed::ProjectedSpaces const& spaces,
                                std::vector<ed::FiniteGroupIrrep> const& irreps, std::map<std::string, double> sector,
                                std::optional<int> count, std::string const& which, std::optional<double> tolerance,
                                std::optional<int> max_iterations) {
        auto const options = finite_options(std::move(sector), count, which, tolerance, max_iterations);
        nb::gil_scoped_release release;
        return ed::finite_spectrum(op, spaces, irreps, options);
    }, nb::arg("op"), nb::arg("spaces"), nb::arg("irreps"), nb::arg("sector") = std::map<std::string, double>{}, PYALPS_SPECTRUM_ARGS);
}

// -- thermodynamics ---------------------------------------------------------------------

void bind_thermodynamics(nb::module_& m) {
    m.def("thermal_weights", [](nb::handle energies, double temperature) {
        auto const levels = vector_from<double>(energies, "energies");
        return vector_to_numpy(ed::thermal_weights(levels, temperature));
    }, nb::arg("energies"), nb::arg("temperature"),
       "Normalized Boltzmann weights of a finite spectrum at a finite nonzero temperature (k_B = 1).");
    m.def("thermal_surface", [](nb::handle energies, nb::handle temperatures, nb::handle fields, std::optional<nb::handle> charges,
                                double field0) {
        auto const levels = vector_from<double>(energies, "energies");
        auto const grid_t = vector_from<double>(temperatures, "temperatures");
        auto const grid_h = vector_from<double>(fields, "fields");
        std::vector<double> const native_charges = charges && !charges->is_none() ? vector_from<double>(*charges, "charges")
                                                                                 : std::vector<double>{};
        std::vector<ed::ThermalAverages> surface;
        {
            nb::gil_scoped_release release;
            surface = ed::thermal_surface(levels, native_charges, grid_t, grid_h, field0);
        }
        nb::dict result;
        auto add = [&](char const* name, double ed::ThermalAverages::*member) {
            std::vector<double> values;
            values.reserve(surface.size());
            for (auto const& point : surface) values.push_back(point.*member);
            result[name] = matrix_to_numpy(std::move(values), grid_t.size(), grid_h.size());
        };
        using A = ed::ThermalAverages;
        add("energy", &A::energy);
        add("free_energy", &A::free_energy);
        add("entropy", &A::entropy);
        add("specific_heat", &A::specific_heat);
        add("charge", &A::charge);
        add("susceptibility", &A::susceptibility);
        add("partition", &A::partition);
        add("energy_reference", &A::energy_reference);
        add("log_partition", &A::log_partition);
        return result;
    }, nb::arg("energies"), nb::arg("temperatures"), nb::arg("fields") = std::vector<double>{0.}, nb::arg("charges") = nb::none(),
       nb::arg("field0") = 0.,
       "Canonical averages on a temperature x field grid as (len(temperatures), len(fields)) arrays: energy, free_energy, "
       "entropy, specific_heat, charge, susceptibility, partition, energy_reference, log_partition. `charges` gives one "
       "conserved charge per level; every level shifts by -charge * (field - field0).");
}

} // namespace

NB_MODULE(ed_c, m) {
    m.doc() = "ALPS::ed: exact diagonalization of compiled operators";
    bind_sparse_assembly(m);
    bind_eigensolvers(m);
    bind_product(m);
    bind_symmetry(m);
    bind_projected(m);
    bind_thermodynamics(m);
}
