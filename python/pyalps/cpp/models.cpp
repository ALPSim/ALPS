// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// pyalps.models: the ALPS::models kit (lattices, local bases, the operator
// algebra, Hamiltonians, compiled operators, physical bases and the Pauli
// export), the lattice-to-XML bridge of the legacy lattice module and the
// XML catalog importer. The classes are bound directly; dense matrices and
// tables cross as NumPy arrays, enumerations also accept their names.
//
// The surface follows the C++ headers of src/alps/models and
// src/alps/model/catalog.hpp; the earlier Python package of the
// skilledwolf/alpsx fork (python/alps/cpp/models_kit_bindings.cpp,
// lattice_bindings.cpp at d37e51227) exchanged dictionaries and tuples with a
// Python-side operator algebra instead and was consulted for the GIL and
// ownership conventions only.
#include "kit_support.hpp"

#include <alps/lattice/kit.hpp>
#include <alps/model/catalog.hpp>
#include <alps/models.hpp>
#include <alps/parameter.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;
namespace lat = alps::lattice;
namespace mdl = alps::models;
using namespace pyalps::kit;

namespace {

// -- legacy alps::Parameters <-> dict ----------------------------------------------

alps::Parameters parameters_from_dict(nb::dict const& values) {
    alps::Parameters parameters;
    for (auto item : values) {
        std::string const key = nb::cast<std::string>(nb::str(item.first));
        nb::handle value = item.second;
        if (nb::isinstance<nb::str>(value)) {
            parameters[key] = nb::cast<std::string>(value);
        } else if (PyBool_Check(value.ptr()) || alps::python::type_name(value).compare(0, 10, "numpy.bool") == 0) {
            parameters[key] = PyObject_IsTrue(value.ptr()) == 1 ? "true" : "false";
        } else if (PyIndex_Check(value.ptr())) {
            nb::object integer = nb::steal(PyNumber_Index(value.ptr()));
            if (!integer.is_valid()) throw nb::python_error();
            parameters[key] = nb::cast<std::string>(nb::str(integer));
        } else if (PyFloat_Check(value.ptr()) || PyObject_HasAttrString(value.ptr(), "__float__")) {
            double const number = PyFloat_AsDouble(value.ptr());
            if (PyErr_Occurred()) throw nb::python_error();
            char buffer[32];
            std::snprintf(buffer, sizeof(buffer), "%.17g", number);
            parameters[key] = buffer;
        } else {
            throw nb::type_error(("parameter '" + key + "' must be a string, bool, int or float").c_str());
        }
    }
    return parameters;
}

nb::dict dict_from_parameters(alps::Parameters const& parameters) {
    nb::dict result;
    for (auto const& item : parameters)
        result[nb::str(item.key().c_str())] = static_cast<std::string>(item.value());
    return result;
}

// -- repr helpers ----------------------------------------------------------------

template <typename T>
std::string join(std::vector<T> const& values) {
    std::string out = "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out += ", ";
        if constexpr (std::is_same_v<T, std::string>) out += "'" + values[i] + "'";
        else out += std::to_string(values[i]);
    }
    return out + "]";
}

std::string join_boundary(std::vector<lat::Boundary> const& values) {
    std::string out = "[";
    for (std::size_t i = 0; i < values.size(); ++i) out += (i ? ", '" : "'") + lat::to_string(values[i]) + "'";
    return out + "]";
}

// -- lattices ------------------------------------------------------------------

void bind_lattice(nb::module_& m) {
    nb::enum_<lat::Boundary>(m, "Boundary", "Per-dimension boundary condition; functions also accept 'periodic'/'open'.")
        .value("periodic", lat::Boundary::periodic)
        .value("open", lat::Boundary::open);
    m.def("boundary_from_string", &lat::boundary_from_string, nb::arg("name"));

    nb::class_<lat::Site>(m, "Site", "A site of the tiled cluster.")
        .def_ro("index", &lat::Site::index)
        .def_ro("site_type", &lat::Site::site_type)
        .def_ro("sublattice", &lat::Site::sublattice, "Vertex index within the unit cell.")
        .def_ro("coordinate", &lat::Site::coordinate)
        .def("__repr__", [](lat::Site const& s) {
            return "Site(index=" + std::to_string(s.index) + ", site_type=" + std::to_string(s.site_type)
                 + ", coordinate=" + join(s.coordinate) + ")";
        });
    nb::class_<lat::Bond>(m, "Bond", "A bond of the tiled cluster; `vector` is target minus source coordinates.")
        .def_ro("source", &lat::Bond::source)
        .def_ro("target", &lat::Bond::target)
        .def_ro("bond_type", &lat::Bond::bond_type)
        .def_ro("vector", &lat::Bond::vector)
        .def_ro("wraps", &lat::Bond::wraps, "True when the periodic modulo folded an endpoint back into the cluster.")
        .def("__repr__", [](lat::Bond const& b) {
            return "Bond(source=" + std::to_string(b.source) + ", target=" + std::to_string(b.target)
                 + ", bond_type=" + std::to_string(b.bond_type) + (b.wraps ? ", wraps=True)" : ")");
        });
    nb::class_<lat::UnitCell::Edge>(m, "Edge", "A unit-cell edge: vertices of (offset) cells with a bond type.")
        .def_ro("source", &lat::UnitCell::Edge::source)
        .def_ro("target", &lat::UnitCell::Edge::target)
        .def_ro("edge_type", &lat::UnitCell::Edge::edge_type)
        .def_ro("source_offset", &lat::UnitCell::Edge::source_offset)
        .def_ro("target_offset", &lat::UnitCell::Edge::target_offset);

    nb::class_<lat::UnitCell>(m, "UnitCell",
        "Repeating motif: typed vertices with fractional coordinates and offset edges with bond types.")
        .def(nb::init<int>(), nb::arg("dimension"))
        .def("add_vertex", [](lat::UnitCell& cell, std::vector<double> coordinate, int vertex_type) {
            return cell.add_vertex(std::move(coordinate), vertex_type);
        }, nb::arg("coordinate"), nb::arg("vertex_type") = 0, "Add a vertex and return its index.")
        .def("add_edge", [](lat::UnitCell& cell, int source, int target, int edge_type,
                            std::vector<int> source_offset, std::vector<int> target_offset) {
            cell.add_edge(source, target, edge_type, std::move(source_offset), std::move(target_offset));
        }, nb::arg("source"), nb::arg("target"), nb::arg("edge_type") = 0,
           nb::arg("source_offset") = std::vector<int>{}, nb::arg("target_offset") = std::vector<int>{},
           "Connect two vertices of (offset) cells; offsets default to the home cell.")
        .def_prop_ro("dimension", &lat::UnitCell::dimension)
        .def_prop_ro("num_vertices", &lat::UnitCell::num_vertices)
        .def_prop_ro("vertex_coords", &lat::UnitCell::vertex_coords)
        .def_prop_ro("vertex_types", &lat::UnitCell::vertex_types)
        .def_prop_ro("edges", &lat::UnitCell::edges)
        .def("tile", [](lat::UnitCell const& cell, std::vector<int> const& extent,
                        std::vector<std::vector<double>> basis_vectors, nb::handle bc) {
            return cell.tile(extent, std::move(basis_vectors), boundaries_of(bc, cell.dimension()));
        }, nb::arg("extent"), nb::arg("basis_vectors") = std::vector<std::vector<double>>{},
           nb::arg("bc") = "periodic",
           "Tile `extent` cells per dimension; `basis_vectors` rows are the primitive vectors "
           "(identity when empty), `bc` one boundary for all dimensions or one per dimension.");

    nb::class_<lat::Lattice>(m, "Lattice",
        "A finite lattice: enumerated sites and bonds plus unit cell, basis, extent and boundary. "
        "Cells enumerate with the last dimension fastest; site = cell * num_vertices + vertex.")
        .def_static("build", [](lat::UnitCell unit_cell, std::vector<std::vector<double>> basis_vectors,
                                std::vector<int> extent, nb::handle boundary, std::string alps_name) {
            auto bc = boundaries_of(boundary, unit_cell.dimension());
            return lat::Lattice::build(std::move(unit_cell), std::move(basis_vectors), std::move(extent),
                                       std::move(bc), std::move(alps_name));
        }, nb::arg("unit_cell"), nb::arg("basis_vectors"), nb::arg("extent"), nb::arg("boundary"),
           nb::arg("alps_name") = "", "The tiler behind the named factories.")
        .def_prop_ro("sites", &lat::Lattice::sites)
        .def_prop_ro("bonds", &lat::Lattice::bonds)
        .def_prop_ro("basis_vectors", &lat::Lattice::basis_vectors)
        .def_prop_ro("extent", &lat::Lattice::extent)
        .def_prop_ro("boundary", &lat::Lattice::boundary)
        .def_prop_ro("unit_cell", &lat::Lattice::unit_cell)
        .def_prop_ro("num_sites", &lat::Lattice::num_sites)
        .def_prop_ro("num_bonds", &lat::Lattice::num_bonds)
        .def_prop_ro("dimension", &lat::Lattice::dimension)
        .def_prop_ro("alps_name", &lat::Lattice::alps_name, "The catalog LATTICEGRAPH base name, '' when custom.")
        .def_prop_ro("is_custom", &lat::Lattice::is_custom)
        .def_prop_ro("coordinates", [](lat::Lattice const& lattice) {
            std::size_t const width = lattice.sites().empty() ? 0 : lattice.sites().front().coordinate.size();
            std::vector<double> values;
            values.reserve(lattice.sites().size() * width);
            for (auto const& site : lattice.sites())
                values.insert(values.end(), site.coordinate.begin(), site.coordinate.end());
            return matrix_to_numpy(std::move(values), lattice.sites().size(), width);
        }, "Real-space site coordinates, one row per site.")
        .def("neighbors", [](lat::Lattice const& lattice, int site) { return lattice.neighbors(site); }, nb::arg("site"))
        .def("sites_by_type", &lat::Lattice::sites_by_type, nb::arg("site_type"))
        .def("bonds_by_type", &lat::Lattice::bonds_by_type, nb::arg("bond_type"))
        .def("periodic_directions", &lat::Lattice::periodic_directions)
        .def("translation_generators", &lat::Lattice::translation_generators,
             "One source-site -> translated-site permutation per periodic direction.")
        .def("momenta", [](lat::Lattice const& lattice) {
            auto rows = lattice.momenta();
            std::size_t const width = rows.empty() ? 0 : rows.front().size();
            std::vector<double> values;
            for (auto const& row : rows) values.insert(values.end(), row.begin(), row.end());
            return matrix_to_numpy(std::move(values), rows.size(), width);
        }, "Allowed total momenta: one row per momentum, one column per periodic direction.")
        .def_prop_ro("num_momenta", &lat::Lattice::num_momenta)
        .def("edges", &lat::Lattice::edges, "The graph half of a CompiledModel: one ModelEdge per bond.")
        .def("site_types", &lat::Lattice::site_types)
        .def("__repr__", [](lat::Lattice const& lattice) {
            return "Lattice(" + (lattice.is_custom() ? std::string("custom") : "'" + lattice.alps_name() + "'")
                 + ", extent=" + join(lattice.extent()) + ", boundary=" + join_boundary(lattice.boundary())
                 + ", sites=" + std::to_string(lattice.num_sites()) + ", bonds=" + std::to_string(lattice.num_bonds()) + ")";
        });

    m.def("chain", [](int L, double a, nb::handle bc) { return lat::chain(L, a, boundary_of(bc)); },
          nb::arg("L"), nb::arg("a") = 1.0, nb::arg("bc") = "periodic");
    m.def("nnn_chain", [](int L, double a, nb::handle bc) { return lat::nnn_chain(L, a, boundary_of(bc)); },
          nb::arg("L"), nb::arg("a") = 1.0, nb::arg("bc") = "periodic",
          "Bond type 0 = nearest, 1 = next-nearest neighbours.");
    m.def("square", [](int L, int W, double a, nb::handle bc) { return lat::square(L, W, a, boundary_of(bc)); },
          nb::arg("L"), nb::arg("W") = 0, nb::arg("a") = 1.0, nb::arg("bc") = "periodic");
    m.def("frustrated_square", [](int L, int W, double a, nb::handle bc) {
        return lat::frustrated_square(L, W, a, boundary_of(bc));
    }, nb::arg("L"), nb::arg("W") = 0, nb::arg("a") = 1.0, nb::arg("bc") = "periodic",
       "Bond type 0 = nearest, 1 = diagonal next-nearest neighbours.");
    m.def("triangular", [](int L, int W, double a, nb::handle bc) { return lat::triangular(L, W, a, boundary_of(bc)); },
          nb::arg("L"), nb::arg("W") = 0, nb::arg("a") = 1.0, nb::arg("bc") = "periodic");
    m.def("honeycomb", [](int L, int W, double a, nb::handle bc) { return lat::honeycomb(L, W, a, boundary_of(bc)); },
          nb::arg("L"), nb::arg("W") = 0, nb::arg("a") = 1.0, nb::arg("bc") = "periodic",
          "Two sublattices; the three inequivalent links carry bond types 0/1/2.");
    m.def("kagome", [](int L, int W, double a, nb::handle bc) { return lat::kagome(L, W, a, boundary_of(bc)); },
          nb::arg("L"), nb::arg("W") = 0, nb::arg("a") = 1.0, nb::arg("bc") = "periodic");
    m.def("simple_cubic", [](int L, int W, int H, double a, nb::handle bc) {
        return lat::simple_cubic(L, W, H, a, boundary_of(bc));
    }, nb::arg("L"), nb::arg("W") = 0, nb::arg("H") = 0, nb::arg("a") = 1.0, nb::arg("bc") = "periodic");
    m.def("ladder", [](int L, double a, nb::handle bc) { return lat::ladder(L, a, boundary_of(bc)); },
          nb::arg("L"), nb::arg("a") = 1.0, nb::arg("bc") = "periodic",
          "Two-leg ladder: legs bond type 0, rungs bond type 1; `bc` applies to the legs.");
    m.def("rectangular", [](int L, int W, double a, double b, nb::handle bc) {
        return lat::rectangular(L, W, a, b, boundary_of(bc));
    }, nb::arg("L"), nb::arg("W") = 0, nb::arg("a") = 1.0, nb::arg("b") = 0.0, nb::arg("bc") = "periodic");
    m.def("face_centered_cubic", [](int L, int W, int H, double a, nb::handle bc) {
        return lat::face_centered_cubic(L, W, H, a, boundary_of(bc));
    }, nb::arg("L"), nb::arg("W") = 0, nb::arg("H") = 0, nb::arg("a") = 1.0, nb::arg("bc") = "periodic");
    m.def("body_centered_cubic", [](int L, int W, int H, double a, nb::handle bc) {
        return lat::body_centered_cubic(L, W, H, a, boundary_of(bc));
    }, nb::arg("L"), nb::arg("W") = 0, nb::arg("H") = 0, nb::arg("a") = 1.0, nb::arg("bc") = "periodic");
    m.def("tetragonal", [](int L, int W, int H, double a, double c, nb::handle bc) {
        return lat::tetragonal(L, W, H, a, c, boundary_of(bc));
    }, nb::arg("L"), nb::arg("W") = 0, nb::arg("H") = 0, nb::arg("a") = 1.0, nb::arg("c") = 1.0,
       nb::arg("bc") = "periodic");
    m.def("orthorhombic", [](int L, int W, int H, double a, double b, double c, nb::handle bc) {
        return lat::orthorhombic(L, W, H, a, b, c, boundary_of(bc));
    }, nb::arg("L"), nb::arg("W") = 0, nb::arg("H") = 0, nb::arg("a") = 1.0, nb::arg("b") = 1.0,
       nb::arg("c") = 1.0, nb::arg("bc") = "periodic");
    m.def("hexagonal", [](int L, int W, int H, double a, double c, nb::handle bc) {
        return lat::hexagonal(L, W, H, a, c, boundary_of(bc));
    }, nb::arg("L"), nb::arg("W") = 0, nb::arg("H") = 0, nb::arg("a") = 1.0, nb::arg("c") = 1.0,
       nb::arg("bc") = "periodic", "The three-dimensional hexagonal lattice.");

    // -- the XML bridge ------------------------------------------------------
    m.def("to_latticegraph_xml", &lat::to_latticegraph_xml, nb::arg("lattice"), nb::arg("name") = "alps_custom_lattice",
          "The <LATTICES> document (basis, unit cell with vertex coordinates, LATTICEGRAPH with extent and boundary).");
    m.def("to_graph_xml", &lat::to_graph_xml, nb::arg("lattice"),
          "The flat <GRAPH> of the tiled cluster with real-space vertex coordinates.");
    m.def("custom_lattice_name", &lat::custom_lattice_name, nb::arg("lattice"),
          "Content-addressed name 'alps_custom_lattice_<16 hex digits>' of the lattice's document.");
    m.def("catalog_graph_name", &lat::catalog_graph_name, nb::arg("lattice"),
          "The catalog LATTICEGRAPH the lattice's boundary selects, '' when there is none.");

    nb::class_<lat::LatticeResolution>(m, "LatticeResolution",
        "How a kit lattice reaches the legacy lattice library: a catalog graph with its L/W/H and bc "
        "parameters, or a custom <LATTICES> document.")
        .def_prop_ro("parameters", [](lat::LatticeResolution const& r) { return dict_from_parameters(r.parameters); },
                     "LATTICE and, on the catalog route, the extent and boundary parameters.")
        .def_ro("graph", &lat::LatticeResolution::graph, "The LATTICEGRAPH name selected.")
        .def_ro("library_xml", &lat::LatticeResolution::library_xml, "The document on the custom route, '' otherwise.")
        .def_prop_ro("from_catalog", &lat::LatticeResolution::from_catalog)
        .def("__repr__", [](lat::LatticeResolution const& r) {
            return "LatticeResolution(graph='" + r.graph + (r.from_catalog() ? "', catalog)" : "', custom)");
        });
    m.def("resolve_lattice", [](lat::Lattice const& lattice, nb::dict const& catalog) {
        return lat::resolve_lattice(lattice, parameters_from_dict(catalog));
    }, nb::arg("lattice"), nb::arg("catalog") = nb::dict(),
       "The catalog route when the library selected by `catalog` (LATTICE_LIBRARY, else lattices.xml) "
       "has the graph and expands to the identical structure; the custom document otherwise.");
    m.def("custom_lattice", &lat::custom_lattice, nb::arg("lattice"), "Always the custom document.");
    m.def("write_lattice_library", [](lat::LatticeResolution const& resolution, std::string const& file) {
        return dict_from_parameters(lat::write_lattice_library(resolution, std::filesystem::path(file)));
    }, nb::arg("resolution"), nb::arg("file"),
       "Write the custom document to `file` and return the parameters with LATTICE_LIBRARY set; "
       "the catalog route writes nothing.");
}

// -- local bases ---------------------------------------------------------------

mdl::LocalOp local_op_from(std::string const& name, nb::handle value) {
    if (nb::isinstance<mdl::LocalOp>(value)) {
        mdl::LocalOp op = nb::cast<mdl::LocalOp>(value);
        if (op.name.empty()) op.name = name;
        return op;
    }
    mdl::LocalOp op;
    op.name = name;
    std::size_t dim = 0;
    op.matrix = square_matrix_from(value, dim, ("operator '" + name + "'").c_str());
    op.dim = static_cast<int>(dim);
    return op;
}

void bind_sites(nb::module_& m) {
    nb::enum_<mdl::Statistics>(m, "Statistics")
        .value("bosonic", mdl::Statistics::bosonic)
        .value("fermionic", mdl::Statistics::fermionic);
    nb::enum_<mdl::FermionConvention>(m, "FermionConvention",
        "'strict': fermionic operators anticommute across sites (Jordan-Wigner strings); "
        "'lattice': every operator is treated as bosonic.")
        .value("strict", mdl::FermionConvention::strict)
        .value("lattice", mdl::FermionConvention::lattice);
    nb::enum_<mdl::Parity>(m, "Parity")
        .value("even", mdl::Parity::even)
        .value("odd", mdl::Parity::odd);

    nb::class_<mdl::LocalOp>(m, "LocalOp", "A named local operator: dense matrix, fermionic parity and charge deltas.")
        .def("__init__", [](mdl::LocalOp* self, std::string name, nb::handle matrix, nb::handle parity,
                            std::map<std::string, double> qn_deltas) {
            std::size_t dim = 0;
            auto values = square_matrix_from(matrix, dim, "LocalOp matrix");
            new (self) mdl::LocalOp{std::move(name), static_cast<int>(dim), std::move(values), parity_of(parity),
                                    std::move(qn_deltas)};
        }, nb::arg("name"), nb::arg("matrix"), nb::arg("parity") = "even",
           nb::arg("qn_deltas") = std::map<std::string, double>{})
        .def_ro("name", &mdl::LocalOp::name)
        .def_ro("dim", &mdl::LocalOp::dim)
        .def_prop_ro("matrix", [](mdl::LocalOp const& op) { return square_matrix(op.matrix, static_cast<std::size_t>(op.dim)); },
                     "Row-major (row = output state, column = input state).")
        .def_ro("parity", &mdl::LocalOp::parity)
        .def_ro("qn_deltas", &mdl::LocalOp::qn_deltas)
        .def("__call__", [](mdl::LocalOp const& op, int row, int column) { return op(row, column); },
             nb::arg("row"), nb::arg("column"))
        .def("is_identity", &mdl::LocalOp::is_identity, nb::arg("tol") = 1e-12)
        .def("dag", &mdl::LocalOp::dag, "Hermitian conjugate; deltas negate, the name gains or loses '_dag'.")
        .def("split_by_qn_delta", &mdl::LocalOp::split_by_qn_delta, nb::arg("quantum_numbers"),
             "Components of definite charge change, in sorted delta order; they sum to the operator.")
        .def("__repr__", [](mdl::LocalOp const& op) {
            return "LocalOp('" + op.name + "', dim=" + std::to_string(op.dim)
                 + (op.parity == mdl::Parity::odd ? ", odd)" : ")");
        });
    m.def("compute_qn_deltas", [](nb::handle matrix, std::map<std::string, std::vector<double>> const& quantum_numbers,
                                  double tol) {
        std::size_t dim = 0;
        auto values = square_matrix_from(matrix, dim, "matrix");
        return mdl::compute_qn_deltas(values, static_cast<int>(dim), quantum_numbers, tol);
    }, nb::arg("matrix"), nb::arg("quantum_numbers"), nb::arg("tol") = 1e-9,
       "The charge change of a matrix for every quantum number it respects.");

    nb::class_<mdl::SiteBasis>(m, "SiteBasis",
        "A local Hilbert space: dimension, statistics, named operators, per-state quantum numbers and the "
        "diagonal parity operator. Built-in sites come from spin(), boson(), hardcore_boson(), "
        "spinless_fermion() and spinful_fermion().")
        .def("__init__", [](mdl::SiteBasis* self, int dim, nb::dict const& operators,
                            std::map<std::string, std::vector<double>> quantum_numbers, nb::handle statistics,
                            std::optional<std::vector<double>> parity_operator, std::vector<std::string> labels,
                            std::vector<std::string> odd_operators) {
            mdl::SiteBasis site;
            site.dim = dim;
            site.statistics = statistics_of(statistics);
            site.quantum_numbers = std::move(quantum_numbers);
            if (parity_operator) site.parity_operator = std::move(*parity_operator);
            site.basis_labels = std::move(labels);
            std::set<std::string> const odd(odd_operators.begin(), odd_operators.end());
            for (auto item : operators) {
                std::string const name = nb::cast<std::string>(nb::str(item.first));
                mdl::LocalOp op = local_op_from(name, item.second);
                if (op.dim != dim)
                    throw nb::value_error(("operator '" + name + "' does not have the site dimension").c_str());
                if (odd.count(name)) op.parity = mdl::Parity::odd;
                site.operators[name] = std::move(op);
            }
            for (auto const& name : odd)
                if (!site.operators.count(name))
                    throw nb::value_error(("odd operator '" + name + "' is not defined").c_str());
            new (self) mdl::SiteBasis(*mdl::finalize_site(std::move(site)));
        }, nb::arg("dim"), nb::arg("operators"), nb::arg("quantum_numbers") = std::map<std::string, std::vector<double>>{},
           nb::arg("statistics") = "bosonic", nb::arg("parity_operator") = nb::none(),
           nb::arg("labels") = std::vector<std::string>{}, nb::arg("odd_operators") = std::vector<std::string>{},
           "A custom site: `operators` maps names to square matrices (or LocalOp values); `odd_operators` "
           "names the fermionic ones. The identity 'I', a +1 parity operator and the charge deltas are "
           "filled in as finalize_site does.")
        .def_ro("dim", &mdl::SiteBasis::dim)
        .def_ro("statistics", &mdl::SiteBasis::statistics)
        .def_ro("operators", &mdl::SiteBasis::operators)
        .def_prop_ro("operator_names", [](mdl::SiteBasis const& site) {
            std::vector<std::string> names;
            for (auto const& item : site.operators) names.push_back(item.first);
            return names;
        })
        .def_ro("quantum_numbers", &mdl::SiteBasis::quantum_numbers)
        .def_ro("parity_operator", &mdl::SiteBasis::parity_operator)
        .def_ro("basis_labels", &mdl::SiteBasis::basis_labels)
        .def("op", &mdl::SiteBasis::op, nb::arg("name"),
             "The operator `name`, or the built-in one the XML catalog spells `name` (Splus, cdag_up, n_down).")
        .def("has_op", &mdl::SiteBasis::has_op, nb::arg("name"))
        .def("adjoint", &mdl::SiteBasis::adjoint, nb::arg("op"),
             "op.dag() under this site's name for it (Sp -> Sm, c_up -> c_dag_up, Sz -> Sz); "
             "TermSum.dag() names its factors this way.")
        .def("__repr__", [](mdl::SiteBasis const& site) {
            std::string names;
            for (auto const& item : site.quantum_numbers) names += (names.empty() ? "" : ", ") + item.first;
            return "SiteBasis(dim=" + std::to_string(site.dim)
                 + (site.statistics == mdl::Statistics::fermionic ? ", fermionic" : ", bosonic")
                 + ", quantum_numbers=[" + names + "])";
        });
    m.def("finalize_site", [](mdl::SiteBasis site) { return mdl::finalize_site(std::move(site)); }, nb::arg("site"));
    {
        nb::dict names;
        for (auto const& [kit, catalog] : mdl::catalog_operator_names)
            names[nb::str(kit.data(), kit.size())] = nb::str(catalog.data(), catalog.size());
        m.attr("catalog_operator_names") = names;  // kit name -> XML catalog name, accepted either way
    }
    m.def("spin", &mdl::spin, nb::arg("S") = 0.5, "Spin S: Sz, Sp, Sm, Sx, Sy, S2; charge Sz; m from +S down.");
    m.def("boson", &mdl::boson, nb::arg("n_max") = 3, "Truncated boson: a, a_dag, n (b aliases); charge N.");
    m.def("hardcore_boson", &mdl::hardcore_boson, "Hard-core boson: a, a_dag, n; charge N.");
    m.def("spinless_fermion", &mdl::spinless_fermion, "Spinless fermion: c, c_dag, n, F; charge N.");
    m.def("spinful_fermion", &mdl::spinful_fermion,
          "Spinful fermion (empty, up, down, double): c_up, c_dn, n_up, n_dn, Sz, ...; charges N and Sz.");

    nb::class_<mdl::LocalBasis>(m, "LocalBasis", "The ordered states of one site: labels, statistics, charges, parities.")
        .def(nb::init<>())
        .def("__init__", [](mdl::LocalBasis* self, std::vector<std::string> labels, nb::handle statistics,
                            std::map<std::string, std::vector<double>> quantum_numbers, std::vector<double> parity) {
            new (self) mdl::LocalBasis{std::move(labels), statistics_of(statistics), std::move(quantum_numbers),
                                       std::move(parity)};
        }, nb::arg("labels"), nb::arg("statistics") = "bosonic",
           nb::arg("quantum_numbers") = std::map<std::string, std::vector<double>>{},
           nb::arg("parity") = std::vector<double>{})
        .def_rw("labels", &mdl::LocalBasis::labels)
        .def_rw("statistics", &mdl::LocalBasis::statistics)
        .def_rw("quantum_numbers", &mdl::LocalBasis::quantum_numbers)
        .def_rw("parity", &mdl::LocalBasis::parity)
        .def_prop_ro("dim", [](mdl::LocalBasis const& b) { return b.labels.size(); })
        .def("__eq__", [](mdl::LocalBasis const& a, mdl::LocalBasis const& b) { return a == b; }, nb::is_operator())
        .def("__ne__", [](mdl::LocalBasis const& a, mdl::LocalBasis const& b) { return a != b; }, nb::is_operator());
    nb::class_<mdl::PhysicalBasis>(m, "PhysicalBasis",
        "The ordered local spaces of a model, its conserved charges (the order of solver charge vectors) "
        "and fermion convention; independent of any Hamiltonian.")
        .def(nb::init<>())
        .def("__init__", [](mdl::PhysicalBasis* self, std::vector<mdl::LocalBasis> sites,
                            std::vector<std::string> conserved, nb::handle convention) {
            new (self) mdl::PhysicalBasis{std::move(sites), std::move(conserved), convention_of(convention)};
        }, nb::arg("sites"), nb::arg("conserved_quantum_numbers") = std::vector<std::string>{},
           nb::arg("fermion_convention") = "strict")
        .def_rw("sites", &mdl::PhysicalBasis::sites)
        .def_rw("conserved_quantum_numbers", &mdl::PhysicalBasis::conserved_quantum_numbers)
        .def_rw("fermion_convention", &mdl::PhysicalBasis::fermion_convention)
        .def_prop_ro("num_sites", [](mdl::PhysicalBasis const& b) { return b.sites.size(); })
        .def_prop_ro("site_dims", [](mdl::PhysicalBasis const& b) {
            std::vector<int> dims;
            for (auto const& site : b.sites) dims.push_back(static_cast<int>(site.labels.size()));
            return dims;
        })
        .def("__eq__", [](mdl::PhysicalBasis const& a, mdl::PhysicalBasis const& b) { return a == b; }, nb::is_operator())
        .def("__ne__", [](mdl::PhysicalBasis const& a, mdl::PhysicalBasis const& b) { return a != b; }, nb::is_operator())
        .def("__repr__", [](mdl::PhysicalBasis const& b) {
            return "PhysicalBasis(" + std::to_string(b.sites.size()) + " sites, conserved="
                 + join(b.conserved_quantum_numbers)
                 + (b.fermion_convention == mdl::FermionConvention::strict ? ", strict)" : ", lattice)");
        });
    m.def("validate_basis", &mdl::validate_basis, nb::arg("basis"));
    m.def("require_same_basis", &mdl::require_same_basis, nb::arg("actual"), nb::arg("expected"));
    m.def("physical_basis", [](std::vector<mdl::SiteBasisPtr> const& sites, std::vector<std::string> const& conserved,
                               nb::handle convention) {
        return mdl::physical_basis(sites, conserved, convention_of(convention));
    }, nb::arg("sites"), nb::arg("conserved_quantum_numbers") = std::vector<std::string>{},
       nb::arg("fermion_convention") = "strict",
       "The complete physical basis of ordered sites; a site lacking a conserved charge is neutral.");
    m.def("physical_basis", [](mdl::Hamiltonian const& model) { return mdl::physical_basis(model); }, nb::arg("model"),
          "The model's local spaces, declared charges and convention.");
    m.def("require_compatible_sites", [](mdl::PhysicalBasis const& reference,
                                         std::vector<std::pair<std::size_t, mdl::SiteBasisPtr>> const& sites) {
        mdl::require_compatible_sites(reference, sites);
    }, nb::arg("reference"), nb::arg("sites"), "Check (site index, SiteBasis) pairs against a complete basis.");
}

// -- the operator algebra ------------------------------------------------------

std::string term_repr(mdl::Term const& term) {
    std::string out = "Term(" + format_complex(term.effective_coeff());
    for (auto const& factor : term.factors) out += " * " + factor.op.name + "@" + std::to_string(factor.site);
    return out + ")";
}

void bind_terms(nb::module_& m) {
    nb::class_<mdl::Factor>(m, "Factor", "A local operator located at a site, with the site's basis.")
        .def("__init__", [](mdl::Factor* self, int site, mdl::LocalOp op, mdl::SiteBasisPtr basis) {
            new (self) mdl::Factor{site, std::move(op), std::move(basis)};
        }, nb::arg("site"), nb::arg("op"), nb::arg("basis"))
        .def_rw("site", &mdl::Factor::site)
        .def_rw("op", &mdl::Factor::op)
        .def_rw("basis", &mdl::Factor::basis)
        .def("__repr__", [](mdl::Factor const& f) { return f.op.name + "@" + std::to_string(f.site); });
    nb::class_<mdl::Term>(m, "Term", "coeff * (op@i1 * op@i2 * ...) with the accumulated fermionic swap sign.")
        .def("__init__", [](mdl::Term* self, Complex coeff, std::vector<mdl::Factor> factors, int sign) {
            new (self) mdl::Term{coeff, std::move(factors), sign};
        }, nb::arg("coeff"), nb::arg("factors"), nb::arg("sign") = 1)
        .def_rw("coeff", &mdl::Term::coeff)
        .def_rw("factors", &mdl::Term::factors)
        .def_rw("sign", &mdl::Term::sign)
        .def_prop_ro("effective_coeff", &mdl::Term::effective_coeff)
        .def("qn_delta", &mdl::Term::qn_delta, "Total charge change for every charge all factors respect.")
        .def("__repr__", &term_repr);
    m.def("canonicalise_term", [](mdl::Term const& term, bool apply_swap_sign, std::optional<double> identity_tolerance) {
        return identity_tolerance ? mdl::canonicalise_term(term, apply_swap_sign, *identity_tolerance)
                                  : mdl::canonicalise_term(term, apply_swap_sign);
    }, nb::arg("term"), nb::arg("apply_swap_sign") = true, nb::arg("identity_tolerance") = nb::none());

    nb::class_<mdl::TermSum>(m, "TermSum",
        "A formal linear combination of terms with distributive arithmetic, Hermitian conjugation and "
        "canonicalisation. Build one with located(basis, name, site) and the +, -, * operators.")
        .def(nb::init<>())
        .def(nb::init<std::vector<mdl::Term>>(), nb::arg("terms"))
        .def_static("zero", &mdl::TermSum::zero)
        .def_static("located", &mdl::TermSum::located, nb::arg("basis"), nb::arg("op"), nb::arg("site"),
                    "The single operator `op` of `basis` at `site`.")
        .def_prop_ro("terms", &mdl::TermSum::terms)
        .def_prop_ro("empty", &mdl::TermSum::empty)
        .def("__len__", [](mdl::TermSum const& s) { return s.terms().size(); })
        .def("__add__", [](mdl::TermSum const& a, mdl::TermSum const& b) { return a + b; }, nb::is_operator())
        .def("__sub__", [](mdl::TermSum const& a, mdl::TermSum const& b) { return a - b; }, nb::is_operator())
        // In place: `H += term` appends without copying H, so accumulating many terms is linear.
        // Like list +=, this changes the object every other name bound to it sees.
        .def(nb::self += nb::self)
        .def(nb::self -= nb::self)
        .def("__neg__", [](mdl::TermSum const& a) { return -a; })
        .def("__mul__", [](mdl::TermSum const& a, mdl::TermSum const& b) { return a * b; }, nb::is_operator())
        .def("__mul__", [](mdl::TermSum const& a, Complex scalar) { return a * scalar; }, nb::is_operator())
        .def("__rmul__", [](mdl::TermSum const& a, Complex scalar) { return scalar * a; }, nb::is_operator())
        .def("__truediv__", [](mdl::TermSum const& a, Complex scalar) { return a * (1. / scalar); }, nb::is_operator())
        .def("dag", &mdl::TermSum::dag, "Hermitian conjugate.")
        .def("add_hc", &mdl::TermSum::add_hc, "The sum plus its Hermitian conjugate.")
        .def("canonicalise", [](mdl::TermSum const& s, bool apply_swap_sign, std::optional<double> identity_tolerance) {
            return identity_tolerance ? s.canonicalise(apply_swap_sign, *identity_tolerance) : s.canonicalise(apply_swap_sign);
        }, nb::arg("apply_swap_sign") = true, nb::arg("identity_tolerance") = nb::none(),
           "Canonicalise every term (site order, merged same-site factors) and combine like terms.")
        .def("split_factors_by_qn_delta", &mdl::TermSum::split_factors_by_qn_delta)
        .def("__repr__", [](mdl::TermSum const& s) {
            if (s.terms().size() > 4) return "TermSum(" + std::to_string(s.terms().size()) + " terms)";
            std::string out = "TermSum(";
            for (std::size_t i = 0; i < s.terms().size(); ++i) out += (i ? " + " : "") + term_repr(s.terms()[i]);
            return out + ")";
        });
    m.def("located", &mdl::TermSum::located, nb::arg("basis"), nb::arg("op"), nb::arg("site"),
          "The single operator `op` of `basis` at `site`, as a TermSum.");
}

// -- compiled operators --------------------------------------------------------

std::size_t factor_dim(mdl::LocalFactor const& factor) {
    auto const dim = static_cast<std::size_t>(std::llround(std::sqrt(static_cast<double>(factor.matrix.size()))));
    if (dim * dim != factor.matrix.size()) throw std::logic_error("local factor matrix is not square");
    return dim;
}

void bind_compiled(nb::module_& m) {
    nb::class_<mdl::LocalFactor>(m, "LocalFactor", "One located row-major matrix (row = output, column = input).")
        .def("__init__", [](mdl::LocalFactor* self, int site, nb::handle matrix) {
            std::size_t dim = 0;
            auto values = square_matrix_from(matrix, dim, "LocalFactor matrix");
            new (self) mdl::LocalFactor{site, std::move(values)};
        }, nb::arg("site"), nb::arg("matrix"))
        .def_rw("site", &mdl::LocalFactor::site)
        .def_prop_ro("dim", &factor_dim)
        .def_prop_ro("matrix", [](mdl::LocalFactor const& f) { return square_matrix(f.matrix, factor_dim(f)); })
        .def("__eq__", [](mdl::LocalFactor const& a, mdl::LocalFactor const& b) { return a == b; }, nb::is_operator())
        .def("__repr__", [](mdl::LocalFactor const& f) {
            return "LocalFactor(site=" + std::to_string(f.site) + ", dim=" + std::to_string(factor_dim(f)) + ")";
        });
    nb::class_<mdl::LocalTerm>(m, "LocalTerm",
        "coefficient * product of local factors * parity strings on the jordan_wigner sites.")
        .def("__init__", [](mdl::LocalTerm* self, Complex coefficient, std::vector<mdl::LocalFactor> factors,
                            std::vector<int> jordan_wigner) {
            new (self) mdl::LocalTerm{coefficient, std::move(factors), std::move(jordan_wigner)};
        }, nb::arg("coefficient"), nb::arg("factors"), nb::arg("jordan_wigner") = std::vector<int>{})
        .def_rw("coefficient", &mdl::LocalTerm::coefficient)
        .def_rw("factors", &mdl::LocalTerm::factors)
        .def_rw("jordan_wigner", &mdl::LocalTerm::jordan_wigner)
        .def_prop_ro("sites", [](mdl::LocalTerm const& t) {
            std::vector<int> sites;
            for (auto const& f : t.factors) sites.push_back(f.site);
            return sites;
        }, "The factor sites, in factor order.")
        .def("__eq__", [](mdl::LocalTerm const& a, mdl::LocalTerm const& b) { return a == b; }, nb::is_operator())
        .def("__repr__", [](mdl::LocalTerm const& t) {
            std::string out = "LocalTerm(" + format_complex(t.coefficient) + ", sites=[";
            for (std::size_t i = 0; i < t.factors.size(); ++i) out += (i ? ", " : "") + std::to_string(t.factors[i].site);
            out += "]";
            if (!t.jordan_wigner.empty()) out += ", jordan_wigner=" + join(t.jordan_wigner);
            return out + ")";
        });
    nb::class_<mdl::CompiledOperator>(m, "CompiledOperator",
        "Owned numeric operator data: the local dimension of every site and the terms. Any consumer of "
        "ALPS::ed or ALPS::classic_dmrg takes one; non-Hermitian observables are allowed.")
        .def(nb::init<>())
        .def("__init__", [](mdl::CompiledOperator* self, std::vector<int> site_dims, std::vector<mdl::LocalTerm> terms) {
            new (self) mdl::CompiledOperator{std::move(site_dims), std::move(terms)};
        }, nb::arg("site_dims"), nb::arg("terms"))
        .def_rw("site_dims", &mdl::CompiledOperator::site_dims)
        .def_rw("terms", &mdl::CompiledOperator::terms)
        .def_prop_ro("num_sites", [](mdl::CompiledOperator const& op) { return op.site_dims.size(); })
        .def("__len__", [](mdl::CompiledOperator const& op) { return op.terms.size(); })
        .def("__eq__", [](mdl::CompiledOperator const& a, mdl::CompiledOperator const& b) { return a == b; }, nb::is_operator())
        .def("__repr__", [](mdl::CompiledOperator const& op) {
            return "CompiledOperator(" + std::to_string(op.site_dims.size()) + " sites, "
                 + std::to_string(op.terms.size()) + " terms)";
        });
    nb::class_<mdl::ModelEdge>(m, "ModelEdge", "A bond of the tiled lattice as the compiled model records it.")
        .def("__init__", [](mdl::ModelEdge* self, int source, int target, int bond_type) {
            new (self) mdl::ModelEdge{source, target, bond_type};
        }, nb::arg("source"), nb::arg("target"), nb::arg("bond_type") = 0)
        .def_rw("source", &mdl::ModelEdge::src)
        .def_rw("target", &mdl::ModelEdge::dst)
        .def_rw("bond_type", &mdl::ModelEdge::bond_type)
        .def("__eq__", [](mdl::ModelEdge const& a, mdl::ModelEdge const& b) { return a == b; }, nb::is_operator())
        .def("__repr__", [](mdl::ModelEdge const& e) {
            return "ModelEdge(" + std::to_string(e.src) + ", " + std::to_string(e.dst) + ", bond_type="
                 + std::to_string(e.bond_type) + ")";
        });
    nb::class_<mdl::CompiledModel>(m, "CompiledModel",
        "A compiled operator with the lattice metadata: site types, edges and per-site-type charge values.")
        .def(nb::init<>())
        .def_rw("op", &mdl::CompiledModel::op)
        .def_rw("site_types", &mdl::CompiledModel::site_types)
        .def_rw("edges", &mdl::CompiledModel::edges)
        .def_rw("spin_basis", &mdl::CompiledModel::spin_basis)
        .def_rw("number_basis", &mdl::CompiledModel::number_basis)
        .def_rw("diag_observables", &mdl::CompiledModel::diag_observables)
        .def("__repr__", [](mdl::CompiledModel const& model) {
            return "CompiledModel(" + std::to_string(model.op.site_dims.size()) + " sites, "
                 + std::to_string(model.op.terms.size()) + " terms, " + std::to_string(model.edges.size()) + " edges)";
        });
}

// -- Hamiltonians ---------------------------------------------------------------

mdl::TermSum term_sum_result(nb::object result, char const* what) {
    if (!nb::isinstance<mdl::TermSum>(result))
        throw nb::type_error((std::string(what) + " must return a TermSum, not "
                              + alps::python::type_name(result)).c_str());
    return nb::cast<mdl::TermSum>(result);
}

void bind_hamiltonian(nb::module_& m) {
    nb::class_<mdl::Hamiltonian>(m, "Hamiltonian",
        "A lattice, one site basis (or one per site type) and the terms added to it. Bond and site "
        "specifications are callables (source, target) -> TermSum and (site) -> TermSum, called on every "
        "bond or site when they are added.")
        .def("__init__", [](mdl::Hamiltonian* self, lat::Lattice lattice, nb::handle sites, nb::handle convention,
                            std::vector<std::string> conserved) {
            auto const conv = convention_of(convention);
            if (nb::isinstance<nb::dict>(sites)) {
                auto by_type = nb::cast<std::map<int, mdl::SiteBasisPtr>>(sites);
                new (self) mdl::Hamiltonian(std::move(lattice), std::move(by_type), conv, std::move(conserved));
            } else {
                new (self) mdl::Hamiltonian(std::move(lattice), nb::cast<mdl::SiteBasisPtr>(sites), conv,
                                            std::move(conserved));
            }
        }, nb::arg("lattice"), nb::arg("sites"), nb::arg("fermion_convention") = "strict",
           nb::arg("conserved_quantum_numbers") = std::vector<std::string>{},
           "`sites` is one SiteBasis for every site or a dict site_type -> SiteBasis.")
        // Python specifications are evaluated here instead of being stored as the
        // C++ model's std::function templates: a callable kept alive from C++ holds
        // its module globals, a reference cycle the garbage collector cannot see
        // (nanobind reports the objects as leaked at exit), and compile() would have
        // to re-enter Python. The expansion is the same as Hamiltonian::expand's.
        .def("add_bond", [](mdl::Hamiltonian& model, Complex coeff, nb::callable spec, std::optional<int> bond_type,
                            bool add_hc) -> mdl::Hamiltonian& {
            for (lat::Bond const& bond : model.lat().bonds()) {
                if (bond_type && bond.bond_type != *bond_type) continue;
                model.add(term_sum_result(spec(bond.source, bond.target), "a bond specification") * coeff, add_hc);
            }
            return model;
        }, nb::arg("coeff"), nb::arg("spec"), nb::arg("bond_type") = nb::none(), nb::arg("add_hc") = false,
           nb::rv_policy::reference_internal,
           "Add coeff * spec(source, target) on every bond (of `bond_type` when given); spec is called now.")
        .def("add_site", [](mdl::Hamiltonian& model, Complex coeff, nb::callable spec, std::optional<int> site_type,
                            bool add_hc) -> mdl::Hamiltonian& {
            for (lat::Site const& site : model.lat().sites()) {
                if (site_type && site.site_type != *site_type) continue;
                model.add(term_sum_result(spec(site.index), "a site specification") * coeff, add_hc);
            }
            return model;
        }, nb::arg("coeff"), nb::arg("spec"), nb::arg("site_type") = nb::none(), nb::arg("add_hc") = false,
           nb::rv_policy::reference_internal,
           "Add coeff * spec(site) on every site (of `site_type` when given); spec is called now.")
        .def("add", [](mdl::Hamiltonian& model, mdl::TermSum term, bool add_hc) -> mdl::Hamiltonian& {
            return model.add(std::move(term), add_hc);
        }, nb::arg("term"), nb::arg("add_hc") = false, nb::rv_policy::reference_internal,
           "Add an explicit located TermSum.")
        .def_prop_ro("lattice", &mdl::Hamiltonian::lat)
        .def_prop_ro("num_sites", &mdl::Hamiltonian::num_sites)
        .def("site_for", &mdl::Hamiltonian::site_for, nb::arg("site"), "The SiteBasis of a lattice site.")
        .def_prop_ro("site_dims", &mdl::Hamiltonian::site_dims)
        .def_prop_ro("conserved_quantum_numbers", &mdl::Hamiltonian::conserved_qns)
        .def_prop_ro("fermion_convention", &mdl::Hamiltonian::fermion_convention)
        .def("expand", &mdl::Hamiltonian::expand, "Every term of the model in one TermSum (not canonicalised).")
        .def("to_local_decomposition", [](mdl::Hamiltonian const& model) {
            nb::gil_scoped_release release;
            return model.to_local_decomposition();
        }, "The Monte Carlo descriptor: charge-split local terms with the lattice metadata.")
        .def("operator_decomposition", [](mdl::Hamiltonian const& model) {
            nb::gil_scoped_release release;
            return model.operator_decomposition();
        }, "The full Hamiltonian lowered without the symmetry declaration.")
        .def("operator_decomposition", [](mdl::Hamiltonian const& model, mdl::TermSum const& op) {
            nb::gil_scoped_release release;
            return model.operator_decomposition(op);
        }, nb::arg("op"), "An arbitrary (charge-changing, non-Hermitian) observable lowered on this model.")
        .def("compile", [](mdl::Hamiltonian const& model) {
            nb::gil_scoped_release release;
            return mdl::compile_hamiltonian(model);
        }, "compile_hamiltonian(self).")
        .def("physical_basis", [](mdl::Hamiltonian const& model) { return mdl::physical_basis(model); })
        .def("__repr__", [](mdl::Hamiltonian const& model) {
            return "Hamiltonian(" + std::to_string(model.num_sites()) + " sites, conserved="
                 + join(model.conserved_qns()) + ")";
        });

    m.def("compile_hamiltonian", [](mdl::Hamiltonian const& model) {
        nb::gil_scoped_release release;
        return mdl::compile_hamiltonian(model);
    }, nb::arg("model"), "Expand and compile a model under its declared charges and convention (a CompiledModel).");
    m.def("compile_hamiltonian", [](mdl::TermSum const& terms, std::vector<mdl::SiteBasisPtr> const& sites,
                                    std::vector<std::string> const& conserved, nb::handle convention) {
        auto const conv = convention_of(convention);
        nb::gil_scoped_release release;
        return mdl::compile_hamiltonian(terms, sites, conserved, conv);
    }, nb::arg("terms"), nb::arg("sites"), nb::arg("conserved_quantum_numbers") = std::vector<std::string>{},
       nb::arg("fermion_convention") = "strict",
       "Compile located terms on ordered sites, validating the declared charges (a CompiledOperator).");
    m.def("compile_operator", [](mdl::TermSum const& terms, std::vector<mdl::SiteBasisPtr> const& sites,
                                 nb::handle convention) {
        auto const conv = convention_of(convention);
        nb::gil_scoped_release release;
        return mdl::compile_operator(terms, sites, conv);
    }, nb::arg("terms"), nb::arg("sites"), nb::arg("fermion_convention") = "strict",
       "Lower an arbitrary observable without Hamiltonian symmetry requirements.");
    m.def("has_sign_problem", [](mdl::Hamiltonian const& model) {
        nb::gil_scoped_release release;
        return mdl::has_sign_problem(model);
    }, nb::arg("model"),
       "The legacy algebraic QMC sign check (sublattice-sign rule and frustration test) on the summed one- and "
       "two-site terms per site and per site pair; raises ValueError for Jordan-Wigner strings, terms on more "
       "than two sites and imaginary matrix elements.");
    m.def("has_sign_problem", [](mdl::CompiledOperator const& op) {
        nb::gil_scoped_release release;
        return mdl::has_sign_problem(op);
    }, nb::arg("op"), "The same for a compiled operator.");
    m.def("compile_operator", [](mdl::TermSum const& terms, mdl::PhysicalBasis const& physical) {
        nb::gil_scoped_release release;
        return mdl::compile_operator(terms, physical);
    }, nb::arg("terms"), nb::arg("physical_basis"),
       "The same on a complete PhysicalBasis (a problem's, a catalog model's): each factor's site must match it.");
    m.def("compile_local_terms", [](mdl::TermSum const& terms, std::vector<mdl::SiteBasisPtr> const& sites,
                                    std::vector<std::string> const& conserved, nb::handle convention) {
        auto const conv = convention_of(convention);
        nb::gil_scoped_release release;
        return mdl::compile_local_terms(terms, sites, conserved, conv);
    }, nb::arg("terms"), nb::arg("sites"), nb::arg("conserved_quantum_numbers") = std::vector<std::string>{},
       nb::arg("fermion_convention") = "strict");
    m.def("compile_local_operator", [](mdl::TermSum const& terms, std::vector<mdl::SiteBasisPtr> const& sites,
                                       std::vector<std::string> const& conserved, nb::handle convention) {
        auto const conv = convention_of(convention);
        nb::gil_scoped_release release;
        return mdl::compile_local_operator(terms, sites, conserved, conv);
    }, nb::arg("terms"), nb::arg("sites"), nb::arg("conserved_quantum_numbers") = std::vector<std::string>{},
       nb::arg("fermion_convention") = "strict");
    m.def("validate_hermiticity", [](mdl::Hamiltonian const& model, double tolerance, std::size_t max_expanded_terms) {
        nb::gil_scoped_release release;
        mdl::validate_hermiticity(model, tolerance, max_expanded_terms);
    }, nb::arg("model"), nb::arg("tolerance") = 1e-11, nb::arg("max_expanded_terms") = 100000,
       "Raise ValueError unless the model is Hermitian within the relative tolerance.");
    m.def("validate_hermiticity", [](mdl::TermSum const& terms, std::vector<mdl::SiteBasisPtr> const& sites,
                                     nb::handle convention, double tolerance, std::size_t max_expanded_terms) {
        auto const conv = convention_of(convention);
        nb::gil_scoped_release release;
        mdl::validate_hermiticity(terms, sites, conv, tolerance, max_expanded_terms);
    }, nb::arg("terms"), nb::arg("sites"), nb::arg("fermion_convention") = "strict", nb::arg("tolerance") = 1e-11,
       nb::arg("max_expanded_terms") = 100000);
    m.def("jw_chain_sites", &mdl::jw_chain_sites, nb::arg("factors"), nb::arg("num_sites"),
          "Sites needing a Jordan-Wigner parity factor for canonically ordered factors.");

    // -- built-in models --
    m.def("heisenberg", [](lat::Lattice const& l, double J, double h, double S, std::optional<int> bond_type) {
        return mdl::heisenberg(l, J, h, S, bond_type);
    }, nb::arg("lattice"), nb::arg("J") = 1.0, nb::arg("h") = 0.0, nb::arg("S") = 0.5, nb::arg("bond_type") = nb::none(),
       "J sum S_i.S_j - h sum Sz_i; charge Sz.");
    m.def("ising", [](lat::Lattice const& l, double J, double h, double S, std::optional<int> bond_type) {
        return mdl::ising(l, J, h, S, bond_type);
    }, nb::arg("lattice"), nb::arg("J") = 1.0, nb::arg("h") = 0.0, nb::arg("S") = 0.5, nb::arg("bond_type") = nb::none());
    m.def("xxz", [](lat::Lattice const& l, double Jxy, double Jz, double h, double S, std::optional<int> bond_type) {
        return mdl::xxz(l, Jxy, Jz, h, S, bond_type);
    }, nb::arg("lattice"), nb::arg("Jxy") = 1.0, nb::arg("Jz") = 1.0, nb::arg("h") = 0.0, nb::arg("S") = 0.5,
       nb::arg("bond_type") = nb::none());
    m.def("xy", [](lat::Lattice const& l, double J, double S, std::optional<int> bond_type) {
        return mdl::xy(l, J, S, bond_type);
    }, nb::arg("lattice"), nb::arg("J") = 1.0, nb::arg("S") = 0.5, nb::arg("bond_type") = nb::none());
    m.def("transverse_field_ising", [](lat::Lattice const& l, double J, double Gamma, double h, double S,
                                       std::optional<int> bond_type) {
        return mdl::transverse_field_ising(l, J, Gamma, h, S, bond_type);
    }, nb::arg("lattice"), nb::arg("J") = 1.0, nb::arg("Gamma") = 1.0, nb::arg("h") = 0.0, nb::arg("S") = 0.5,
       nb::arg("bond_type") = nb::none());
    m.def("heisenberg_d", [](lat::Lattice const& l, double J, double D, double h, double S, std::optional<int> bond_type) {
        return mdl::heisenberg_d(l, J, D, h, S, bond_type);
    }, nb::arg("lattice"), nb::arg("J") = 1.0, nb::arg("D") = 0.0, nb::arg("h") = 0.0, nb::arg("S") = 1.0,
       nb::arg("bond_type") = nb::none(), "Heisenberg with single-ion anisotropy D sum (Sz_i)^2.");
    m.def("biquadratic", [](lat::Lattice const& l, double J, double K, double S, std::optional<int> bond_type) {
        return mdl::biquadratic(l, J, K, S, bond_type);
    }, nb::arg("lattice"), nb::arg("J") = 1.0, nb::arg("K") = 0.0, nb::arg("S") = 1.0, nb::arg("bond_type") = nb::none());
    m.def("j1j2", [](lat::Lattice const& l, double J1, double J2, double S) { return mdl::j1j2(l, J1, J2, S); },
          nb::arg("lattice"), nb::arg("J1") = 1.0, nb::arg("J2") = 0.5, nb::arg("S") = 0.5,
          "J1 on bond type 0, J2 on bond type 1 (nnn_chain, frustrated_square).");
    m.def("hubbard", [](lat::Lattice const& l, double t, double U, double mu, std::optional<int> bond_type) {
        return mdl::hubbard(l, t, U, mu, bond_type);
    }, nb::arg("lattice"), nb::arg("t") = 1.0, nb::arg("U") = 1.0, nb::arg("mu") = 0.0, nb::arg("bond_type") = nb::none(),
       "-t hopping + U n_up n_dn - mu n; charges N and Sz.");
    m.def("spinless_fermion_model", [](lat::Lattice const& l, double t, double mu, double V, std::optional<int> bond_type) {
        return mdl::spinless_fermion_model(l, t, mu, V, bond_type);
    }, nb::arg("lattice"), nb::arg("t") = 1.0, nb::arg("mu") = 0.0, nb::arg("V") = 0.0, nb::arg("bond_type") = nb::none());
    m.def("t_j", [](lat::Lattice const& l, double t, double J, std::optional<int> bond_type) {
        return mdl::t_j(l, t, J, bond_type);
    }, nb::arg("lattice"), nb::arg("t") = 1.0, nb::arg("J") = 0.4, nb::arg("bond_type") = nb::none());
    m.def("kondo_lattice", [](lat::Lattice const& l, double t, double J, double mu) { return mdl::kondo_lattice(l, t, J, mu); },
          nb::arg("lattice"), nb::arg("t") = 1.0, nb::arg("J") = 1.0, nb::arg("mu") = 0.0,
          "Site type 0 itinerant spinful fermions, site type 1 localized spins; hopping on bond type 0, "
          "Kondo exchange on bond type 1.");
    m.def("hardcore_boson_model", [](lat::Lattice const& l, double t, double mu, double V, std::optional<int> bond_type) {
        return mdl::hardcore_boson_model(l, t, mu, V, bond_type);
    }, nb::arg("lattice"), nb::arg("t") = 1.0, nb::arg("mu") = 0.0, nb::arg("V") = 0.0, nb::arg("bond_type") = nb::none());
    m.def("bose_hubbard", [](lat::Lattice const& l, double t, double U, double mu, double V, int n_max,
                             std::optional<int> bond_type) {
        return mdl::bose_hubbard(l, t, U, mu, V, n_max, bond_type);
    }, nb::arg("lattice"), nb::arg("t") = 1.0, nb::arg("U") = 1.0, nb::arg("mu") = 0.0, nb::arg("V") = 0.0,
       nb::arg("n_max") = 3, nb::arg("bond_type") = nb::none());
    m.def("ab_initio", [](nb::handle h1, nb::handle eri, double tol) {
        auto const one = alps::python::as_contiguous<double>(h1);
        auto const two = alps::python::as_contiguous<double>(eri);
        if (one.ndim() != 2 || one.shape(0) != one.shape(1))
            throw nb::value_error("h1 must be a square two-dimensional array");
        std::size_t const n = one.shape(0);
        if (two.ndim() != 4 || two.shape(0) != n || two.shape(1) != n || two.shape(2) != n || two.shape(3) != n)
            throw nb::value_error("eri must be an (n, n, n, n) array matching h1");
        std::vector<double> h(one.data(), one.data() + n * n), g(two.data(), two.data() + n * n * n * n);
        nb::gil_scoped_release release;
        return mdl::ab_initio(h, g, static_cast<int>(n), tol);
    }, nb::arg("h1"), nb::arg("eri"), nb::arg("tol") = 1e-12,
       "The molecular Hamiltonian of real one- and two-electron integrals (chemist's notation): one "
       "spinful_fermion() site per orbital, strict convention, N and Sz conserved; integrals that are zero or "
       "below `tol` in magnitude are dropped.");
}

// -- qubits and Pauli strings -------------------------------------------------------

mdl::PauliExportOptions pauli_options(double tolerance, std::size_t max_terms, bool require_hermitian, bool require_even,
                                      bool require_conservation, std::size_t max_local_matrix_bytes,
                                      std::optional<double> validation_tolerance) {
    mdl::PauliExportOptions options;
    options.tolerance = tolerance;
    options.max_terms = max_terms;
    options.require_hermitian = require_hermitian;
    options.require_even = require_even;
    options.require_conservation = require_conservation;
    options.max_local_matrix_bytes = max_local_matrix_bytes;
    options.validation_tolerance = validation_tolerance;
    return options;
}

#define PYALPS_PAULI_OPTION_ARGS                                                                               \
    nb::arg("tolerance") = 1e-12, nb::arg("max_terms") = std::size_t{100000}, nb::arg("require_hermitian") = true, \
    nb::arg("require_even") = true, nb::arg("require_conservation") = true,                                    \
    nb::arg("max_local_matrix_bytes") = mdl::default_local_matrix_bytes, nb::arg("validation_tolerance") = nb::none()

void bind_qubits(nb::module_& m) {
    nb::enum_<mdl::QubitEncodingScheme>(m, "QubitEncodingScheme")
        .value("binary", mdl::QubitEncodingScheme::binary)
        .value("gray", mdl::QubitEncodingScheme::gray);
    m.attr("default_local_matrix_bytes") = mdl::default_local_matrix_bytes;
    m.attr("default_state_amplitudes") = mdl::default_state_amplitudes;

    nb::class_<mdl::QubitEncoding>(m, "QubitEncoding",
        "The qubit register of a PhysicalBasis: every site occupies consecutive qubits, site 0 the least "
        "significant ones, and its local index is written as a binary or Gray code word.")
        .def("__init__", [](mdl::QubitEncoding* self, mdl::PhysicalBasis const& basis, nb::handle scheme,
                            std::size_t max_local_matrix_bytes) {
            new (self) mdl::QubitEncoding(basis, scheme_of(scheme), max_local_matrix_bytes);
        }, nb::arg("basis"), nb::arg("scheme") = "binary", nb::arg("max_local_matrix_bytes") = mdl::default_local_matrix_bytes)
        .def_prop_ro("num_qubits", &mdl::QubitEncoding::num_qubits)
        .def_prop_ro("scheme", &mdl::QubitEncoding::scheme)
        .def_prop_ro("site_qubits", &mdl::QubitEncoding::site_qubits)
        .def_prop_ro("site_offsets", &mdl::QubitEncoding::site_offsets)
        .def_prop_ro("site_dimensions", &mdl::QubitEncoding::site_dimensions)
        .def_prop_ro("sites", &mdl::QubitEncoding::sites)
        .def("codeword", &mdl::QubitEncoding::codeword, nb::arg("site"), nb::arg("local_index"))
        .def("encode_basis_bits", &mdl::QubitEncoding::encode_basis_bits, nb::arg("local_indices"))
        .def("decode_basis_bits", &mdl::QubitEncoding::decode_basis_bits, nb::arg("bits"),
             "The local indices of a bit string, None for an unphysical code word.")
        .def("state_index", &mdl::QubitEncoding::state_index, nb::arg("local_indices"),
             "Register index of a physical product state (at most 63 qubits).")
        .def("require_same_sites", &mdl::QubitEncoding::require_same_sites, nb::arg("basis"))
        .def("encode_statevector", [](mdl::QubitEncoding const& e, nb::handle physical, std::size_t max_amplitudes) {
            auto const state = vector_from<Complex>(physical, "physical state");
            return vector_to_numpy(e.encode_statevector(state, max_amplitudes));
        }, nb::arg("physical"), nb::arg("max_amplitudes") = mdl::default_state_amplitudes,
           "The register amplitudes of a physical state vector (site 0 most significant).")
        .def("decode_statevector", [](mdl::QubitEncoding const& e, nb::handle encoded, std::size_t max_amplitudes) {
            auto const state = vector_from<Complex>(encoded, "encoded state");
            auto decoded = e.decode_statevector(state, max_amplitudes);
            return std::make_pair(vector_to_numpy(std::move(decoded.amplitudes)), decoded.discarded_norm_squared);
        }, nb::arg("encoded"), nb::arg("max_amplitudes") = mdl::default_state_amplitudes,
           "(physical amplitudes, discarded squared norm) of a register state; no renormalisation.")
        .def("physical_weight", [](mdl::QubitEncoding const& e, nb::handle encoded, std::size_t max_amplitudes) {
            return e.physical_weight(vector_from<Complex>(encoded, "encoded state"), max_amplitudes);
        }, nb::arg("encoded"), nb::arg("max_amplitudes") = mdl::default_state_amplitudes)
        .def("__repr__", [](mdl::QubitEncoding const& e) {
            return "QubitEncoding(" + std::to_string(e.num_qubits()) + " qubits, "
                 + (e.scheme() == mdl::QubitEncodingScheme::gray ? "gray" : "binary") + ")";
        });
    m.def("make_qubit_encoding", [](mdl::Hamiltonian const& model, nb::handle scheme, std::size_t max_local_matrix_bytes) {
        return mdl::make_qubit_encoding(model, scheme_of(scheme), max_local_matrix_bytes);
    }, nb::arg("model"), nb::arg("scheme") = "binary", nb::arg("max_local_matrix_bytes") = mdl::default_local_matrix_bytes,
       "The encoding of physical_basis(model).");

    nb::class_<mdl::PauliOperator>(m, "PauliOperator",
        "A sum of Pauli strings on a qubit register; `terms` maps strings (qubit 0 at the right end) to coefficients.")
        .def_ro("num_qubits", &mdl::PauliOperator::num_qubits)
        .def_ro("site_qubits", &mdl::PauliOperator::site_qubits)
        .def_ro("fermion_mode_qubits", &mdl::PauliOperator::fermion_mode_qubits,
                "Site-major order of the fermionic modes (down, up per spinful site); empty for bosons.")
        .def_ro("terms", &mdl::PauliOperator::terms)
        .def_ro("discarded_l1", &mdl::PauliOperator::discarded_l1, "Operator-norm bound on the pruning error.")
        .def_ro("encoding", &mdl::PauliOperator::encoding)
        .def_ro("local_decompositions", &mdl::PauliOperator::local_decompositions)
        .def_ro("local_cache_hits", &mdl::PauliOperator::local_cache_hits)
        .def("coefficient_l1", &mdl::PauliOperator::coefficient_l1, nb::arg("include_identity") = true)
        .def("state_index", &mdl::PauliOperator::state_index, nb::arg("local_indices"))
        .def("__len__", [](mdl::PauliOperator const& p) { return p.terms.size(); })
        .def("__repr__", [](mdl::PauliOperator const& p) {
            return "PauliOperator(" + std::to_string(p.num_qubits) + " qubits, " + std::to_string(p.terms.size()) + " strings)";
        });
    m.def("export_pauli", [](mdl::CompiledOperator const& op, mdl::PhysicalBasis const& basis,
                             std::optional<mdl::QubitEncoding> const& encoding, double tolerance, std::size_t max_terms,
                             bool require_hermitian, bool require_even, bool require_conservation,
                             std::size_t max_local_matrix_bytes, std::optional<double> validation_tolerance) {
        auto const options = pauli_options(tolerance, max_terms, require_hermitian, require_even, require_conservation,
                                           max_local_matrix_bytes, validation_tolerance);
        nb::gil_scoped_release release;
        return encoding ? mdl::export_pauli(op, basis, *encoding, options) : mdl::export_pauli(op, basis, options);
    }, nb::arg("op"), nb::arg("basis"), nb::arg("encoding") = nb::none(), PYALPS_PAULI_OPTION_ARGS,
       "Pauli strings of a compiled operator on its physical basis (binary encoding unless one is given). "
       "Switch the three requirements off for an arbitrary observable.");
    m.def("export_pauli", [](mdl::Hamiltonian const& model, std::optional<mdl::QubitEncoding> const& encoding,
                             double tolerance, std::size_t max_terms, bool require_hermitian, bool require_even,
                             bool require_conservation, std::size_t max_local_matrix_bytes,
                             std::optional<double> validation_tolerance) {
        auto const options = pauli_options(tolerance, max_terms, require_hermitian, require_even, require_conservation,
                                           max_local_matrix_bytes, validation_tolerance);
        nb::gil_scoped_release release;
        return encoding ? mdl::export_pauli(model, *encoding, options) : mdl::export_pauli(model, options);
    }, nb::arg("model"), nb::arg("encoding") = nb::none(), PYALPS_PAULI_OPTION_ARGS,
       "Pauli strings of compile_hamiltonian(model) on physical_basis(model).");
    m.def("local_valid_projector", [](mdl::QubitEncoding const& encoding, std::size_t site, double tolerance,
                                      std::size_t max_terms, bool require_hermitian, bool require_even,
                                      bool require_conservation, std::size_t max_local_matrix_bytes,
                                      std::optional<double> validation_tolerance) {
        auto const options = pauli_options(tolerance, max_terms, require_hermitian, require_even, require_conservation,
                                           max_local_matrix_bytes, validation_tolerance);
        nb::gil_scoped_release release;
        return mdl::local_valid_projector(encoding, site, options);
    }, nb::arg("encoding"), nb::arg("site"), PYALPS_PAULI_OPTION_ARGS,
       "The projector onto the physical code words of one site.");
    m.def("invalid_state_penalty", [](mdl::QubitEncoding const& encoding, double strength, double tolerance,
                                      std::size_t max_terms, bool require_hermitian, bool require_even,
                                      bool require_conservation, std::size_t max_local_matrix_bytes,
                                      std::optional<double> validation_tolerance) {
        auto const options = pauli_options(tolerance, max_terms, require_hermitian, require_even, require_conservation,
                                           max_local_matrix_bytes, validation_tolerance);
        nb::gil_scoped_release release;
        return mdl::invalid_state_penalty(encoding, strength, options);
    }, nb::arg("encoding"), nb::arg("strength"), PYALPS_PAULI_OPTION_ARGS,
       "strength * sum_i (I - P_i) over the sites' valid-state projectors.");
}

// -- the XML catalog importer ---------------------------------------------------------

using Catalog = mdl::CatalogModel<>;

void bind_catalog(nb::module_& m) {
    nb::class_<Catalog>(m, "CatalogModel",
        "An XML catalog model (MODEL, LATTICE and their parameters) resolved through the legacy lattice and "
        "model libraries and compiled to CompiledOperator values; the catalogs come from ALPS_XML_PATH or "
        "LATTICE_LIBRARY/MODEL_LIBRARY. The second constructor runs a catalog model on a kit lattice.")
        .def("__init__", [](Catalog* self, nb::dict const& parameters) {
            auto const p = parameters_from_dict(parameters);
            nb::gil_scoped_release release;
            new (self) Catalog(p);
        }, nb::arg("parameters"))
        .def("__init__", [](Catalog* self, lat::Lattice const& lattice, nb::dict const& parameters) {
            auto const p = parameters_from_dict(parameters);
            nb::gil_scoped_release release;
            new (self) Catalog(lattice, p);
        }, nb::arg("lattice"), nb::arg("parameters"))
        .def("compile_hamiltonian", [](Catalog const& model) {
            nb::gil_scoped_release release;
            return model.compile_hamiltonian();
        }, "The catalog Hamiltonian as a CompiledOperator (no Hermiticity or conservation check).")
        .def("compile_operator", [](Catalog const& model, std::string const& name, std::optional<int> site,
                                    std::optional<int> bond) {
            if (site && bond) throw nb::value_error("give a site or a bond, not both");
            nb::gil_scoped_release release;
            if (site) {
                if (*site < 0 || *site >= static_cast<int>(model.graph().num_sites()))
                    throw std::out_of_range("site index outside the lattice");
                return model.compile_operator(name, static_cast<Catalog::site_descriptor>(*site));
            }
            if (bond) {
                if (*bond < 0 || *bond >= static_cast<int>(model.graph().num_bonds()))
                    throw std::out_of_range("bond index outside the lattice");
                auto it = model.graph().bonds().first;
                std::advance(it, *bond);
                return model.compile_operator(name, *it);
            }
            return model.compile_operator(name);
        }, nb::arg("name"), nb::arg("site") = nb::none(), nb::arg("bond") = nb::none(),
           "A named catalog site, bond or global operator (or a site operator expression such as 'Sz') "
           "summed over the lattice, or located on one site or one bond.")
        .def_prop_ro("site_dims", [](Catalog const& model) {
            std::vector<int> dims;
            for (auto const& local : model.basis()) dims.push_back(static_cast<int>(local.size()));
            return dims;
        })
        .def_prop_ro("num_sites", [](Catalog const& model) { return static_cast<int>(model.graph().num_sites()); })
        .def_prop_ro("num_bonds", [](Catalog const& model) { return static_cast<int>(model.graph().num_bonds()); })
        .def_prop_ro("site_types", [](Catalog const& model) {
            std::vector<int> types;
            for (auto it = model.graph().sites().first; it != model.graph().sites().second; ++it)
                types.push_back(static_cast<int>(model.graph().site_type(*it)));
            return types;
        })
        .def_prop_ro("edges", [](Catalog const& model) {
            std::vector<mdl::ModelEdge> edges;
            for (auto it = model.graph().bonds().first; it != model.graph().bonds().second; ++it)
                edges.push_back({static_cast<int>(model.graph().source(*it)), static_cast<int>(model.graph().target(*it)),
                                 static_cast<int>(model.graph().bond_type(*it))});
            return edges;
        }, "The lattice bonds as ModelEdge values.")
        .def_prop_ro("coordinates", [](Catalog const& model) {
            auto const& graph = model.graph();
            std::size_t const width = graph.dimension();
            std::vector<double> values;
            for (auto it = graph.sites().first; it != graph.sites().second; ++it) {
                auto const& coordinate = graph.coordinate(*it);
                for (std::size_t d = 0; d < width; ++d) values.push_back(d < coordinate.size() ? coordinate[d] : 0.);
            }
            return matrix_to_numpy(std::move(values), graph.num_sites(), width);
        })
        .def("local_parity", &Catalog::local_parity, "One +1/-1 diagonal per site.")
        .def("conserved_quantum_numbers", &Catalog::conserved_quantum_numbers,
             "CONSERVED_QUANTUMNUMBERS plus every constrained total (Sz_total, N_total, ...), sorted.")
        .def("sector", &Catalog::sector, "The evaluated constraints as totals keyed by quantum number.")
        .def("local_charge", &Catalog::local_charge, nb::arg("name"), "Per site, per local state values of a quantum number.")
        .def("physical_basis", [](Catalog const& model, std::optional<std::vector<std::string>> conserved) {
            return conserved ? model.physical_basis(*conserved) : model.physical_basis();
        }, nb::arg("conserved_quantum_numbers") = nb::none(),
           "The catalog's local spaces, parities and charges as a PhysicalBasis (strict convention); the "
           "conserved charges default to conserved_quantum_numbers().")
        .def_prop_ro("parameters", [](Catalog const& model) { return dict_from_parameters(model.parameters()); })
        .def("set_parameters", [](Catalog& model, nb::dict const& parameters) {
            model.set_parameters(parameters_from_dict(parameters));
        }, nb::arg("parameters"), "Merge evaluation and sector parameters (not the lattice or model definition).")
        .def("__repr__", [](Catalog const& model) {
            return "CatalogModel(MODEL='" + model.parameters().value_or_default("MODEL", "")
                 + "', LATTICE='" + model.parameters().value_or_default("LATTICE", "") + "', sites="
                 + std::to_string(model.graph().num_sites()) + ")";
        });
    m.def("has_sign_problem", [](Catalog const& model) {
        nb::gil_scoped_release release;
        return mdl::has_sign_problem(model.compile_hamiltonian());
    }, nb::arg("model"), "The same for a catalog model's compiled Hamiltonian.");
    m.def("compile_hamiltonian", [](Catalog const& model) {
        mdl::CompiledModel compiled;
        {
            nb::gil_scoped_release release;
            compiled.op = model.compile_hamiltonian();
        }
        auto const& graph = model.graph();
        for (auto it = graph.sites().first; it != graph.sites().second; ++it)
            compiled.site_types.push_back(static_cast<int>(graph.site_type(*it)));
        for (auto it = graph.bonds().first; it != graph.bonds().second; ++it)
            compiled.edges.push_back({static_cast<int>(graph.source(*it)), static_cast<int>(graph.target(*it)),
                                      static_cast<int>(graph.bond_type(*it))});
        return compiled;
    }, nb::arg("model"),
       "A catalog model as a CompiledModel: op, site_types and edges (no basis tables; no Hermiticity or "
       "conservation check, as CatalogModel.compile_hamiltonian).");
}

} // namespace

NB_MODULE(models_c, m) {
    m.doc() = "ALPS::models: lattices, local bases, operator algebra, Hamiltonians, compiled operators, "
              "Pauli export and the XML catalog importer";
    nb::module_ lattice = m.def_submodule("lattice",
        "Programmatic lattices: unit cells, tiling, the named lattices and the legacy XML bridge");
    bind_lattice(lattice);
    bind_sites(m);
    bind_terms(m);
    bind_compiled(m);
    bind_hamiltonian(m);
    bind_qubits(m);
    bind_catalog(m);
}
