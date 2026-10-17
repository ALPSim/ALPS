// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// pyalps.classic_dmrg: the classic (two-site, finite-system) DMRG of
// ALPS::classic_dmrg on compiled models. Options and measurement requests are
// plain objects, the solve releases the GIL and reports every iteration to an
// optional Python callback (which may raise to abort; Ctrl-C aborts too), and
// the result carries NumPy arrays.
//
// The surface follows src/alps/classic_dmrg/include/alps/dmrg/classic.hpp;
// the skilledwolf/alpsx binding (python/alps/cpp/classic_dmrg_bindings.cpp at
// d37e51227) took catalog models and returned dictionaries; its callback and
// GIL handling are kept. solve_molecular follows alpsx's binding of the same
// name (branch WIP/classic-dmrg-chemistry at d4384a7f8), returning a
// MolecularResult object with NumPy arrays instead of a dictionary; its
// observables, ci_vector, spin blocks and CI fields are added here. evolve
// (time-dependent DMRG) has no alpsx counterpart on this engine; it is bound
// as solve is.
#include "kit_support.hpp"

#include <alps/dmrg/classic.hpp>
#include <alps/dmrg/molecular.hpp>
#include <alps/models.hpp>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;
namespace classic = alps::dmrg::classic;
namespace mdl = alps::models;
using namespace pyalps::kit;

namespace {

classic::Arithmetic arithmetic_of(nb::handle value) {
    return enum_or_name<classic::Arithmetic>(value, {{"real", classic::Arithmetic::real},
                                                     {"complex", classic::Arithmetic::complex}}, "arithmetic");
}

char const* direction_name(classic::Direction value) {
    return value == classic::Direction::left_to_right ? "left_to_right" : "right_to_left";
}

std::vector<classic::SweepStates> sweeps_of(nb::handle value) {
    std::vector<classic::SweepStates> sweeps;
    for (nb::handle item : nb::borrow<nb::iterable>(value)) {
        if (nb::isinstance<classic::SweepStates>(item)) {
            sweeps.push_back(nb::cast<classic::SweepStates>(item));
        } else {
            auto const pair = nb::cast<std::pair<std::size_t, std::size_t>>(item);
            sweeps.push_back({pair.first, pair.second});
        }
    }
    return sweeps;
}

classic::Result run(mdl::CompiledOperator const* op, mdl::PhysicalBasis const* basis, mdl::Hamiltonian const* model,
                    classic::Options options, classic::Measurements const& measurements, nb::object callback,
                    std::optional<std::string> const& log_file) {
    if (!callback.is_none() && !PyCallable_Check(callback.ptr()))
        throw nb::type_error("callback must be callable or None");
    std::ofstream log;
    if (log_file) {
        log.open(*log_file);
        if (!log) throw std::runtime_error("cannot open the DMRG log file " + *log_file);
        options.log = &log;
    }
    classic::IterationCallback observe = [&](classic::Iteration const& iteration) {
        // Waiting for the engine's invocation mutex and the sweeps run without
        // the GIL; the callback takes it back and lets Ctrl-C abort the solve.
        nb::gil_scoped_acquire acquire;
        if (PyErr_CheckSignals() != 0) throw nb::python_error();
        if (!callback.is_none()) callback(iteration);
    };
    nb::gil_scoped_release release;
    if (model) return classic::solve(*model, options, measurements, observe);
    return classic::solve(*op, *basis, options, measurements, observe);
}

classic::EvolutionResult run_evolution(mdl::CompiledOperator const* initial, mdl::CompiledOperator const* final_op,
                                       mdl::PhysicalBasis const* basis, mdl::Hamiltonian const* initial_model,
                                       mdl::Hamiltonian const* final_model, classic::Evolution const& evolution,
                                       classic::Options options, classic::Measurements const& measurements,
                                       nb::object callback, std::optional<std::string> const& log_file) {
    if (!callback.is_none() && !PyCallable_Check(callback.ptr()))
        throw nb::type_error("callback must be callable or None");
    std::ofstream log;
    if (log_file) {
        log.open(*log_file);
        if (!log) throw std::runtime_error("cannot open the DMRG log file " + *log_file);
        options.log = &log;
    }
    classic::IterationCallback observe = [&](classic::Iteration const& iteration) {
        nb::gil_scoped_acquire acquire;
        if (PyErr_CheckSignals() != 0) throw nb::python_error();
        if (!callback.is_none()) callback(iteration);
    };
    nb::gil_scoped_release release;
    if (initial_model) return classic::evolve(*initial_model, *final_model, evolution, options, measurements, observe);
    return classic::evolve(*initial, *final_op, *basis, evolution, options, measurements, observe);
}

using Integrals2 = nb::ndarray<nb::numpy, double const, nb::ndim<2>, nb::c_contig>;
using Integrals4 = nb::ndarray<nb::numpy, double const, nb::ndim<4>, nb::c_contig>;

classic::MolecularResult run_molecular(Integrals2 one, Integrals4 two, std::size_t n_alpha, std::size_t n_beta,
                                       double core_energy, double integral_cutoff,
                                       std::map<std::string, mdl::TermSum> observables, bool ci_vector,
                                       classic::Options options, nb::object callback,
                                       std::optional<std::string> const& log_file) {
    auto const n = one.shape(0);
    if (one.shape(1) != n || two.shape(0) != n || two.shape(1) != n || two.shape(2) != n || two.shape(3) != n)
        throw std::invalid_argument("molecular integral shapes must be (n, n) and (n, n, n, n)");
    if (!callback.is_none() && !PyCallable_Check(callback.ptr()))
        throw nb::type_error("callback must be callable or None");
    classic::MolecularInput input;
    input.orbitals = n;
    input.alpha_electrons = n_alpha;
    input.beta_electrons = n_beta;
    input.core_energy = core_energy;
    input.integral_cutoff = integral_cutoff;
    input.observables = std::move(observables);
    input.ci_vector = ci_vector;
    if (one.size()) input.one_body.assign(one.data(), one.data() + one.size());
    if (two.size()) input.two_body.assign(two.data(), two.data() + two.size());
    std::ofstream log;
    if (log_file) {
        log.open(*log_file);
        if (!log) throw std::runtime_error("cannot open the DMRG log file " + *log_file);
        options.log = &log;
    }
    classic::IterationCallback observe = [&](classic::Iteration const& iteration) {
        nb::gil_scoped_acquire acquire;
        if (PyErr_CheckSignals() != 0) throw nb::python_error();
        if (!callback.is_none()) callback(iteration);
    };
    nb::gil_scoped_release release;
    return classic::solve_molecular(input, options, observe);
}

} // namespace

NB_MODULE(classic_dmrg_c, m) {
    m.doc() = "ALPS::classic_dmrg: classic two-site finite-system DMRG on compiled models";

    nb::enum_<classic::Arithmetic>(m, "Arithmetic")
        .value("real", classic::Arithmetic::real)
        .value("complex", classic::Arithmetic::complex);
    nb::enum_<classic::Direction>(m, "Direction")
        .value("left_to_right", classic::Direction::left_to_right)
        .value("right_to_left", classic::Direction::right_to_left);

    nb::class_<classic::SweepStates>(m, "SweepStates", "Kept states of one sweep; the right-to-left half sweep runs first.")
        .def("__init__", [](classic::SweepStates* self, std::size_t right_to_left, std::size_t left_to_right) {
            new (self) classic::SweepStates{right_to_left, left_to_right};
        }, nb::arg("right_to_left"), nb::arg("left_to_right"))
        .def_rw("right_to_left", &classic::SweepStates::right_to_left)
        .def_rw("left_to_right", &classic::SweepStates::left_to_right)
        .def("__eq__", [](classic::SweepStates const& a, classic::SweepStates const& b) { return a == b; }, nb::is_operator())
        .def("__repr__", [](classic::SweepStates const& s) {
            return "SweepStates(" + std::to_string(s.right_to_left) + ", " + std::to_string(s.left_to_right) + ")";
        });
    m.def("fixed_sweeps", &classic::fixed_sweeps, nb::arg("sweeps"), nb::arg("states"), "`sweeps` sweeps keeping `states` states.");
    m.def("ramped_sweeps", &classic::ramped_sweeps, nb::arg("sweeps"), nb::arg("maximum"),
          "The historical MAXSTATES ramp: half sweep k keeps floor((k + 1) * maximum / (2 * sweeps)) states.");

    nb::class_<classic::Options>(m, "Options",
        "DMRG controls: the sweep schedule, warm-up states, targeted states, optional truncation-error "
        "control, Lanczos tolerance, arithmetic, the conserved sector and the scratch directory.")
        .def(nb::init<>())
        .def_prop_rw("sweeps", [](classic::Options const& o) { return o.sweeps; },
                     [](classic::Options& o, nb::handle value) { o.sweeps = sweeps_of(value); },
                     "A list of SweepStates or (right_to_left, left_to_right) pairs; see fixed_sweeps and ramped_sweeps.")
        .def_rw("warmup_states", &classic::Options::warmup_states)
        .def_rw("warmup_configuration", &classic::Options::warmup_configuration,
                "One local-state index per site: the warm-up targets this product configuration's prefix charges "
                "instead of uniform shares of the sector (for unevenly distributed ground states). [] keeps "
                "the uniform shares.")
        .def_rw("num_eigenvalues", &classic::Options::num_eigenvalues, "Lowest states targeted together.")
        .def_rw("truncation_error", &classic::Options::truncation_error)
        .def_rw("truncation_max_states", &classic::Options::truncation_max_states)
        .def_rw("lanczos_tolerance", &classic::Options::lanczos_tolerance)
        .def_prop_rw("arithmetic", [](classic::Options const& o) { return o.arithmetic; },
                     [](classic::Options& o, nb::handle value) { o.arithmetic = arithmetic_of(value); },
                     "'real' (real coefficients and factors) or 'complex'.")
        .def_rw("sector", &classic::Options::sector, "Fixed totals of declared charges (name -> value).")
        .def_prop_rw("scratch_parent", [](classic::Options const& o) { return o.scratch_parent.string(); },
                     [](classic::Options& o, std::string const& path) { o.scratch_parent = path; },
                     "Parent of the private scratch directory ('' selects the system temporary directory).")
        .def_rw("verbose", &classic::Options::verbose)
        .def_rw("keep_mps", &classic::Options::keep_mps,
                "Also return every targeted state as a MatrixProductState (Result.mps).")
        .def("__repr__", [](classic::Options const& o) {
            return "Options(" + std::to_string(o.sweeps.size()) + " sweeps, warmup_states=" + std::to_string(o.warmup_states)
                 + ", num_eigenvalues=" + std::to_string(o.num_eigenvalues)
                 + (o.arithmetic == classic::Arithmetic::complex ? ", complex)" : ", real)");
        });

    nb::class_<classic::Measurements>(m, "Measurements",
        "Operators measured on every targeted state: `expectation` (name -> operator, one value <psi|O|psi>) "
        "and `terms` (name -> operator, one value per term: a column per site, per pair, ...).")
        .def(nb::init<>())
        .def("__init__", [](classic::Measurements* self, std::map<std::string, mdl::CompiledOperator> expectation,
                            std::map<std::string, mdl::CompiledOperator> terms) {
            new (self) classic::Measurements{std::move(expectation), std::move(terms)};
        }, nb::arg("expectation") = std::map<std::string, mdl::CompiledOperator>{},
           nb::arg("terms") = std::map<std::string, mdl::CompiledOperator>{})
        .def_rw("expectation", &classic::Measurements::expectation)
        .def_rw("terms", &classic::Measurements::terms);

    nb::enum_<classic::MeasurementKind>(m, "MeasurementKind")
        .value("expectation", classic::MeasurementKind::expectation)
        .value("terms", classic::MeasurementKind::terms);
    nb::class_<classic::Measurement>(m, "Measurement", "Measured values: rows follow the targeted states, columns the terms.")
        .def_ro("kind", &classic::Measurement::kind)
        .def_ro("sites", &classic::Measurement::sites, "The factor sites of every column (empty for an expectation value).")
        .def_prop_ro("values", [](classic::Measurement const& measurement) {
            std::size_t const rows = measurement.values.size();
            std::size_t const columns = rows ? measurement.values.front().size() : measurement.sites.size();
            std::vector<double> values;
            values.reserve(rows * columns);
            for (auto const& row : measurement.values) {
                if (row.size() != columns) throw std::logic_error("DMRG measurement rows differ in length");
                values.insert(values.end(), row.begin(), row.end());
            }
            return matrix_to_numpy(std::move(values), rows, columns);
        }, "(targets, columns) array of real parts.")
        .def("__repr__", [](classic::Measurement const& measurement) {
            return "Measurement(" + std::string(measurement.kind == classic::MeasurementKind::terms ? "terms" : "expectation")
                 + ", " + std::to_string(measurement.values.size()) + "x" + std::to_string(measurement.sites.size()) + ")";
        });

    nb::class_<classic::Iteration>(m, "Iteration", "One completed DMRG iteration.")
        .def_prop_ro("direction", [](classic::Iteration const& it) { return direction_name(it.direction); })
        .def_ro("index", &classic::Iteration::index)
        .def_ro("energy", &classic::Iteration::energy)
        .def_ro("truncation_error", &classic::Iteration::truncation_error)
        .def_ro("entropy", &classic::Iteration::entropy)
        .def("__repr__", [](classic::Iteration const& it) {
            return "Iteration(" + std::string(direction_name(it.direction)) + ", " + std::to_string(it.index) + ", energy="
                 + std::to_string(it.energy) + ", truncation_error=" + std::to_string(it.truncation_error) + ")";
        });

    nb::class_<classic::ResolvedOptions>(m, "ResolvedOptions", "The scientific controls that took effect.")
        .def_ro("sweeps", &classic::ResolvedOptions::sweeps)
        .def_ro("warmup_states", &classic::ResolvedOptions::warmup_states)
        .def_ro("warmup_configuration", &classic::ResolvedOptions::warmup_configuration)
        .def_ro("num_eigenvalues", &classic::ResolvedOptions::num_eigenvalues)
        .def_ro("truncation_error", &classic::ResolvedOptions::truncation_error)
        .def_ro("truncation_max_states", &classic::ResolvedOptions::truncation_max_states)
        .def_ro("lanczos_tolerance", &classic::ResolvedOptions::lanczos_tolerance)
        .def_ro("arithmetic", &classic::ResolvedOptions::arithmetic)
        .def_ro("sector", &classic::ResolvedOptions::sector);

    nb::class_<classic::MatrixProductState>(m, "MatrixProductState",
        "A targeted state after the final sweep as the engine holds it: site tensors indexed "
        "[left bond, local state, right bond], local states in each site's PhysicalBasis order, mixed "
        "canonical with the norm at `center` (left-orthonormal before it, right-orthonormal after it). "
        "With fermions the product states are creators in site order applied to the vacuum, as in pyalps.ed.")
        .def_prop_ro("tensors", [](classic::MatrixProductState const& mps) {
            nb::list out;
            for (auto const& site : mps.sites) {
                std::vector<std::size_t> shape{site.left, site.local, site.right};
                if (site.complex.empty())
                    out.append(alps::python::make_numpy_array<double>(std::vector<double>(site.real), shape));
                else
                    out.append(alps::python::make_numpy_array<Complex>(std::vector<Complex>(site.complex), shape));
            }
            return out;
        }, "The site tensors as NumPy arrays (float64 in real arithmetic, complex128 in complex), copied on access.")
        .def_ro("center", &classic::MatrixProductState::center)
        .def_prop_ro("bond_dimensions", [](classic::MatrixProductState const& mps) {
            std::vector<std::size_t> out;
            for (std::size_t i = 0; i + 1 < mps.sites.size(); ++i) out.push_back(mps.sites[i].right);
            return out;
        }, "The dimension of each inner bond, between sites i and i + 1.")
        .def_prop_ro("local_dimensions", [](classic::MatrixProductState const& mps) {
            std::vector<std::size_t> out;
            for (auto const& site : mps.sites) out.push_back(site.local);
            return out;
        })
        .def_prop_ro("is_complex", [](classic::MatrixProductState const& mps) {
            return !mps.sites.empty() && !mps.sites.front().complex.empty();
        })
        .def("__len__", [](classic::MatrixProductState const& mps) { return mps.sites.size(); })
        .def("__repr__", [](classic::MatrixProductState const& mps) {
            std::size_t largest = 1;
            for (auto const& site : mps.sites) largest = std::max(largest, site.right);
            return "MatrixProductState(" + std::to_string(mps.sites.size()) + " sites, max bond "
                 + std::to_string(largest) + ", center " + std::to_string(mps.center) + ")";
        });

    nb::class_<classic::Result>(m, "Result", "A finished solve: energies, measurements, the iterations observed and the options.")
        .def_prop_ro("energies", [](classic::Result const& r) { return vector_to_numpy(std::vector<double>(r.energies)); },
                     "Ascending energies of the targeted states.")
        .def_ro("measurements", &classic::Result::measurements)
        .def_ro("iterations", &classic::Result::iterations)
        .def_ro("final_truncation_error", &classic::Result::final_truncation_error)
        .def_ro("options", &classic::Result::options)
        .def_ro("mps", &classic::Result::mps,
                "With keep_mps: one MatrixProductState per targeted state, in energy order; otherwise empty.")
        .def_prop_ro("energy", [](classic::Result const& r) { return r.energies.empty() ? 0. : r.energies.front(); },
                     "The ground-state energy.")
        .def_prop_ro("iteration_energies", [](classic::Result const& r) {
            std::vector<double> values;
            for (auto const& it : r.iterations) values.push_back(it.energy);
            return vector_to_numpy(std::move(values));
        })
        .def_prop_ro("iteration_truncation_errors", [](classic::Result const& r) {
            std::vector<double> values;
            for (auto const& it : r.iterations) values.push_back(it.truncation_error);
            return vector_to_numpy(std::move(values));
        })
        .def_prop_ro("iteration_entropies", [](classic::Result const& r) {
            std::vector<double> values;
            for (auto const& it : r.iterations) values.push_back(it.entropy);
            return vector_to_numpy(std::move(values));
        })
        .def("__repr__", [](classic::Result const& r) {
            return "Result(energy=" + (r.energies.empty() ? std::string("?") : std::to_string(r.energies.front())) + ", "
                 + std::to_string(r.energies.size()) + " states, " + std::to_string(r.iterations.size()) + " iterations)";
        });

    m.def("solve", [](mdl::Hamiltonian const& model, classic::Options const& options, classic::Measurements const& measurements,
                      nb::object callback, std::optional<std::string> const& log_file) {
        return run(nullptr, nullptr, &model, options, measurements, std::move(callback), log_file);
    }, nb::arg("model"), nb::arg("options") = classic::Options{}, nb::arg("measurements") = classic::Measurements{},
       nb::arg("callback") = nb::none(), nb::arg("log_file") = nb::none(),
       "Validate Hermiticity, compile the model under its declared charges and solve it on its physical basis. "
       "`callback(iteration)` observes every completed iteration and may raise to abort.");
    m.def("solve", [](mdl::CompiledOperator const& hamiltonian, mdl::PhysicalBasis const& basis, classic::Options const& options,
                      classic::Measurements const& measurements, nb::object callback, std::optional<std::string> const& log_file) {
        return run(&hamiltonian, &basis, nullptr, options, measurements, std::move(callback), log_file);
    }, nb::arg("hamiltonian"), nb::arg("basis"), nb::arg("options") = classic::Options{},
       nb::arg("measurements") = classic::Measurements{}, nb::arg("callback") = nb::none(), nb::arg("log_file") = nb::none(),
       "Solve a compiled Hermitian operator on its physical basis.");
    nb::class_<classic::MolecularResult>(m, "MolecularResult",
        "A finished molecular solve: the energy (with the core energy), the final-state RDMs in the PySCF "
        "convention and the iterations observed.")
        .def_ro("orbitals", &classic::MolecularResult::orbitals)
        .def_ro("energy", &classic::MolecularResult::energy)
        .def_prop_ro("one_rdm_alpha", [](classic::MolecularResult const& r) {
            return matrix_to_numpy(std::vector<double>(r.one_rdm_alpha), r.orbitals, r.orbitals);
        }, "dm1[p, q] = <c+[q, alpha] c[p, alpha]>.")
        .def_prop_ro("one_rdm_beta", [](classic::MolecularResult const& r) {
            return matrix_to_numpy(std::vector<double>(r.one_rdm_beta), r.orbitals, r.orbitals);
        }, "dm1[p, q] = <c+[q, beta] c[p, beta]>.")
        .def_prop_ro("one_rdm", [](classic::MolecularResult const& r) {
            std::vector<double> sum(r.one_rdm_alpha);
            for (std::size_t i = 0; i < sum.size(); ++i) sum[i] += r.one_rdm_beta[i];
            return matrix_to_numpy(std::move(sum), r.orbitals, r.orbitals);
        }, "The spin-summed one-particle density matrix.")
        .def_prop_ro("two_rdm", [](classic::MolecularResult const& r) {
            auto const n = r.orbitals;
            return alps::python::make_numpy_array<double>(std::vector<double>(r.two_rdm), {n, n, n, n});
        }, "dm2[p, q, r, s] = sum_ab <c+[p, a] c+[r, b] c[s, b] c[q, a]> (PySCF make_rdm12).")
        .def_prop_ro("two_rdm_aa", [](classic::MolecularResult const& r) {
            auto const n = r.orbitals;
            return alps::python::make_numpy_array<double>(std::vector<double>(r.two_rdm_aa), {n, n, n, n});
        }, "dm2aa[p, q, r, s] = <a+[p] a+[r] a[s] a[q]> (PySCF make_rdm12s).")
        .def_prop_ro("two_rdm_ab", [](classic::MolecularResult const& r) {
            auto const n = r.orbitals;
            return alps::python::make_numpy_array<double>(std::vector<double>(r.two_rdm_ab), {n, n, n, n});
        }, "dm2ab[p, q, r, s] = <a+[p] b+[r] b[s] a[q]> (PySCF make_rdm12s).")
        .def_prop_ro("two_rdm_bb", [](classic::MolecularResult const& r) {
            auto const n = r.orbitals;
            return alps::python::make_numpy_array<double>(std::vector<double>(r.two_rdm_bb), {n, n, n, n});
        }, "dm2bb[p, q, r, s] = <b+[p] b+[r] b[s] b[q]> (PySCF make_rdm12s).")
        .def_ro("observables", &classic::MolecularResult::observables,
                "The requested observables' final-state expectation values, by name.")
        .def_prop_ro("ci_alpha", [](classic::MolecularResult const& r) {
            return vector_to_numpy(std::vector<std::uint64_t>(r.ci_alpha));
        }, "Alpha occupation word of each CI entry (bit p = orbital p); empty unless ci_vector was requested.")
        .def_prop_ro("ci_beta", [](classic::MolecularResult const& r) {
            return vector_to_numpy(std::vector<std::uint64_t>(r.ci_beta));
        }, "Beta occupation word of each CI entry.")
        .def_prop_ro("ci_coefficients", [](classic::MolecularResult const& r) {
            return vector_to_numpy(std::vector<double>(r.ci_coefficients));
        }, "CI amplitudes in the ALPS mode order (orbitals ascending, down before up), sorted by alpha then beta "
           "word, largest amplitude positive.")
        .def_ro("iterations", &classic::MolecularResult::iterations)
        .def_ro("final_truncation_error", &classic::MolecularResult::final_truncation_error)
        .def_ro("mps", &classic::MolecularResult::mps,
                "With keep_mps: the final state as one MatrixProductState over the orbitals, on spinful_fermion "
                "sites (states 0, up, down, up-down); otherwise empty.")
        .def_prop_ro("iteration_energies", [](classic::MolecularResult const& r) {
            std::vector<double> values;
            for (auto const& it : r.iterations) values.push_back(it.energy);
            return vector_to_numpy(std::move(values));
        }, "Electronic energies of the iterations (without the core energy).")
        .def("__repr__", [](classic::MolecularResult const& r) {
            return "MolecularResult(energy=" + std::to_string(r.energy) + ", " + std::to_string(r.orbitals)
                 + " orbitals, " + std::to_string(r.iterations.size()) + " iterations)";
        });
    m.def("solve_molecular", &run_molecular, nb::arg("one_body").noconvert(), nb::arg("two_body").noconvert(),
          nb::arg("n_alpha"), nb::arg("n_beta"), nb::arg("core_energy") = 0., nb::arg("integral_cutoff") = 0.,
          nb::arg("observables") = std::map<std::string, mdl::TermSum>{}, nb::arg("ci_vector") = false,
          nb::arg("options") = classic::Options{}, nb::arg("callback") = nb::none(), nb::arg("log_file") = nb::none(),
          "Classic block DMRG of a restricted real molecular Hamiltonian (chemist-notation integrals with the "
          "real-orbital symmetries) in the sector of n_alpha and n_beta electrons; returns a MolecularResult.");
    nb::class_<classic::Drive>(m, "Drive",
        "A time-dependent term of the final Hamiltonian: coefficients[n] * op during step n.")
        .def("__init__", [](classic::Drive* self, mdl::CompiledOperator op, std::vector<double> coefficients) {
            new (self) classic::Drive{std::move(op), std::move(coefficients)};
        }, nb::arg("op"), nb::arg("coefficients"))
        .def_rw("op", &classic::Drive::op)
        .def_rw("coefficients", &classic::Drive::coefficients, "One coefficient per step.");
    nb::class_<classic::Evolution>(m, "Evolution",
        "Controls of a time evolution: the step, the step count, real or imaginary time, the kept states or "
        "truncation error, the measured steps, the drives and whether to keep the states.")
        .def(nb::init<>())
        .def_rw("time_step", &classic::Evolution::time_step)
        .def_rw("steps", &classic::Evolution::steps)
        .def_rw("imaginary", &classic::Evolution::imaginary, "exp(-tau H), renormalized, instead of exp(-i t H).")
        .def_rw("states", &classic::Evolution::states, "The most states a truncation keeps.")
        .def_rw("truncation_error", &classic::Evolution::truncation_error,
                "Keep as many states as this discarded weight needs, up to `states`.")
        .def_rw("measure_steps", &classic::Evolution::measure_steps, "Measured steps (0 is the initial state); empty: all.")
        .def_rw("drives", &classic::Evolution::drives)
        .def_rw("keep_mps", &classic::Evolution::keep_mps, "Also return the state at every measured step.")
        .def("__repr__", [](classic::Evolution const& e) {
            return "Evolution(time_step=" + std::to_string(e.time_step) + ", steps=" + std::to_string(e.steps)
                 + (e.imaginary ? ", imaginary" : ", real") + " time, states=" + std::to_string(e.states) + ")";
        });
    nb::class_<classic::EvolutionResult>(m, "EvolutionResult",
        "A finished time evolution: the preparation's Result, then per measured step its time, the final "
        "Hamiltonian's energy, the measurements and the kept states; per step the discarded weight.")
        .def_ro("initial", &classic::EvolutionResult::initial,
                "The DMRG that prepared the initial state (energy, iterations, resolved options).")
        .def_prop_ro("steps", [](classic::EvolutionResult const& r) {
            return vector_to_numpy(std::vector<std::uint64_t>(r.steps.begin(), r.steps.end()));
        }, "The measured steps, ascending.")
        .def_prop_ro("times", [](classic::EvolutionResult const& r) { return vector_to_numpy(std::vector<double>(r.times)); },
                     "steps * time_step (tau in imaginary time).")
        .def_prop_ro("energies", [](classic::EvolutionResult const& r) {
            return vector_to_numpy(std::vector<double>(r.energies));
        }, "<psi|H_final|psi> at the measured steps, without the drives.")
        .def_ro("measurements", &classic::EvolutionResult::measurements, "One row per measured step.")
        .def_prop_ro("truncation_errors", [](classic::EvolutionResult const& r) {
            return vector_to_numpy(std::vector<double>(r.truncation_errors));
        }, "Discarded weight of each step, steps + 1 entries (entry 0: the initial state's measurement).")
        .def_ro("mps", &classic::EvolutionResult::mps, "With keep_mps: the state at every measured step.")
        .def("__repr__", [](classic::EvolutionResult const& r) {
            return "EvolutionResult(" + std::to_string(r.steps.size()) + " measured steps, "
                 + std::to_string(r.truncation_errors.empty() ? 0 : r.truncation_errors.size() - 1) + " steps)";
        });
    m.def("evolve", [](mdl::Hamiltonian const& initial, mdl::Hamiltonian const& final_model,
                       classic::Evolution const& evolution, classic::Options const& options,
                       classic::Measurements const& measurements, nb::object callback,
                       std::optional<std::string> const& log_file) {
        return run_evolution(nullptr, nullptr, nullptr, &initial, &final_model, evolution, options, measurements,
                             std::move(callback), log_file);
    }, nb::arg("initial"), nb::arg("final"), nb::arg("evolution"), nb::arg("options") = classic::Options{},
       nb::arg("measurements") = classic::Measurements{}, nb::arg("callback") = nb::none(),
       nb::arg("log_file") = nb::none(),
       "Evolve the DMRG ground state of `initial` (in options.sector) under `final`, two models on one physical "
       "basis; returns an EvolutionResult.");
    m.def("evolve", [](mdl::CompiledOperator const& initial, mdl::CompiledOperator const& final_op,
                       mdl::PhysicalBasis const& basis, classic::Evolution const& evolution,
                       classic::Options const& options, classic::Measurements const& measurements,
                       nb::object callback, std::optional<std::string> const& log_file) {
        return run_evolution(&initial, &final_op, &basis, nullptr, nullptr, evolution, options, measurements,
                             std::move(callback), log_file);
    }, nb::arg("initial"), nb::arg("final"), nb::arg("basis"), nb::arg("evolution"),
       nb::arg("options") = classic::Options{}, nb::arg("measurements") = classic::Measurements{},
       nb::arg("callback") = nb::none(), nb::arg("log_file") = nb::none(),
       "Evolve under a compiled final Hamiltonian, both operators on `basis`.");
    m.def("site_terms", &classic::site_terms, nb::arg("model"), nb::arg("op"),
          "One term per site carrying the operator `op`, in site order (one measurement column per site).");
    m.def("pair_terms", &classic::pair_terms, nb::arg("model"), nb::arg("first"), nb::arg("second"),
          "One term per ordered site pair (i, j), i-major, with the on-site product at i == j.");
}
