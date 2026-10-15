// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Adapted from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, files
// src/solvers/tensor_network/src/classic_dmrg/dmrg.hpp (the dmtk driver of
// the ALPS dmrg application, DMRG<T>::run, create_site_operator,
// build_*_operator and save_results) and classic.ipp (option validation,
// sector feasibility and solve). Adaptation: the engine is fed from an
// ALPS::models CompiledOperator and PhysicalBasis instead of graph_helper,
// model_helper and alps::SiteOperator; local factors are split into
// definite-charge components and like terms merged before dmtk sees them;
// measurements are compiled operators; no HDF5, parameters or catalog
// snapshot. The dmtk calls themselves follow the alpsx driver line by line.
#include <alps/dmrg/classic.hpp>
#include "invocation.hpp"

#include <alps/models/terms.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace alps::dmrg::classic {
namespace {

using Complex = std::complex<double>;
using Key = std::vector<std::pair<double, double>>; // exact matrix entries, orderable

Key key_of(std::vector<Complex> const& matrix) {
    Key key;
    key.reserve(matrix.size());
    for (auto const& v : matrix) key.emplace_back(v.real(), v.imag());
    return key;
}

// ------------------------------------------------------------------ options

void bounded_count(std::size_t value, char const* name) {
    if (!value || value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument(std::string("classic DMRG ") + name + " must be a positive 32-bit integer");
}

// Twice a half-integer charge value, or an error naming what was not one.
long twice(double value, std::string const& what) {
    double const doubled = 2. * value;
    if (!std::isfinite(doubled) || std::abs(doubled - std::round(doubled)) > 1e-9
        || doubled <= std::numeric_limits<short>::min() || doubled >= std::numeric_limits<short>::max())
        throw std::invalid_argument("classic DMRG " + what + " must be a representable half-integer");
    return static_cast<long>(std::llround(doubled));
}

void validate_options(Options const& options) {
    bounded_count(options.sweeps.size(), "sweep count");
    if (options.sweeps.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) / 2)
        throw std::invalid_argument("classic DMRG half-sweep count exceeds the integer limit");
    for (auto const& sweep : options.sweeps) {
        bounded_count(sweep.right_to_left, "right_to_left states");
        bounded_count(sweep.left_to_right, "left_to_right states");
    }
    bounded_count(options.warmup_states, "warmup_states");
    bounded_count(options.num_eigenvalues, "num_eigenvalues");
    if (options.arithmetic != Arithmetic::real && options.arithmetic != Arithmetic::complex)
        throw std::invalid_argument("unknown classic DMRG arithmetic");
    if (options.truncation_error && (!std::isfinite(*options.truncation_error)
                                     || *options.truncation_error <= 0. || *options.truncation_error >= 1.))
        throw std::invalid_argument("classic DMRG truncation_error must be finite and between zero and one");
    if (options.lanczos_tolerance && (!std::isfinite(*options.lanczos_tolerance) || *options.lanczos_tolerance <= 0.))
        throw std::invalid_argument("classic DMRG lanczos_tolerance must be finite and positive");
    if (options.truncation_max_states) {
        bounded_count(*options.truncation_max_states, "truncation_max_states");
        if (!options.truncation_error)
            throw std::invalid_argument("classic DMRG truncation_max_states requires truncation_error");
    }
    for (auto const& [name, value] : options.sector) {
        if (name.empty() || name.find_first_of(" ,\t\r\n") != std::string::npos)
            throw std::invalid_argument("classic DMRG sector charge names must be nonempty without separators");
        twice(value, "sector target " + name);
    }
    if (options.verbose < 0) throw std::invalid_argument("classic DMRG verbose must not be negative");
}

// ------------------------------------------------------------------ problem

struct Charge {
    std::string name;
    bool fermionic = false;
    bool conserved = false;
    long target = 0; // twice the target value of a conserved charge
};

// One distinct local space; lattice sites with equal spaces share it.
struct Site {
    int dim = 0;
    std::vector<double> parity;           // +1/-1 per state
    std::vector<std::vector<long>> twice; // [charge][state]
};

// One definite-charge, definite-parity component of a local factor.
struct Factor {
    int site = 0;
    std::vector<Complex> matrix; // row-major, row = output state
    std::vector<long> delta;     // [charge]: twice (output - input)
    bool odd = false;
};

struct ReadyTerm {
    Complex coefficient;
    std::vector<Factor> factors;
};

struct Column {
    std::vector<int> sites;
    std::vector<ReadyTerm> components;
    Complex scalar;
};

struct Plan {
    std::string name;
    MeasurementKind kind = MeasurementKind::expectation;
    std::vector<Column> columns;
};

struct Problem {
    int num_sites = 0;
    std::vector<Site> sites;
    std::vector<int> site_index; // lattice site -> sites
    std::vector<Charge> charges;
    std::vector<ReadyTerm> hamiltonian;
    Complex scalar; // identity part of the Hamiltonian
    std::vector<Plan> measurements;
};

bool is_odd_integer(double value) { return std::llround(value) % 2 != 0; }

std::vector<Charge> select_charges(models::PhysicalBasis const& basis, Options const& options) {
    std::vector<Charge> charges;
    auto const& declared = basis.conserved_quantum_numbers;
    for (auto const& [name, value] : options.sector) {
        if (std::find(declared.begin(), declared.end(), name) == declared.end())
            throw std::invalid_argument("classic DMRG sector names a charge the physical basis does not declare conserved: " + name);
        charges.push_back(Charge{name, false, true, twice(value, "sector target " + name)});
    }
    bool const strict = basis.fermion_convention == models::FermionConvention::strict;
    bool fermionic = false;
    if (strict)
        for (auto const& site : basis.sites)
            for (double p : site.parity) fermionic = fermionic || p < 0.;
    if (!fermionic) return charges;
    // The engine derives fermion signs from integer quantum numbers flagged
    // fermionic: find one whose parity is the state parity everywhere.
    std::vector<std::string> candidates{"N"};
    for (auto const& site : basis.sites)
        for (auto const& [name, values] : site.quantum_numbers)
            if (std::find(candidates.begin(), candidates.end(), name) == candidates.end()) candidates.push_back(name);
    for (auto const& name : candidates) {
        bool valid = false;
        bool ok = true;
        for (auto const& site : basis.sites) {
            auto found = site.quantum_numbers.find(name);
            for (std::size_t state = 0; state < site.parity.size() && ok; ++state) {
                double const q = found == site.quantum_numbers.end() ? 0. : found->second.at(state);
                if (!std::isfinite(q) || std::abs(q - std::round(q)) > 1e-9) ok = false;
                else if (is_odd_integer(q) != (site.parity[state] < 0.)) ok = false;
                else valid = true;
            }
            if (!ok) break;
        }
        if (!ok || !valid) continue;
        auto existing = std::find_if(charges.begin(), charges.end(), [&](Charge const& c) { return c.name == name; });
        if (existing != charges.end()) existing->fermionic = true;
        else charges.push_back(Charge{name, true, false, 0});
        return charges;
    }
    throw std::invalid_argument("classic DMRG needs an integer charge whose parity equals the fermion parity "
                                "on every site (for example the particle number N)");
}

// Components of one local factor with definite charge deltas and parity change.
std::vector<Factor> split(models::LocalFactor const& factor, Site const& site, std::size_t num_charges) {
    int const dim = site.dim;
    if (factor.matrix.size() != static_cast<std::size_t>(dim) * static_cast<std::size_t>(dim))
        throw std::invalid_argument("classic DMRG local factor matrix does not match its site dimension");
    std::map<std::pair<std::vector<long>, int>, std::vector<Complex>> groups;
    for (int out = 0; out < dim; ++out)
        for (int in = 0; in < dim; ++in) {
            Complex const v = factor.matrix[static_cast<std::size_t>(out) * dim + in];
            if (!std::isfinite(v.real()) || !std::isfinite(v.imag()))
                throw std::invalid_argument("classic DMRG local operator contains a nonfinite matrix element");
            if (v == Complex{}) continue;
            std::vector<long> delta(num_charges);
            for (std::size_t c = 0; c < num_charges; ++c) delta[c] = site.twice[c][out] - site.twice[c][in];
            int const change = site.parity[out] * site.parity[in] < 0. ? -1 : 1;
            auto& matrix = groups[{delta, change}];
            if (matrix.empty()) matrix.assign(factor.matrix.size(), Complex{});
            matrix[static_cast<std::size_t>(out) * dim + in] = v;
        }
    std::vector<Factor> components;
    for (auto& [key, matrix] : groups)
        components.push_back(Factor{factor.site, std::move(matrix), key.first, key.second < 0});
    return components;
}

// Expand one compiled term into definite-component terms. Scalar terms go to
// `scalar`; vanishing terms contribute nothing.
void expand(models::LocalTerm const& term, Problem const& problem, std::vector<ReadyTerm>& out, Complex& scalar,
            char const* what) {
    if (!std::isfinite(term.coefficient.real()) || !std::isfinite(term.coefficient.imag()))
        throw std::invalid_argument(std::string("classic DMRG ") + what + " coefficient is nonfinite");
    if (term.factors.empty()) {
        scalar += term.coefficient;
        return;
    }
    if (term.factors.size() > 2)
        throw std::invalid_argument(std::string("classic DMRG ") + what
                                    + " terms may act on at most two sites; got " + std::to_string(term.factors.size()));
    for (std::size_t k = 0; k < term.factors.size(); ++k) {
        int const site = term.factors[k].site;
        if (site < 0 || site >= problem.num_sites)
            throw std::invalid_argument(std::string("classic DMRG ") + what + " factor site is out of range");
        if (k && site <= term.factors[k - 1].site)
            throw std::invalid_argument(std::string("classic DMRG ") + what + " factors must be on distinct ascending sites");
    }
    std::vector<std::vector<Factor>> components;
    for (auto const& factor : term.factors) {
        components.push_back(split(factor, problem.sites[problem.site_index[factor.site]], problem.charges.size()));
        if (components.back().empty()) return; // an exactly zero factor
    }
    auto emit = [&](std::vector<Factor> factors) {
        std::size_t odd = 0;
        for (auto const& f : factors) odd += f.odd ? 1 : 0;
        if (odd % 2)
            throw std::invalid_argument(std::string("classic DMRG ") + what
                                        + " contains an odd-parity term, which the engine cannot represent");
        std::vector<int> expected;
        if (odd == 2)
            for (int s = factors[0].site; s < factors[1].site; ++s) expected.push_back(s);
        if (term.jordan_wigner != expected)
            throw std::invalid_argument(std::string("classic DMRG ") + what
                                        + " term carries a Jordan-Wigner string other than the strict convention's [i, j) "
                                          "between two odd factors");
        out.push_back(ReadyTerm{term.coefficient, std::move(factors)});
    };
    for (auto const& first : components[0]) {
        if (components.size() == 1) emit({first});
        else
            for (auto const& second : components[1]) emit({first, second});
    }
}

using TermKey = std::vector<std::pair<int, Key>>;

TermKey term_key(ReadyTerm const& term) {
    TermKey key;
    for (auto const& f : term.factors) key.emplace_back(f.site, key_of(f.matrix));
    return key;
}

// Sum like components, drop the vanished ones, then require that every
// surviving term conserves the sector's charges and that the sum is Hermitian.
std::vector<ReadyTerm> merge_hamiltonian(std::vector<ReadyTerm> const& terms, Problem const& problem) {
    std::vector<ReadyTerm> merged;
    std::map<TermKey, std::size_t> index;
    for (auto const& term : terms) {
        auto [pos, inserted] = index.emplace(term_key(term), merged.size());
        if (inserted) merged.push_back(term);
        else merged[pos->second].coefficient += term.coefficient;
    }
    double scale = 0.;
    for (auto const& term : merged) scale = std::max(scale, std::abs(term.coefficient));
    std::vector<ReadyTerm> kept;
    for (auto& term : merged) {
        if (term.coefficient == Complex{}) continue;
        for (std::size_t c = 0; c < problem.charges.size(); ++c) {
            if (!problem.charges[c].conserved) continue;
            long delta = 0;
            for (auto const& f : term.factors) delta += f.delta[c];
            if (delta != 0 && std::abs(term.coefficient) > 1e-12 * scale)
                throw std::invalid_argument("classic DMRG Hamiltonian does not conserve the requested charge: "
                                            + problem.charges[c].name);
            if (delta != 0) { term.coefficient = Complex{}; break; }
        }
        if (term.coefficient != Complex{}) kept.push_back(std::move(term));
    }
    if (kept.empty()) throw std::invalid_argument("classic DMRG Hamiltonian has no term");
    // Hermiticity: the adjoint of c * F_i F_j is conj(c) * F_j^+ F_i^+, which in
    // ascending site order is -conj(c) * F_i^+ F_j^+ for two odd factors.
    index.clear();
    for (std::size_t i = 0; i < kept.size(); ++i) index.emplace(term_key(kept[i]), i);
    for (auto const& term : kept) {
        ReadyTerm adjoint = term;
        adjoint.coefficient = std::conj(term.coefficient);
        if (term.factors.size() == 2 && term.factors[0].odd && term.factors[1].odd) adjoint.coefficient = -adjoint.coefficient;
        for (auto& f : adjoint.factors) {
            int const dim = problem.sites[problem.site_index[f.site]].dim;
            std::vector<Complex> dagger(f.matrix.size());
            for (int out = 0; out < dim; ++out)
                for (int in = 0; in < dim; ++in)
                    dagger[static_cast<std::size_t>(out) * dim + in] = std::conj(f.matrix[static_cast<std::size_t>(in) * dim + out]);
            f.matrix = std::move(dagger);
        }
        auto found = index.find(term_key(adjoint));
        Complex const present = found == index.end() ? Complex{} : kept[found->second].coefficient;
        if (std::abs(present - adjoint.coefficient) > 1e-10 * scale) {
            std::string where;
            for (auto const& f : term.factors) where += (where.empty() ? "" : ", ") + std::to_string(f.site);
            throw std::invalid_argument("classic DMRG Hamiltonian is not Hermitian (term on sites " + where + ")");
        }
    }
    return kept;
}

Problem prepare(models::CompiledOperator const& hamiltonian, models::PhysicalBasis const& basis,
                Options const& options, Measurements const& measurements) {
    models::validate_basis(basis);
    Problem problem;
    problem.num_sites = static_cast<int>(basis.sites.size());
    if (hamiltonian.site_dims.size() != basis.sites.size())
        throw std::invalid_argument("classic DMRG operator and physical basis disagree on the number of sites");
    if (problem.num_sites < 4) throw std::invalid_argument("classic DMRG requires at least four sites");
    problem.charges = select_charges(basis, options);
    if (problem.charges.size() > QN_MAX_SIZE)
        throw std::invalid_argument("classic DMRG tracks at most " + std::to_string(QN_MAX_SIZE)
                                    + " quantum numbers (the sector's charges plus the fermionic charge)");
    bool const strict = basis.fermion_convention == models::FermionConvention::strict;
    std::map<std::tuple<int, std::vector<double>, std::vector<std::vector<long>>>, int> distinct;
    for (std::size_t i = 0; i < basis.sites.size(); ++i) {
        auto const& local = basis.sites[i];
        int const dim = static_cast<int>(local.labels.size());
        if (hamiltonian.site_dims[i] != dim)
            throw std::invalid_argument("classic DMRG operator site dimension disagrees with the physical basis");
        if (dim < 1 || dim >= std::numeric_limits<short>::max())
            throw std::invalid_argument("classic DMRG needs a finite nonempty local basis on every site");
        Site site;
        site.dim = dim;
        site.parity = strict ? local.parity : std::vector<double>(local.parity.size(), 1.);
        for (auto const& charge : problem.charges) {
            auto found = local.quantum_numbers.find(charge.name);
            std::vector<long> row(static_cast<std::size_t>(dim), 0);
            if (found != local.quantum_numbers.end())
                for (int s = 0; s < dim; ++s) row[static_cast<std::size_t>(s)] = twice(found->second.at(static_cast<std::size_t>(s)), "charge " + charge.name);
            site.twice.push_back(std::move(row));
        }
        auto [pos, inserted] = distinct.emplace(std::make_tuple(dim, site.parity, site.twice), static_cast<int>(problem.sites.size()));
        if (inserted) problem.sites.push_back(std::move(site));
        problem.site_index.push_back(pos->second);
    }
    // Accumulated charges must stay within the engine's half-integer range.
    for (std::size_t c = 0; c < problem.charges.size(); ++c) {
        long low = 0, high = 0;
        for (int i = 0; i < problem.num_sites; ++i) {
            auto const& row = problem.sites[problem.site_index[i]].twice[c];
            low += *std::min_element(row.begin(), row.end());
            high += *std::max_element(row.begin(), row.end());
        }
        if (low <= std::numeric_limits<short>::min() || high >= std::numeric_limits<short>::max())
            throw std::invalid_argument("classic DMRG accumulated quantum numbers exceed the supported half-integer range");
        if (problem.charges[c].conserved && (problem.charges[c].target < low || problem.charges[c].target > high))
            throw std::invalid_argument("classic DMRG sector target is outside the local-basis charge bounds: "
                                        + problem.charges[c].name);
    }
    std::vector<ReadyTerm> expanded;
    for (auto const& term : hamiltonian.terms) expand(term, problem, expanded, problem.scalar, "Hamiltonian");
    problem.hamiltonian = merge_hamiltonian(expanded, problem);
    auto plan = [&](std::map<std::string, models::CompiledOperator> const& entries, MeasurementKind kind) {
        for (auto const& [name, op] : entries) {
            if (name.empty()) throw std::invalid_argument("classic DMRG measurement names must be nonempty");
            if (op.site_dims != hamiltonian.site_dims)
                throw std::invalid_argument("classic DMRG measurement operator is not on the Hamiltonian's sites: " + name);
            Plan entry{name, kind, {}};
            std::string const what = "measurement " + name;
            auto column = [&](std::vector<models::LocalTerm const*> const& terms, std::vector<int> sites) {
                Column col{std::move(sites), {}, Complex{}};
                std::vector<ReadyTerm> components;
                for (auto const* term : terms) expand(*term, problem, components, col.scalar, what.c_str());
                for (auto& component : components) {
                    bool leaves_sector = false;
                    for (std::size_t c = 0; c < problem.charges.size() && !leaves_sector; ++c) {
                        if (!problem.charges[c].conserved) continue;
                        long delta = 0;
                        for (auto const& f : component.factors) delta += f.delta[c];
                        leaves_sector = delta != 0; // zero expectation value in a definite sector
                    }
                    if (!leaves_sector) col.components.push_back(std::move(component));
                }
                entry.columns.push_back(std::move(col));
            };
            if (kind == MeasurementKind::expectation) {
                std::vector<models::LocalTerm const*> all;
                for (auto const& term : op.terms) all.push_back(&term);
                column(all, {});
            } else {
                for (auto const& term : op.terms) {
                    std::vector<int> sites;
                    for (auto const& f : term.factors) sites.push_back(f.site);
                    column({&term}, std::move(sites));
                }
            }
            problem.measurements.push_back(std::move(entry));
        }
    };
    plan(measurements.expectation, MeasurementKind::expectation);
    plan(measurements.terms, MeasurementKind::terms);
    return problem;
}

// Count product states in the requested sector, stopping at the number of
// targets, without constructing the product basis (alpsx validate_owned_sector).
void check_sector(Problem const& problem, std::size_t targets) {
    using Tuple = std::array<long, QN_MAX_SIZE>;
    std::vector<std::size_t> conserved;
    Tuple target{};
    for (std::size_t c = 0; c < problem.charges.size(); ++c)
        if (problem.charges[c].conserved) {
            target[conserved.size()] = problem.charges[c].target;
            conserved.push_back(c);
        }
    if (conserved.empty()) {
        std::size_t bound = 1;
        for (int i = 0; i < problem.num_sites && bound < targets; ++i)
            bound *= static_cast<std::size_t>(problem.sites[problem.site_index[i]].dim);
        if (bound < targets) throw std::invalid_argument("classic DMRG target count exceeds the Hilbert-space dimension");
        return;
    }
    using Choices = std::vector<std::pair<Tuple, std::size_t>>;
    std::vector<Choices> choices;
    for (auto const& site : problem.sites) {
        std::map<Tuple, std::size_t> counts;
        for (int s = 0; s < site.dim; ++s) {
            Tuple q{};
            for (std::size_t j = 0; j < conserved.size(); ++j) q[j] = site.twice[conserved[j]][static_cast<std::size_t>(s)];
            auto& count = counts[q];
            if (count < targets) ++count;
        }
        choices.emplace_back(counts.begin(), counts.end());
    }
    std::size_t const n = static_cast<std::size_t>(problem.num_sites);
    std::vector<Tuple> suffix_min(n + 1), suffix_max(n + 1);
    for (std::size_t i = n; i-- > 0;) {
        suffix_min[i] = suffix_min[i + 1];
        suffix_max[i] = suffix_max[i + 1];
        auto const& local = choices[static_cast<std::size_t>(problem.site_index[i])];
        for (std::size_t j = 0; j < conserved.size(); ++j) {
            long lo = std::numeric_limits<long>::max(), hi = std::numeric_limits<long>::min();
            for (auto const& [q, count] : local) { lo = std::min(lo, q[j]); hi = std::max(hi, q[j]); }
            suffix_min[i][j] += lo;
            suffix_max[i][j] += hi;
        }
    }
    using State = std::pair<std::size_t, Tuple>;
    struct Frame { State key; std::size_t next = 0, count = 0; };
    std::vector<Frame> stack{Frame{State{0, target}}};
    std::map<State, std::size_t> memo;
    auto finish = [&](std::size_t count) { memo.emplace(stack.back().key, count); stack.pop_back(); };
    while (!stack.empty()) {
        auto& frame = stack.back();
        auto const& [i, remaining] = frame.key;
        bool feasible = true;
        for (std::size_t j = 0; j < conserved.size(); ++j)
            if (remaining[j] < suffix_min[i][j] || remaining[j] > suffix_max[i][j]) feasible = false;
        if (!feasible) { finish(0); continue; }
        if (i == n) { finish(1); continue; }
        auto const& local = choices[static_cast<std::size_t>(problem.site_index[i])];
        if (frame.count == targets || frame.next == local.size()) { finish(frame.count); continue; }
        auto const& [q, multiplicity] = local[frame.next];
        Tuple residual = remaining;
        for (std::size_t j = 0; j < conserved.size(); ++j) residual[j] -= q[j];
        State child{i + 1, residual};
        auto found = memo.find(child);
        if (found == memo.end()) { stack.push_back(Frame{child}); continue; }
        auto const room = targets - frame.count;
        frame.count += found->second > room / multiplicity ? room : found->second * multiplicity;
        ++frame.next;
    }
    auto const dimension = memo.at(State{0, target});
    if (!dimension) throw std::invalid_argument("classic DMRG requested quantum-number sector is empty");
    if (dimension < targets) throw std::invalid_argument("classic DMRG target count exceeds the selected-sector dimension");
}

// ------------------------------------------------------------------- engine

template<class T> struct Engine;

template<class T>
bool end_iteration(dmtk::System<T>& system, std::size_t signal, void* data) {
    if (signal == dmtk::SYSTEM_SIGNAL_END_ITER) static_cast<Engine<T>*>(data)->observe(system);
    return false;
}

template<class T>
struct Engine {
    Problem const& problem;
    Options const& options;
    IterationCallback callback;
    std::vector<Iteration> rows;
    dmtk::Lattice lattice;
    std::vector<dmtk::Block<T>> blocks;
    dmtk::Hami<T> hami;
    dmtk::Hami<T> measured;
    std::vector<std::size_t> qn_index; // per charge
    int conserved_mask = 0;
    std::map<std::pair<int, Key>, std::string> operator_names;
    std::vector<std::vector<std::vector<std::size_t>>> layout; // [plan][column] -> measured term indices
    bool two_site_measurements = false;

    Engine(Problem const& p, Options const& o, IterationCallback cb)
        : problem(p), options(o), callback(std::move(cb)), lattice(p.num_sites, dmtk::OBC), hami(lattice) {}

    void observe(dmtk::System<T> const& s) {
        Iteration row{s.get_dir() == dmtk::RIGHT2LEFT ? Direction::right_to_left : Direction::left_to_right,
                      static_cast<std::size_t>(s.get_iter()), s.energy[0], s.truncation_error(), s.entropy()};
        rows.push_back(row);
        if (callback) callback(row);
    }

    static T scalar(Complex v, char const* what) {
        if constexpr (std::is_same_v<T, double>) {
            if (std::abs(v.imag()) > 1e-12 * std::max(1., std::abs(v)))
                throw std::invalid_argument(std::string("classic DMRG real arithmetic needs a real ") + what
                                            + "; select complex arithmetic");
            return v.real();
        } else {
            return v;
        }
    }

    void register_charges() {
        dmtk::QN::init();
        int mask = 0;
        for (auto const& charge : problem.charges) {
            auto const index = dmtk::QN::add_qn_index(charge.name, charge.fermionic);
            qn_index.push_back(index);
            mask |= 1 << index;
            if (charge.conserved) conserved_mask |= 1 << index;
        }
        dmtk::QN::set_qn_mask(mask);
    }

    void build_blocks() {
        blocks.resize(problem.sites.size());
        for (std::size_t b = 0; b < problem.sites.size(); ++b) {
            auto const& site = problem.sites[b];
            dmtk::Basis basis(static_cast<std::size_t>(site.dim));
            for (int s = 0; s < site.dim; ++s) {
                dmtk::QN qn;
                for (std::size_t c = 0; c < problem.charges.size(); ++c)
                    qn[static_cast<int>(qn_index[c])].set_half(static_cast<short>(site.twice[c][static_cast<std::size_t>(s)]));
                basis[static_cast<std::size_t>(s)] = dmtk::State(static_cast<std::size_t>(s), qn);
            }
            basis.reorder();
            blocks[b].resize(basis);
            blocks[b].set_lattice(dmtk::Lattice(1, dmtk::OBC));
        }
        for (int i = 0; i < problem.num_sites; ++i)
            hami.sites(static_cast<std::size_t>(i)) = &blocks[static_cast<std::size_t>(problem.site_index[static_cast<std::size_t>(i)])];
    }

    // The block operator of one component, registered once per (block, matrix).
    dmtk::BasicOp<T> local_operator(Factor const& factor) {
        int const b = problem.site_index[static_cast<std::size_t>(factor.site)];
        auto& block = blocks[static_cast<std::size_t>(b)];
        int const dim = problem.sites[static_cast<std::size_t>(b)].dim;
        auto const key = std::make_pair(b, key_of(factor.matrix));
        auto found = operator_names.find(key);
        if (found != operator_names.end()) {
            dmtk::BasicOp<T>* op = block(found->second.c_str(), 0);
            if (!op) throw std::logic_error("classic DMRG lost a registered block operator");
            return *op;
        }
        std::string const name = "o" + std::to_string(operator_names.size());
        // The block basis is sorted by quantum numbers (Basis::reorder); each
        // sorted state keeps its kit index in i1, which is where the dense
        // kit matrix is addressed. dest is a dmtk matrix: [row = out][col = in].
        auto const& basis = block.basis();
        dmtk::Matrix<T> dest(static_cast<std::size_t>(dim), static_cast<std::size_t>(dim));
        dmtk::QN dqn;
        bool first = true;
        for (std::size_t out = 0; out < static_cast<std::size_t>(dim); ++out)
            for (std::size_t in = 0; in < static_cast<std::size_t>(dim); ++in) {
                Complex const v = factor.matrix[basis[out].i1 * static_cast<std::size_t>(dim) + basis[in].i1];
                dest[out][in] = scalar(v, "local operator");
                if (std::abs(v) > 1.e-10) {
                    dmtk::QN const change = basis[out].qn() - basis[in].qn();
                    if (!first && change != dqn)
                        throw std::logic_error("classic DMRG matrix elements are inconsistent with the change in quantum numbers");
                    dqn = change;
                    first = false;
                }
            }
        dmtk::BasicOp<T> op(name.c_str(), 0);
        op.dqn = dqn;
        // Jordan-Wigner phases are gated on the explicit fermion flag; derive
        // it from the parity change of the fermionic quantum number (alpsx
        // dmrg.hpp, upstream ALPSim/ALPS#135).
        op.set_fermion(dqn.fermion_sign() == -1);
        op.resize(basis);
        op = dest;
        block.push_back(op);
        operator_names.emplace(key, name);
        return op;
    }

    dmtk::Term<T> make_term(ReadyTerm const& term) {
        dmtk::BasicOp<T> first = local_operator(term.factors[0]);
        first.set_site(term.factors[0].site);
        dmtk::Term<T> product;
        if (term.factors.size() == 1) {
            product = dmtk::Term<T>(first);
        } else {
            dmtk::BasicOp<T> second = local_operator(term.factors[1]);
            second.set_site(term.factors[1].site);
            product = first * second;
        }
        product.coef() = scalar(term.coefficient, "term coefficient");
        return product;
    }

    void build_hamiltonian() {
        for (auto const& term : problem.hamiltonian) hami += make_term(term);
        hami.reorder_terms();
    }

    void build_measurements() {
        measured = dmtk::Hami<T>(hami);
        measured.clear();
        layout.resize(problem.measurements.size());
        for (std::size_t p = 0; p < problem.measurements.size(); ++p) {
            auto const& plan = problem.measurements[p];
            layout[p].resize(plan.columns.size());
            for (std::size_t c = 0; c < plan.columns.size(); ++c) {
                auto const& column = plan.columns[c];
                if (plan.kind == MeasurementKind::terms) {
                    for (auto const& component : column.components) {
                        dmtk::Term<T> term = make_term(component);
                        // On-site terms are measured as weighted terms, two-site
                        // products as composite operators (alpsx dmrg.hpp).
                        if (component.factors.size() == 1) {
                            measured += term;
                        } else {
                            measured += dmtk::BasicOp<T>(term);
                            two_site_measurements = true;
                        }
                        layout[p][c].push_back(measured.size() - 1);
                    }
                } else if (!column.components.empty()) {
                    dmtk::Hami<T> sum;
                    for (auto const& component : column.components) {
                        sum += make_term(component);
                        two_site_measurements = two_site_measurements || component.factors.size() == 2;
                    }
                    sum.set_name(("m" + std::to_string(p)).c_str());
                    measured += dmtk::BasicOp<T>(sum);
                    layout[p][c].push_back(measured.size() - 1);
                }
            }
        }
    }

    Result run() {
        register_charges();
        build_blocks();
        build_hamiltonian();
        build_measurements();

        dmtk::System<T> S(hami, lattice, "ALPS");
        // A two-site product inside one block is re-measured for a second
        // target only from a stored block operator; otherwise its value stays
        // the one taken when the pair was central, during the ground-state
        // sweep. Keep the lighter unstored path when it is correct.
        S.set_store_products(options.num_eigenvalues > 1 && two_site_measurements);
        S.signal_connect(&end_iteration<T>, dmtk::SYSTEM_SIGNAL_END_ITER, this);
        if (options.truncation_error)
            S.set_error(*options.truncation_error,
                        options.truncation_max_states ? static_cast<int>(*options.truncation_max_states) : -1);
        if (options.lanczos_tolerance) S.set_lanczos_tolerance(*options.lanczos_tolerance);
        S.set_calc_gap(static_cast<int>(options.num_eigenvalues) - 1);
        std::size_t const sweeps = options.sweeps.size();
        dmtk::Matrix<std::size_t> nstates(2, sweeps);
        dmtk::QN target;
        for (std::size_t c = 0; c < problem.charges.size(); ++c)
            if (problem.charges[c].conserved)
                target[static_cast<int>(qn_index[c])].set_half(static_cast<short>(problem.charges[c].target));
        S.qnt = target;
        S.set_qn_mask(conserved_mask);
        S.set_use_hc(false);
        S.set_grow_symmetric(false);
        S.set_verbose(options.verbose);
        for (std::size_t i = 0; i < sweeps; ++i) {
            nstates(0, i) = options.sweeps[i].right_to_left;
            nstates(1, i) = options.sweeps[i].left_to_right;
        }
        S.start(sweeps, nstates);
        S.run(static_cast<int>(options.warmup_states));
        if (S.store_products()) S.set_full_sweep(true);
        S.corr += measured;
        S.final_sweep(options.sweeps.back().left_to_right, dmtk::RIGHT2LEFT, 1, true);
        // The unstored path measures while sweeping and measure() refreshes
        // what the final configuration can reach; the stored path measures
        // nothing during the sweep, so this call takes every value.
        S.measure();

        Result result;
        for (auto const& plan : problem.measurements)
            result.measurements.emplace(plan.name, Measurement{plan.kind, {}, {}});
        auto read = [&]() {
            for (std::size_t p = 0; p < problem.measurements.size(); ++p) {
                auto const& plan = problem.measurements[p];
                auto& measurement = result.measurements.at(plan.name);
                std::vector<double> row;
                for (std::size_t c = 0; c < plan.columns.size(); ++c) {
                    double value = plan.columns[c].scalar.real();
                    for (auto index : layout[p][c]) value += std::real(S.corr[index].value());
                    row.push_back(value);
                }
                measurement.values.push_back(std::move(row));
            }
        };
        read();
        for (std::size_t i = 1; i < S._target.size(); ++i) {
            S.gs = S._target[i];
            S.measure();
            read();
        }
        for (auto const& plan : problem.measurements) {
            auto& measurement = result.measurements.at(plan.name);
            for (auto const& column : plan.columns) measurement.sites.push_back(column.sites);
        }
        for (std::size_t i = 0; i < S.energy.size(); ++i) result.energies.push_back(S.energy[i] + problem.scalar.real());
        result.iterations = std::move(rows);
        result.final_truncation_error = S.truncation_error();
        auto& resolved = result.options;
        resolved.sweeps = options.sweeps;
        resolved.warmup_states = options.warmup_states;
        resolved.num_eigenvalues = options.num_eigenvalues;
        resolved.truncation_error = options.truncation_error;
        resolved.truncation_max_states = options.truncation_max_states;
        resolved.lanczos_tolerance = S.lanczos_tolerance();
        resolved.arithmetic = options.arithmetic;
        resolved.sector = options.sector;
        return result;
    }
};

// Lower one located product and keep exactly one column for it.
void append_column(models::CompiledOperator& out, models::Hamiltonian const& model, models::TermSum const& product) {
    auto lowered = model.operator_decomposition(product).op;
    if (lowered.terms.empty()) out.terms.push_back(models::LocalTerm{Complex{}, {}, {}});
    else if (lowered.terms.size() == 1) out.terms.push_back(std::move(lowered.terms.front()));
    else throw std::logic_error("classic DMRG measurement builder lowered one product into several terms");
}

} // namespace

std::vector<SweepStates> fixed_sweeps(std::size_t sweeps, std::size_t states) {
    bounded_count(sweeps, "sweep count");
    bounded_count(states, "state count");
    if (sweeps > static_cast<std::size_t>(std::numeric_limits<int>::max()) / 2)
        throw std::invalid_argument("classic DMRG half-sweep count exceeds the integer limit");
    return std::vector<SweepStates>(sweeps, SweepStates{states, states});
}

std::vector<SweepStates> ramped_sweeps(std::size_t sweeps, std::size_t maximum) {
    auto result = fixed_sweeps(sweeps, maximum);
    if (maximum < 2 * sweeps) throw std::invalid_argument("classic DMRG ramp would retain zero states");
    for (std::size_t i = 0; i < sweeps; ++i)
        result[i] = {static_cast<std::size_t>((std::uint64_t(2 * i + 1) * maximum) / (2 * sweeps)),
                     static_cast<std::size_t>((std::uint64_t(2 * i + 2) * maximum) / (2 * sweeps))};
    return result;
}

Result solve(models::CompiledOperator const& hamiltonian, models::PhysicalBasis const& basis,
             Options const& options, Measurements const& measurements, IterationCallback callback) {
    validate_options(options);
    Problem const problem = prepare(hamiltonian, basis, options, measurements);
    check_sector(problem, options.num_eigenvalues);
    dmtk::InvocationGuard invocation(options.scratch_parent, options.log);
    // Every engine object is destroyed before the guard restores the engine state.
    if (options.arithmetic == Arithmetic::complex) {
        Engine<Complex> engine(problem, options, std::move(callback));
        return engine.run();
    }
    Engine<double> engine(problem, options, std::move(callback));
    return engine.run();
}

Result solve(models::Hamiltonian const& model, Options const& options, Measurements const& measurements,
             IterationCallback callback) {
    models::validate_hermiticity(model);
    auto const compiled = models::compile_hamiltonian(model);
    return solve(compiled.op, models::physical_basis(model), options, measurements, std::move(callback));
}

models::CompiledOperator site_terms(models::Hamiltonian const& model, std::string const& op) {
    models::CompiledOperator out;
    out.site_dims = model.site_dims();
    for (int site = 0; site < model.num_sites(); ++site) {
        auto const basis = model.site_for(site);
        if (!basis->has_op(op)) continue;
        append_column(out, model, models::TermSum::located(basis, op, site));
    }
    return out;
}

models::CompiledOperator pair_terms(models::Hamiltonian const& model, std::string const& first,
                                    std::string const& second) {
    models::CompiledOperator out;
    out.site_dims = model.site_dims();
    for (int i = 0; i < model.num_sites(); ++i) {
        auto const a = model.site_for(i);
        if (!a->has_op(first)) continue;
        for (int j = 0; j < model.num_sites(); ++j) {
            auto const b = model.site_for(j);
            if (!b->has_op(second)) continue;
            append_column(out, model, models::TermSum::located(a, first, i) * models::TermSum::located(b, second, j));
        }
    }
    return out;
}

} // namespace alps::dmrg::classic
