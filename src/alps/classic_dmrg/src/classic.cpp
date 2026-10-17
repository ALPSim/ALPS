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
// The time evolution (evolve, Engine::evolve and the bond matrices) is
// written here; its engine side is src/tdmrg.hpp.
// Written here, not taken from alpsx (it drove dmtk from catalog operators
// through model_helper): the search for a set of fermionic carrier charges,
// accepting two-site factors in either order, normalizing factor matrices,
// combining terms that differ by a scale or by their projector expansion
// (simplify), the dense Hermiticity check, and measuring every component of
// a column on its own instead of applying the summed operator.
#include <alps/dmrg/classic.hpp>
#include "invocation.hpp"
#include "mps.hpp"
#include "tdmrg.hpp"

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
// Matrix entries on a grid of 2^-40, orderable. Keys identify normalized
// factors (largest entry one), so entries that differ by rounding share a key;
// a key that misses only keeps two equal terms apart.
using Key = std::vector<std::pair<long long, long long>>;

Key key_of(std::vector<Complex> const& matrix) {
    constexpr double grid = 1099511627776.; // 2^40
    Key key;
    key.reserve(matrix.size());
    for (auto const& v : matrix) key.emplace_back(std::llround(v.real() * grid), std::llround(v.imag() * grid));
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
    std::vector<Complex> matrix; // row-major, row = output state; leading entry one
    std::vector<long> delta;     // [charge]: twice (output - input)
    bool odd = false;
    Complex scale{1.};           // taken out of the matrix into the term coefficient
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
    // Time evolution: the final Hamiltonian and the drives, merged on the same
    // charges, and the index in `measurements` of the final Hamiltonian's
    // expectation value, which is reported as the energy.
    std::vector<ReadyTerm> final_hamiltonian;
    Complex final_scalar; // the final Hamiltonian's identity part (a phase in real time)
    std::vector<std::vector<ReadyTerm>> drives;
    std::optional<std::size_t> energy_plan;
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
    // The engine derives fermion signs from the integer quantum numbers it
    // flags fermionic: the parity of their sum is the sign. Use one charge
    // whose parity is the state parity everywhere ("N" first), else the
    // smallest set of charges whose summed parity is (the catalog's fermion
    // basis carries Nup and Ndown and no N); the dmrg application registered
    // every fermionic catalog quantum number this way.
    std::vector<std::string> names{"N"};
    for (auto const& site : basis.sites)
        for (auto const& [name, values] : site.quantum_numbers)
            if (std::find(names.begin(), names.end(), name) == names.end()) names.push_back(name);
    std::vector<std::string> candidates;
    for (auto const& name : names) {
        bool present = false, integral = true;
        for (auto const& site : basis.sites) {
            auto found = site.quantum_numbers.find(name);
            if (found == site.quantum_numbers.end()) continue;
            present = true;
            for (double q : found->second)
                if (!std::isfinite(q) || std::abs(q - std::round(q)) > 1e-9) integral = false;
        }
        if (present && integral) candidates.push_back(name);
    }
    auto carries_parity = [&](std::vector<std::size_t> const& chosen) {
        for (auto const& site : basis.sites)
            for (std::size_t state = 0; state < site.parity.size(); ++state) {
                long long sum = 0;
                for (auto c : chosen) {
                    auto found = site.quantum_numbers.find(candidates[c]);
                    if (found != site.quantum_numbers.end()) sum += std::llround(found->second.at(state));
                }
                if ((sum % 2 != 0) != (site.parity[state] < 0.)) return false;
            }
        return true;
    };
    std::size_t const limit = std::min<std::size_t>(candidates.size(), QN_MAX_SIZE);
    for (std::size_t size = 1; size <= limit; ++size) {
        // Combinations of `size` candidates in lexicographic order.
        std::vector<std::size_t> chosen(size);
        for (std::size_t k = 0; k < size; ++k) chosen[k] = k;
        while (true) {
            if (carries_parity(chosen)) {
                for (auto c : chosen) {
                    auto const& name = candidates[c];
                    auto existing = std::find_if(charges.begin(), charges.end(), [&](Charge const& q) { return q.name == name; });
                    if (existing != charges.end()) existing->fermionic = true;
                    else charges.push_back(Charge{name, true, false, 0});
                }
                return charges;
            }
            std::size_t k = size;
            while (k > 0 && chosen[k - 1] == candidates.size() - size + k - 1) --k;
            if (k == 0) break;
            ++chosen[k - 1];
            for (std::size_t m = k; m < size; ++m) chosen[m] = chosen[m - 1] + 1;
        }
    }
    throw std::invalid_argument("classic DMRG needs integer charges whose summed parity equals the fermion parity "
                                "on every site (for example the particle number N, or Nup and Ndown)");
}

// Scale a matrix so its leading entry (the first, in row-major order, whose
// magnitude reaches the largest) is one, and return the scale. Proportional
// factors then share one matrix however a producer distributes scalars
// between coefficients and factors (the catalog keeps fermion signs in the
// factor matrices, ALPS::models in the coefficients).
Complex normalize(std::vector<Complex>& matrix) {
    double largest = 0.;
    for (auto const& v : matrix) largest = std::max(largest, std::abs(v));
    if (largest == 0.) return Complex{};
    Complex pivot;
    for (auto const& v : matrix)
        if (std::abs(v) >= largest * (1. - 1e-12)) { pivot = v; break; }
    for (auto& v : matrix) v /= pivot;
    return pivot;
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
    for (auto& [key, matrix] : groups) {
        Complex const scale = normalize(matrix);
        components.push_back(Factor{factor.site, std::move(matrix), key.first, key.second < 0, scale});
    }
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
    // Local factors on distinct sites form a tensor product: their order in the
    // term carries no sign (the parity string is read from the input state),
    // so take them in ascending site order, as the engine needs them.
    std::vector<models::LocalFactor const*> ordered;
    for (auto const& factor : term.factors) {
        if (factor.site < 0 || factor.site >= problem.num_sites)
            throw std::invalid_argument(std::string("classic DMRG ") + what + " factor site is out of range");
        ordered.push_back(&factor);
    }
    std::sort(ordered.begin(), ordered.end(), [](auto const* a, auto const* b) { return a->site < b->site; });
    if (ordered.size() == 2 && ordered[0]->site == ordered[1]->site)
        throw std::invalid_argument(std::string("classic DMRG ") + what + " factors must be on distinct sites");
    std::vector<std::vector<Factor>> components;
    for (auto const* factor : ordered) {
        components.push_back(split(*factor, problem.sites[problem.site_index[factor->site]], problem.charges.size()));
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
        Complex coefficient = term.coefficient;
        for (auto& f : factors) {
            coefficient *= f.scale;
            f.scale = 1.;
        }
        out.push_back(ReadyTerm{coefficient, std::move(factors)});
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

// Sum like terms, then combine terms that share one factor (or, on one site,
// a charge change): c1 A1 x B + c2 A2 x B = (c1 A1 + c2 A2) x B when A1 and A2
// change the charges and the parity alike. Catalog operators arrive as one
// term per local matrix element of the first site (projector times matrix);
// this restores their compact form, so the engine applies as many terms as
// the operator's natural decomposition has.
std::vector<ReadyTerm> simplify(std::vector<ReadyTerm> const& terms) {
    std::vector<ReadyTerm> current;
    {
        std::map<TermKey, std::size_t> index;
        for (auto const& term : terms) {
            auto [pos, inserted] = index.emplace(term_key(term), current.size());
            if (inserted) current.push_back(term);
            else current[pos->second].coefficient += term.coefficient;
        }
    }
    using Bucket = std::tuple<std::vector<int>, Key, std::vector<long>, bool>;
    // Combine within buckets that share everything but the summed factor.
    auto combine = [](std::vector<ReadyTerm> const& input, auto bucket_of, std::size_t summed) {
        std::vector<ReadyTerm> output;
        std::map<Bucket, std::size_t> index;
        std::vector<std::vector<std::size_t>> members;
        for (std::size_t t = 0; t < input.size(); ++t) {
            auto [pos, inserted] = index.emplace(bucket_of(input[t]), members.size());
            if (inserted) members.emplace_back();
            members[pos->second].push_back(t);
        }
        for (auto const& group : members) {
            ReadyTerm const& first = input[group.front()];
            if (group.size() == 1) { output.push_back(first); continue; }
            ReadyTerm merged = first;
            auto& matrix = merged.factors[summed].matrix;
            std::fill(matrix.begin(), matrix.end(), Complex{});
            double magnitude = 0.;
            for (auto t : group) {
                auto const& term = input[t];
                magnitude = std::max(magnitude, std::abs(term.coefficient));
                auto const& other = term.factors[summed].matrix;
                for (std::size_t e = 0; e < matrix.size(); ++e) matrix[e] += term.coefficient * other[e];
            }
            double largest = 0.;
            for (auto const& v : matrix) largest = std::max(largest, std::abs(v));
            if (largest <= 1e-14 * magnitude) continue; // the members cancel
            merged.coefficient = normalize(matrix);
            output.push_back(std::move(merged));
        }
        return output;
    };
    for (int pass = 0; pass < 8; ++pass) {
        std::size_t const before = current.size();
        auto sites_of = [](ReadyTerm const& t) {
            std::vector<int> sites;
            for (auto const& f : t.factors) sites.push_back(f.site);
            return sites;
        };
        std::vector<ReadyTerm> singles, pairs;
        for (auto& t : current) (t.factors.size() == 2 ? pairs : singles).push_back(std::move(t));
        // On one site: sum every component with the same charge change.
        singles = combine(singles, [&](ReadyTerm const& t) {
            return Bucket{sites_of(t), Key{}, t.factors[0].delta, t.factors[0].odd};
        }, 0);
        // On two sites: sum the factors that multiply one shared factor.
        for (std::size_t shared : {std::size_t{1}, std::size_t{0}}) {
            std::size_t const summed = 1 - shared;
            pairs = combine(pairs, [&](ReadyTerm const& t) {
                return Bucket{sites_of(t), key_of(t.factors[shared].matrix), t.factors[summed].delta, t.factors[summed].odd};
            }, summed);
        }
        current = std::move(singles);
        current.insert(current.end(), std::make_move_iterator(pairs.begin()), std::make_move_iterator(pairs.end()));
        if (current.size() == before) break;
    }
    std::vector<ReadyTerm> kept;
    for (auto& term : current)
        if (term.coefficient != Complex{}) kept.push_back(std::move(term));
    return kept;
}

// Hermiticity, checked numerically on every site and site pair: the even part
// of the operator must equal its adjoint and the odd-odd part its negated
// adjoint (the adjoint of c F_i F_j with a string over [i, j) is
// -conj(c) F_i^+ F_j^+ with the same string).
void check_hermitian(std::vector<ReadyTerm> const& terms, Problem const& problem, double scale,
                     std::string const& what) {
    std::map<std::vector<int>, std::vector<ReadyTerm const*>> groups;
    for (auto const& term : terms) {
        std::vector<int> sites;
        for (auto const& f : term.factors) sites.push_back(f.site);
        groups[sites].push_back(&term);
    }
    for (auto const& [sites, members] : groups) {
        std::vector<std::size_t> dims;
        for (int s : sites) dims.push_back(static_cast<std::size_t>(problem.sites[problem.site_index[s]].dim));
        std::size_t const second = dims.size() == 2 ? dims[1] : 1;
        std::size_t const dim = dims[0] * second;
        std::vector<Complex> even(dim * dim), odd(dim * dim);
        for (auto const* term : members) {
            auto& target = term->factors.size() == 2 && term->factors[0].odd ? odd : even;
            auto const& a = term->factors[0].matrix;
            for (std::size_t o1 = 0; o1 < dims[0]; ++o1)
                for (std::size_t i1 = 0; i1 < dims[0]; ++i1) {
                    Complex const x = term->coefficient * a[o1 * dims[0] + i1];
                    if (x == Complex{}) continue;
                    if (dims.size() == 1) { target[o1 * dim + i1] += x; continue; }
                    auto const& b = term->factors[1].matrix;
                    for (std::size_t o2 = 0; o2 < second; ++o2)
                        for (std::size_t i2 = 0; i2 < second; ++i2)
                            target[(o1 * second + o2) * dim + (i1 * second + i2)] += x * b[o2 * second + i2];
                }
        }
        for (std::size_t r = 0; r < dim; ++r)
            for (std::size_t c = 0; c < dim; ++c)
                if (std::abs(even[r * dim + c] - std::conj(even[c * dim + r])) > 1e-10 * scale
                    || std::abs(odd[r * dim + c] + std::conj(odd[c * dim + r])) > 1e-10 * scale) {
                    std::string where;
                    for (int s : sites) where += (where.empty() ? "" : ", ") + std::to_string(s);
                    throw std::invalid_argument("classic DMRG " + what + " is not Hermitian (term on sites " + where + ")");
                }
    }
}

// Simplify the expanded terms, require that every surviving term conserves
// the sector's charges and that the sum is Hermitian. `what` names the
// operator in errors; a time evolution's operators may have no term.
std::vector<ReadyTerm> merge_hamiltonian(std::vector<ReadyTerm> const& terms, Problem const& problem,
                                         std::string const& what = "Hamiltonian", bool allow_empty = false) {
    auto merged = simplify(terms);
    double scale = 0.;
    for (auto const& term : merged) scale = std::max(scale, std::abs(term.coefficient));
    std::vector<ReadyTerm> kept;
    for (auto& term : merged) {
        if (std::abs(term.coefficient) <= 1e-14 * scale) continue;
        for (std::size_t c = 0; c < problem.charges.size(); ++c) {
            if (!problem.charges[c].conserved) continue;
            long delta = 0;
            for (auto const& f : term.factors) delta += f.delta[c];
            if (delta != 0 && std::abs(term.coefficient) > 1e-12 * scale)
                throw std::invalid_argument("classic DMRG " + what + " does not conserve the requested charge: "
                                            + problem.charges[c].name);
            if (delta != 0) { term.coefficient = Complex{}; break; }
        }
        if (term.coefficient != Complex{}) kept.push_back(std::move(term));
    }
    if (kept.empty() && !allow_empty) throw std::invalid_argument("classic DMRG " + what + " has no term");
    check_hermitian(kept, problem, scale, what);
    return kept;
}

// One requested operator as measurement columns: the whole operator in one
// column (expectation) or one column per term (terms).
Plan lower_plan(std::string const& name, models::CompiledOperator const& op, MeasurementKind kind,
                Problem const& problem, std::vector<int> const& site_dims) {
    if (op.site_dims != site_dims)
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
        col.components = simplify(col.components); // the column reports their sum
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
    return entry;
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
            problem.measurements.push_back(lower_plan(name, op, kind, problem, hamiltonian.site_dims));
        }
    };
    plan(measurements.expectation, MeasurementKind::expectation);
    plan(measurements.terms, MeasurementKind::terms);
    return problem;
}

// A warm-up configuration names a state of every site, and its conserved
// charges add up to the sector.
void check_warmup_configuration(Problem const& problem, Options const& options) {
    auto const& configuration = options.warmup_configuration;
    if (configuration.empty()) return;
    if (configuration.size() != static_cast<std::size_t>(problem.num_sites))
        throw std::invalid_argument("classic DMRG warmup_configuration needs one state per site: "
                                    + std::to_string(configuration.size()) + " for "
                                    + std::to_string(problem.num_sites) + " sites");
    std::vector<long> sum(problem.charges.size(), 0);
    for (int i = 0; i < problem.num_sites; ++i) {
        auto const& site = problem.sites[static_cast<std::size_t>(problem.site_index[static_cast<std::size_t>(i)])];
        auto const state = configuration[static_cast<std::size_t>(i)];
        if (state >= static_cast<std::size_t>(site.dim))
            throw std::invalid_argument("classic DMRG warmup_configuration names state " + std::to_string(state)
                                        + " of site " + std::to_string(i) + ", which has "
                                        + std::to_string(site.dim) + " states");
        for (std::size_t c = 0; c < problem.charges.size(); ++c) sum[c] += site.twice[c][state];
    }
    auto half = [](long doubled) {
        return doubled % 2 == 0 ? std::to_string(doubled / 2) : std::to_string(doubled) + "/2";
    };
    for (std::size_t c = 0; c < problem.charges.size(); ++c)
        if (problem.charges[c].conserved && sum[c] != problem.charges[c].target)
            throw std::invalid_argument("classic DMRG warmup_configuration has " + problem.charges[c].name + " = "
                                        + half(sum[c]) + ", but the sector asks for "
                                        + half(problem.charges[c].target));
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

// ----------------------------------------------------------- time evolution

void validate_evolution(Evolution const& evolution, Options const& options) {
    if (!std::isfinite(evolution.time_step) || evolution.time_step <= 0.)
        throw std::invalid_argument("classic DMRG time_step must be finite and positive");
    if (evolution.steps > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("classic DMRG time evolution step count exceeds the integer limit");
    bounded_count(evolution.states, "evolution states");
    if (evolution.truncation_error && (!std::isfinite(*evolution.truncation_error)
                                       || *evolution.truncation_error <= 0. || *evolution.truncation_error >= 1.))
        throw std::invalid_argument("classic DMRG evolution truncation_error must be finite and between zero and one");
    for (auto step : evolution.measure_steps)
        if (step > evolution.steps)
            throw std::invalid_argument("classic DMRG measure step " + std::to_string(step) + " is after the last step");
    for (std::size_t k = 0; k < evolution.drives.size(); ++k) {
        auto const& coefficients = evolution.drives[k].coefficients;
        if (coefficients.size() != evolution.steps)
            throw std::invalid_argument("classic DMRG drive " + std::to_string(k) + " needs one coefficient per step");
        for (double c : coefficients)
            if (!std::isfinite(c))
                throw std::invalid_argument("classic DMRG drive " + std::to_string(k) + " has a nonfinite coefficient");
    }
    if (evolution.method != EvolutionMethod::trotter && evolution.method != EvolutionMethod::targeting)
        throw std::invalid_argument("unknown classic DMRG evolution method");
    if (evolution.method == EvolutionMethod::targeting) {
        bounded_count(evolution.half_sweeps, "half_sweeps");
        if (!evolution.drives.empty())
            throw std::invalid_argument("classic DMRG time-step targeting takes no drives (its block Hamiltonians "
                                        "are built once); use the Trotter method");
    }
    if (options.num_eigenvalues != 1)
        throw std::invalid_argument("classic DMRG time evolution evolves one state: num_eigenvalues must be one");
    if (options.keep_mps)
        throw std::invalid_argument("classic DMRG time evolution keeps its states through Evolution::keep_mps, "
                                    "not Options::keep_mps");
}

// Trotter gates act on neighbouring sites in chain order only.
void require_bonds(std::vector<ReadyTerm> const& terms, std::string const& what) {
    for (auto const& term : terms)
        if (term.factors.size() == 2 && term.factors[1].site != term.factors[0].site + 1)
            throw std::invalid_argument("classic DMRG time evolution applies nearest-neighbour bond gates in chain "
                                        "order; the " + what + " couples sites " + std::to_string(term.factors[0].site)
                                        + " and " + std::to_string(term.factors[1].site));
}

// Merge the final Hamiltonian and the drives on the initial problem's
// charges (each must conserve the sector and be Hermitian, and may be empty)
// and add the final Hamiltonian's expectation value as the last measurement.
void prepare_evolution(Problem& problem, models::CompiledOperator const& final_hamiltonian,
                       Evolution const& evolution, std::vector<int> const& site_dims) {
    auto merged = [&](models::CompiledOperator const& op, std::string const& what, Complex& scalar) {
        if (op.site_dims != site_dims)
            throw std::invalid_argument("classic DMRG " + what + " is not on the initial Hamiltonian's sites");
        std::vector<ReadyTerm> expanded;
        for (auto const& term : op.terms) expand(term, problem, expanded, scalar, what.c_str());
        auto terms = merge_hamiltonian(expanded, problem, what, true);
        if (evolution.method == EvolutionMethod::trotter) require_bonds(terms, what);
        return terms;
    };
    // The identity parts are a phase in real time; neither method applies them.
    problem.final_hamiltonian = merged(final_hamiltonian, "final Hamiltonian", problem.final_scalar);
    for (std::size_t k = 0; k < evolution.drives.size(); ++k) {
        Complex scalar;
        problem.drives.push_back(merged(evolution.drives[k].op, "drive " + std::to_string(k), scalar));
    }
    problem.energy_plan = problem.measurements.size();
    problem.measurements.push_back(
        lower_plan("energy", final_hamiltonian, MeasurementKind::expectation, problem, site_dims));
}

// The two-site matrix of every bond (b, b + 1) in the physical basis order,
// [out][in] with the pair index s * d(b + 1) + t: the bond's two-site terms
// plus each site's one-site terms divided evenly among its bonds (an end
// site's go wholly to its one bond). A product of two odd factors carries
// the parity of site b in the input state, the strict convention's string
// over [b, b + 1), as ALPS::ed applies it.
std::vector<std::vector<Complex>> bond_matrices(std::vector<ReadyTerm> const& terms, Problem const& problem) {
    int const length = problem.num_sites;
    auto dim = [&](int site) {
        return static_cast<std::size_t>(problem.sites[static_cast<std::size_t>(problem.site_index[static_cast<std::size_t>(site)])].dim);
    };
    std::vector<std::vector<Complex>> bonds(static_cast<std::size_t>(length - 1));
    for (int b = 0; b + 1 < length; ++b) {
        std::size_t const n = dim(b) * dim(b + 1);
        bonds[static_cast<std::size_t>(b)].assign(n * n, Complex{});
    }
    for (auto const& term : terms) {
        int const i = term.factors[0].site;
        auto const& a = term.factors[0].matrix;
        if (term.factors.size() == 1) {
            std::vector<int> touching;
            if (i > 0) touching.push_back(i - 1);
            if (i + 1 < length) touching.push_back(i);
            Complex const weight = term.coefficient / static_cast<double>(touching.size());
            for (int b : touching) {
                std::size_t const d1 = dim(b), d2 = dim(b + 1), n = d1 * d2;
                auto& m = bonds[static_cast<std::size_t>(b)];
                for (std::size_t s = 0; s < d1; ++s)
                    for (std::size_t t = 0; t < d2; ++t)
                        if (b == i) {
                            for (std::size_t out = 0; out < d1; ++out)
                                m[(out * d2 + t) * n + s * d2 + t] += weight * a[out * d1 + s];
                        } else {
                            for (std::size_t out = 0; out < d2; ++out)
                                m[(s * d2 + out) * n + s * d2 + t] += weight * a[out * d2 + t];
                        }
            }
            continue;
        }
        auto const& c = term.factors[1].matrix;
        std::size_t const d1 = dim(i), d2 = dim(i + 1), n = d1 * d2;
        auto const& parity = problem.sites[static_cast<std::size_t>(problem.site_index[static_cast<std::size_t>(i)])].parity;
        bool const string = term.factors[0].odd && term.factors[1].odd;
        auto& m = bonds[static_cast<std::size_t>(i)];
        for (std::size_t s = 0; s < d1; ++s)
            for (std::size_t t = 0; t < d2; ++t) {
                Complex const sign = term.coefficient * (string ? parity[s] : 1.);
                for (std::size_t s_out = 0; s_out < d1; ++s_out) {
                    Complex const x = sign * a[s_out * d1 + s];
                    if (x == Complex{}) continue;
                    for (std::size_t t_out = 0; t_out < d2; ++t_out)
                        m[(s_out * d2 + t_out) * n + s * d2 + t] += x * c[t_out * d2 + t];
                }
            }
    }
    return bonds;
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
    dmtk::Hami<T> final_hami; // time-step targeting's final Hamiltonian
    std::vector<std::size_t> qn_index; // per charge
    int conserved_mask = 0;
    std::map<std::pair<int, Key>, std::string> operator_names;
    std::vector<std::vector<std::vector<std::size_t>>> layout; // [plan][column] -> measured term indices
    bool two_site_measurements = false;
    bool evolving = false; // the time evolution's iterations solve nothing

    Engine(Problem const& p, Options const& o, IterationCallback cb)
        : problem(p), options(o), callback(std::move(cb)), lattice(p.num_sites, dmtk::OBC), hami(lattice) {}

    void observe(dmtk::System<T> const& s) {
        // An evolution's iterations report the final Hamiltonian's energy
        // where time-step targeting computes it, else NaN.
        double const energy = evolving ? s.energy[0] + problem.final_scalar.real() : s.energy[0];
        Iteration row{s.get_dir() == dmtk::RIGHT2LEFT ? Direction::right_to_left : Direction::left_to_right,
                      static_cast<std::size_t>(s.get_iter()), energy, s.truncation_error(), s.entropy()};
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
                // Every component is measured on its own and the column reports
                // their sum: on-site terms as weighted terms, two-site products
                // as composite operators (alpsx dmrg.hpp). A composite built
                // from a whole sum would be applied to the state at every
                // measurement, as the engine applies a Hamiltonian.
                for (auto const& component : plan.columns[c].components) {
                    dmtk::Term<T> term = make_term(component);
                    if (component.factors.size() == 1) {
                        measured += term;
                    } else {
                        measured += dmtk::BasicOp<T>(term);
                        two_site_measurements = true;
                    }
                    layout[p][c].push_back(measured.size() - 1);
                }
            }
        }
    }

    // The controls of the initial DMRG, up to the system's start (shared by
    // solve and evolve).
    void configure(dmtk::System<T>& S, bool store_products) {
        S.set_store_products(store_products);
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
    }

    // The warm-up targets of Options::warmup_configuration: entry k holds the
    // charges of its first k sites (System::warmup_loop reads entry k for
    // its k-site chain).
    dmtk::Vector<dmtk::QN> warmup_sectors() const {
        dmtk::Vector<dmtk::QN> sectors(static_cast<std::size_t>(problem.num_sites) + 1);
        std::vector<long> sum(problem.charges.size(), 0);
        for (int i = 0; i < problem.num_sites; ++i) {
            auto const& site = problem.sites[static_cast<std::size_t>(problem.site_index[static_cast<std::size_t>(i)])];
            auto const state = options.warmup_configuration[static_cast<std::size_t>(i)];
            for (std::size_t c = 0; c < problem.charges.size(); ++c) {
                sum[c] += site.twice[c][state];
                sectors[static_cast<std::size_t>(i) + 1][static_cast<int>(qn_index[c])].set_half(static_cast<short>(sum[c]));
            }
        }
        return sectors;
    }

    ResolvedOptions resolved(dmtk::System<T> const& S) const {
        ResolvedOptions out;
        out.sweeps = options.sweeps;
        out.warmup_states = options.warmup_states;
        out.warmup_configuration = options.warmup_configuration;
        out.num_eigenvalues = options.num_eigenvalues;
        out.truncation_error = options.truncation_error;
        out.truncation_max_states = options.truncation_max_states;
        out.lanczos_tolerance = S.lanczos_tolerance();
        out.arithmetic = std::is_same_v<T, double> ? Arithmetic::real : Arithmetic::complex;
        out.sector = options.sector;
        return out;
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
        configure(S, options.num_eigenvalues > 1 && two_site_measurements);
        if (options.warmup_configuration.empty()) {
            S.run(static_cast<int>(options.warmup_states));
        } else {
            // System::run on the configuration's prefix charges: a new System
            // reads custom warm-up targets.
            S.warmup_loop(options.warmup_states, warmup_sectors());
            S.resume(1, dmtk::RIGHT2LEFT, 1);
        }
        if (S.store_products()) S.set_full_sweep(true);
        S.corr += measured;
        S.final_sweep(options.sweeps.back().left_to_right, dmtk::RIGHT2LEFT, 1, true);
        // The unstored path measures while sweeping and measure() refreshes
        // what the final configuration can reach; the stored path measures
        // nothing during the sweep, so this call takes every value.
        S.measure();

        Result result;
        // The tensors come from the scratch files, which the invocation
        // removes when the solve returns, and from the closing wave
        // functions, before the reading below overwrites S.gs.
        if (options.keep_mps) {
            detail::MpsExport<T> const exporter(S);
            result.mps.push_back(exporter.state(S.gs));
            for (std::size_t i = 1; i < S._target.size(); ++i) result.mps.push_back(exporter.state(S._target[i]));
        }
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
        result.options = resolved(S);
        return result;
    }

    // A bond matrix of bond_matrices() in the engine's charge-sorted site
    // bases (State::i1 is the kit index of a sorted state, as in
    // local_operator), real when the engine is.
    std::vector<Complex> sorted_bond(std::vector<Complex> const& kit, int b) const {
        auto const& first = blocks[static_cast<std::size_t>(problem.site_index[static_cast<std::size_t>(b)])].basis();
        auto const& second = blocks[static_cast<std::size_t>(problem.site_index[static_cast<std::size_t>(b + 1)])].basis();
        std::size_t const d1 = first.size(), d2 = second.size(), n = d1 * d2;
        std::vector<Complex> out(n * n);
        for (std::size_t po = 0; po < d1; ++po)
            for (std::size_t qo = 0; qo < d2; ++qo)
                for (std::size_t pi = 0; pi < d1; ++pi)
                    for (std::size_t qi = 0; qi < d2; ++qi) {
                        Complex const v = kit[(first[po].i1 * d2 + second[qo].i1) * n + first[pi].i1 * d2 + second[qi].i1];
                        out[(po * d2 + qo) * n + pi * d2 + qi] = scalar(v, "evolution Hamiltonian");
                    }
        return out;
    }

    // Prepare the ground state of the initial Hamiltonian, then apply the
    // second-order Trotter steps exp(-dt A/2) exp(-dt B) exp(-dt A/2) (A the
    // even bonds, B the odd ones; dt imaginary in real time) as half sweeps
    // of bond gates (detail::EvolvingSystem). Adjacent half steps of A are
    // merged into one half sweep of their gate products, which is exact,
    // drives included; a measured step closes its A half step and measures
    // in an identity half sweep.
    EvolutionResult evolve(Evolution const& evolution) {
        register_charges();
        build_blocks();
        build_hamiltonian();
        build_measurements();
        bool const targeting = evolution.method == EvolutionMethod::targeting;
        if (targeting) {
            // Named apart from the initial Hamiltonian ("H"), whose block
            // Hamiltonians the preparation leaves in the blocks.
            final_hami = dmtk::Hami<T>(hami);
            final_hami.clear();
            final_hami.set_name("HF");
            for (auto const& term : problem.final_hamiltonian) final_hami += make_term(term);
            final_hami.reorder_terms();
        }

        detail::EvolvingSystem<T> S(hami, lattice);
        configure(S, false);
        if (options.warmup_configuration.empty())
            S.prepare(static_cast<int>(options.warmup_states));
        else
            S.prepare(static_cast<int>(options.warmup_states), warmup_sectors());

        EvolutionResult result;
        auto& initial = result.initial;
        for (std::size_t i = 0; i < S.energy.size(); ++i) initial.energies.push_back(S.energy[i] + problem.scalar.real());
        initial.iterations = std::move(rows);
        rows.clear();
        initial.final_truncation_error = S.truncation_error();
        initial.options = resolved(S);

        S.corr += measured;
        S.init_correlations();
        if (targeting)
            S.begin_targeting(final_hami, evolution.time_step, evolution.imaginary, evolution.states,
                              evolution.truncation_error);
        else
            S.begin_evolution(evolution.states, evolution.truncation_error);
        evolving = true;

        // The Trotter method's bond Hamiltonians of the final Hamiltonian and
        // of every drive.
        int const length = problem.num_sites;
        std::size_t const bonds = static_cast<std::size_t>(length - 1);
        auto sorted = [&](std::vector<ReadyTerm> const& terms) {
            auto kit = bond_matrices(terms, problem);
            for (std::size_t b = 0; b < bonds; ++b) kit[b] = sorted_bond(kit[b], static_cast<int>(b));
            return kit;
        };
        std::vector<std::vector<Complex>> final_bonds;
        std::vector<std::vector<std::vector<Complex>>> drive_bonds;
        if (!targeting) {
            final_bonds = sorted(problem.final_hamiltonian);
            for (auto const& drive : problem.drives) drive_bonds.push_back(sorted(drive));
        }
        auto dim = [&](int site) {
            return static_cast<std::size_t>(problem.sites[static_cast<std::size_t>(problem.site_index[static_cast<std::size_t>(site)])].dim);
        };
        using Gates = std::vector<detail::BondGate<T>>; // one per bond; empty is the identity
        // The gates of step `step` over a time t on the bonds b with b % 2 == parity.
        std::map<std::pair<double, int>, Gates> cache; // without drives every step has the same gates
        auto gates = [&](std::size_t step, double t, int parity) {
            if (drive_bonds.empty()) {
                auto found = cache.find({t, parity});
                if (found != cache.end()) return found->second;
            }
            Gates out(bonds);
            for (std::size_t b = static_cast<std::size_t>(parity); b < bonds; b += 2) {
                auto h = final_bonds[b];
                for (std::size_t k = 0; k < drive_bonds.size(); ++k) {
                    double const c = evolution.drives[k].coefficients[step];
                    for (std::size_t x = 0; x < h.size(); ++x) h[x] += c * drive_bonds[k][b][x];
                }
                out[b] = detail::exponential<T>(h, dim(static_cast<int>(b)), dim(static_cast<int>(b) + 1), t,
                                                evolution.imaginary);
            }
            if (drive_bonds.empty()) cache.emplace(std::make_pair(t, parity), out);
            return out;
        };
        Gates const identity(bonds);

        std::size_t direction = dmtk::RIGHT2LEFT;
        auto sweep = [&](Gates const& g, bool measure) {
            std::vector<detail::BondGate<T> const*> pointers;
            for (auto const& gate : g) pointers.push_back(gate.identity() ? nullptr : &gate);
            if (measure) // every value must be taken again in this half sweep
                for (auto& term : S.corr) term.set_value(T(std::numeric_limits<double>::quiet_NaN()));
            double const discarded = S.half_sweep(direction, pointers, measure);
            direction = direction == dmtk::RIGHT2LEFT ? dmtk::LEFT2RIGHT : dmtk::RIGHT2LEFT;
            return discarded;
        };

        for (std::size_t p = 0; p < problem.measurements.size(); ++p)
            if (p != problem.energy_plan)
                result.measurements.emplace(problem.measurements[p].name,
                                            Measurement{problem.measurements[p].kind, {}, {}});
        auto record = [&](std::size_t step) {
            result.steps.push_back(step);
            result.times.push_back(static_cast<double>(step) * evolution.time_step);
            for (std::size_t p = 0; p < problem.measurements.size(); ++p) {
                auto const& plan = problem.measurements[p];
                std::vector<double> row;
                for (std::size_t c = 0; c < plan.columns.size(); ++c) {
                    double value = plan.columns[c].scalar.real();
                    for (auto index : layout[p][c]) {
                        T const v = S.corr[index].value();
                        if (!std::isfinite(std::real(v)) || !std::isfinite(std::imag(v)))
                            throw std::logic_error("classic DMRG time evolution did not measure " + plan.name
                                                   + " in its measuring half sweep");
                        value += std::real(v);
                    }
                    row.push_back(value);
                }
                if (p == problem.energy_plan) result.energies.push_back(row.at(0));
                else result.measurements.at(plan.name).values.push_back(std::move(row));
            }
            // The tensors come from this half sweep's rotations, on disk until
            // the invocation returns, and from the closing wave function.
            if (evolution.keep_mps) {
                detail::MpsExport<T> const exporter(S);
                result.mps.push_back(exporter.state(S.gs));
            }
        };

        std::set<std::size_t> measured_steps(evolution.measure_steps.begin(), evolution.measure_steps.end());
        if (measured_steps.empty())
            for (std::size_t n = 0; n <= evolution.steps; ++n) measured_steps.insert(n);
        double const dt = evolution.time_step;
        result.truncation_errors.assign(evolution.steps + 1, 0.);
        if (targeting) {
            auto targeted = [&](bool advance, bool measure) {
                if (measure)
                    for (auto& term : S.corr) term.set_value(T(std::numeric_limits<double>::quiet_NaN()));
                double const lost = S.targeting_half_sweep(direction, advance, measure);
                direction = direction == dmtk::RIGHT2LEFT ? dmtk::LEFT2RIGHT : dmtk::RIGHT2LEFT;
                return lost;
            };
            // Two half sweeps build the blocks of the final Hamiltonian and of
            // the measured operators, one side each; the second measures the
            // initial state.
            bool const initial_measured = measured_steps.count(0) != 0;
            double lost = targeted(false, false);
            lost += targeted(false, initial_measured);
            result.truncation_errors[0] = lost;
            if (initial_measured) record(0);
            for (std::size_t n = 0; n < evolution.steps; ++n) {
                lost = 0.;
                for (std::size_t h = 0; h < evolution.half_sweeps; ++h)
                    lost += targeted(h + 1 == evolution.half_sweeps, false);
                bool const measure = measured_steps.count(n + 1) != 0;
                if (measure) lost += targeted(false, true);
                result.truncation_errors[n + 1] = lost;
                if (measure) record(n + 1);
            }
        } else {
            if (measured_steps.count(0)) {
                // The left blocks of the preparation hold no measured operator;
                // a first half sweep rebuilds the right blocks with them and the
                // measuring one the left blocks.
                double discarded = sweep(identity, false);
                discarded += sweep(identity, true);
                result.truncation_errors[0] = discarded;
                record(0);
            }
            Gates pending; // the closing even half step of the previous step, not yet applied
            for (std::size_t n = 0; n < evolution.steps; ++n) {
                Gates half = gates(n, 0.5 * dt, 0);
                if (pending.empty()) {
                    result.truncation_errors[n + 1] += sweep(half, false);
                } else {
                    Gates merged(bonds);
                    for (std::size_t b = 0; b < bonds; ++b) merged[b] = detail::compose(half[b], pending[b]);
                    result.truncation_errors[n] += sweep(merged, false);
                }
                double discarded = sweep(gates(n, dt, 1), false);
                pending = std::move(half);
                bool const measure = measured_steps.count(n + 1) != 0;
                if (measure || n + 1 == evolution.steps) {
                    discarded += sweep(pending, false);
                    pending.clear();
                }
                if (measure) discarded += sweep(identity, true);
                result.truncation_errors[n + 1] += discarded;
                if (measure) record(n + 1);
            }
        }
        for (std::size_t p = 0; p < problem.measurements.size(); ++p) {
            if (p == problem.energy_plan) continue;
            auto& measurement = result.measurements.at(problem.measurements[p].name);
            for (auto const& column : problem.measurements[p].columns) measurement.sites.push_back(column.sites);
        }
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
    check_warmup_configuration(problem, options);
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

EvolutionResult evolve(models::CompiledOperator const& initial, models::CompiledOperator const& final_hamiltonian,
                       models::PhysicalBasis const& basis, Evolution const& evolution, Options const& options,
                       Measurements const& measurements, IterationCallback callback) {
    validate_options(options);
    validate_evolution(evolution, options);
    Problem problem = prepare(initial, basis, options, measurements);
    prepare_evolution(problem, final_hamiltonian, evolution, initial.site_dims);
    check_sector(problem, 1);
    check_warmup_configuration(problem, options);
    dmtk::InvocationGuard invocation(options.scratch_parent, options.log);
    // Every engine object is destroyed before the guard restores the engine state.
    if (evolution.imaginary && options.arithmetic == Arithmetic::real) {
        Engine<double> engine(problem, options, std::move(callback));
        return engine.evolve(evolution);
    }
    Engine<Complex> engine(problem, options, std::move(callback));
    return engine.evolve(evolution);
}

EvolutionResult evolve(models::Hamiltonian const& initial, models::Hamiltonian const& final_hamiltonian,
                       Evolution const& evolution, Options const& options, Measurements const& measurements,
                       IterationCallback callback) {
    models::validate_hermiticity(initial);
    models::validate_hermiticity(final_hamiltonian);
    auto const basis = models::physical_basis(initial);
    models::require_same_basis(models::physical_basis(final_hamiltonian), basis);
    auto const first = models::compile_hamiltonian(initial);
    auto const second = models::compile_hamiltonian(final_hamiltonian);
    return evolve(first.op, second.op, basis, evolution, options, measurements, std::move(callback));
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

// The molecular solver shares this translation unit's engine helpers.
#include "molecular.ipp"
