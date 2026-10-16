// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// WIP/classic-dmrg-chemistry at d4384a7f8, file src/dmrg/molecular.ipp.
// Included at the end of classic.cpp, as alpsx included it in its classic
// translation unit, so the engine and its invocation-global state stay
// private to the component. Adaptations: ALPS::models headers replace
// alpsx's alps/solve/models; the iteration callback is called directly (this
// port's InvocationGuard has no caller-context wrapper, as in solve()); the
// charges are registered on a fresh registry and the mask is built from the
// indices they receive; the site basis is built from its N and Sz values
// only; the solver names molecular_detail::Site, since classic.cpp has a
// Site of its own. The validation, Hamiltonian lowering, measurement plan
// and engine calls follow alpsx line by line.
#include <alps/dmrg/molecular.hpp>
#include <alps/models/hamiltonian.hpp>
#include <alps/models/site.hpp>
#include <alps/models/terms.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace alps::dmrg::classic {
namespace molecular_detail {

std::size_t product_size(std::size_t a, std::size_t b) {
    if (a && b > std::vector<double>().max_size() / a)
        throw std::length_error("classic molecular tensor is too large");
    return a * b;
}

void validate(MolecularInput const& input, Options const& options) {
    auto const n = input.orbitals;
    if (!n || n > static_cast<std::size_t>(std::numeric_limits<short>::max()) / 4)
        throw std::invalid_argument("classic molecular orbital count must be positive and representable");
    if (input.alpha_electrons > n || input.beta_electrons > n)
        throw std::invalid_argument("classic molecular electron counts must lie between zero and the orbital count");
    auto const n2 = product_size(n, n), n4 = product_size(n2, n2);
    if (input.one_body.size() != n2 || input.two_body.size() != n4)
        throw std::invalid_argument("classic molecular integral shapes must be (n,n) and (n,n,n,n)");
    auto finite = [](double x) { return std::isfinite(x); };
    if (!std::all_of(input.one_body.begin(), input.one_body.end(), finite) ||
        !std::all_of(input.two_body.begin(), input.two_body.end(), finite) ||
        !finite(input.core_energy) || !finite(input.integral_cutoff) || input.integral_cutoff < 0.)
        throw std::invalid_argument("classic molecular integrals and core energy must be finite; cutoff must be finite and nonnegative");
    // A relative comparison cannot let a large energy shift hide a broken
    // hopping term. Check the operator actually retained by the explicit
    // cutoff: discarded roundoff entries do not impose false asymmetries.
    auto equal = [cutoff = input.integral_cutoff](double a, double b) {
        if (std::abs(a) < cutoff) a = 0.;
        if (std::abs(b) < cutoff) b = 0.;
        return std::abs(a - b) <= 1.e-12 * std::max(std::abs(a), std::abs(b));
    };
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < n; ++q) {
            if (!equal(input.one_body[p * n + q], input.one_body[q * n + p]))
                throw std::invalid_argument("classic molecular one-body integrals must be symmetric");
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t s = 0; s < n; ++s) {
                    double a = input.two_body[((p * n + q) * n + r) * n + s];
                    if (!equal(a, input.two_body[((q * n + p) * n + r) * n + s]) ||
                        !equal(a, input.two_body[((p * n + q) * n + s) * n + r]) ||
                        !equal(a, input.two_body[((r * n + s) * n + p) * n + q]))
                        throw std::invalid_argument("classic molecular two-body integrals require real-orbital eightfold symmetry");
                }
        }
    bounded_count(options.sweeps.size(), "sweep count");
    if (options.sweeps.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()) / 2)
        throw std::invalid_argument("classic molecular half-sweep count exceeds the integer limit");
    for (auto const& sweep : options.sweeps) {
        bounded_count(sweep.right_to_left, "right_to_left states");
        bounded_count(sweep.left_to_right, "left_to_right states");
    }
    bounded_count(options.warmup_states, "warmup_states");
    if (options.num_eigenvalues != 1 || options.arithmetic != Arithmetic::real || !options.sector.empty())
        throw std::invalid_argument("classic molecular solves require one real target and electron counts from MolecularInput, with options.sector empty");
    if (options.lanczos_tolerance && (!finite(*options.lanczos_tolerance) || *options.lanczos_tolerance <= 0.))
        throw std::invalid_argument("classic molecular lanczos_tolerance must be finite and positive");
    if (options.truncation_error && (!finite(*options.truncation_error) ||
                                    *options.truncation_error <= 0. || *options.truncation_error >= 1.))
        throw std::invalid_argument("classic molecular truncation_error must be between zero and one");
    if (options.truncation_max_states) {
        bounded_count(*options.truncation_max_states, "truncation_max_states");
        if (!options.truncation_error)
            throw std::invalid_argument("classic molecular truncation_max_states requires truncation_error");
    }
}

struct Site {
    models::SiteBasisPtr source;
    dmtk::Basis basis;
    dmtk::Block<double> block;
    // Numeric matrices, not display names, define local operator identity.
    std::map<std::vector<double>, std::size_t> indices;
    std::vector<dmtk::BasicOp<double>> operators;

    Site(models::SiteBasisPtr input, std::size_t n_index, std::size_t sz_index)
        : source(std::move(input)), basis(static_cast<std::size_t>(source->dim)) {
        auto const& n = source->quantum_numbers.at("N");
        auto const& sz = source->quantum_numbers.at("Sz");
        for (int i = 0; i < source->dim; ++i) {
            dmtk::QN q;
            q[static_cast<int>(n_index)] = n[static_cast<std::size_t>(i)];
            q[static_cast<int>(sz_index)] = sz[static_cast<std::size_t>(i)];
            basis[static_cast<std::size_t>(i)] = dmtk::State(static_cast<std::size_t>(i), q);
        }
        basis.reorder();
        block.resize(basis);
        block.set_lattice(dmtk::Lattice(1, dmtk::OBC));
    }

    std::size_t intern(models::LocalOp const& op) {
        std::vector<double> key;
        for (auto value : op.matrix) {
            if (value.imag() != 0. || !std::isfinite(value.real()))
                throw std::invalid_argument("classic molecular operator must be finite and real");
            key.push_back(value.real());
        }
        if (auto found = indices.find(key); found != indices.end()) return found->second;
        auto const id = operators.size();
        std::string name = op.is_identity(0.) ? "I" : "m" + std::to_string(id);
        dmtk::BasicOp<double> native(name.c_str(), 0);
        auto const dim = static_cast<std::size_t>(source->dim);
        dmtk::Matrix<double> matrix(dim, dim);
        bool first = true, diagonal = true;
        for (std::size_t row = 0; row < dim; ++row)
            for (std::size_t col = 0; col < dim; ++col) {
                // The block basis is charge-sorted; i1 is the site-basis index.
                double value = key[basis[row].i1 * dim + basis[col].i1];
                matrix[row][col] = value; // Matrix's [] is row-first; () is column-first.
                if (value == 0.) continue;
                auto delta = basis[row].qn() - basis[col].qn();
                if (!first && delta != native.dqn)
                    throw std::logic_error("classic molecular local operator mixes charge changes");
                native.dqn = delta;
                first = false;
                if (row != col) diagonal = false;
            }
        native.set_fermion(op.parity == models::Parity::odd).set_diagonal(diagonal);
        native.resize(basis);
        native = matrix;
        block.push_back(native);
        operators.push_back(native);
        indices.emplace(std::move(key), id);
        return id;
    }
};

using Word = std::vector<std::pair<int, std::size_t>>;
struct Lowered {
    dmtk::Term<double> term;
    Word word;
};

Lowered lower(models::Term const& input, Site& site) {
    auto term = models::canonicalise_term(input, true, 0.);
    Lowered output;
    if (term.effective_coeff().imag() != 0.)
        throw std::logic_error("classic molecular coefficient is not real");
    output.term.coef() = term.effective_coeff().real();
    if (output.term.coef() == 0.) return output;
    // DMTK stores rightmost factors first; operator*= itself only appends.
    for (auto factor = term.factors.rbegin(); factor != term.factors.rend(); ++factor) {
        auto const id = site.intern(factor->op);
        auto op = site.operators[id];
        op.set_site(factor->site);
        output.term *= op;
        output.word.emplace_back(factor->site, id);
    }
    if (term.factors.empty()) {
        auto const id = site.intern(site.source->op("I"));
        output.term *= site.operators[id];
        output.word.emplace_back(0, id);
    }
    return output;
}

models::Term expression(models::SiteBasisPtr const& site,
                        std::initializer_list<std::pair<std::size_t, char const*>> factors,
                        double coefficient = 1.) {
    models::Term term;
    term.coeff = coefficient;
    for (auto const& [position, name] : factors)
        term.factors.push_back({static_cast<int>(position), site->op(name), site});
    return term;
}

struct Observer {
    MolecularResult& result;
    IterationCallback const& callback;
    static bool observe(dmtk::System<double>& system, std::size_t signal, void* data) {
        if (signal != dmtk::SYSTEM_SIGNAL_END_ITER) return false;
        auto& observer = *static_cast<Observer*>(data);
        Iteration row{system.get_dir() == dmtk::RIGHT2LEFT ? Direction::right_to_left : Direction::left_to_right,
                      static_cast<std::size_t>(system.get_iter()), system.energy[0],
                      system.truncation_error(), system.entropy()};
        observer.result.iterations.push_back(row);
        if (observer.callback) observer.callback(row);
        return false;
    }
};

} // namespace molecular_detail

MolecularResult solve_molecular(MolecularInput const& input, Options const& options, IterationCallback callback) {
    using namespace molecular_detail;
    // Copy controls before callbacks can mutate the caller's instance.
    const Options controls = options;
    validate(input, controls);
    auto const n = input.orbitals, n2 = n * n, n4 = n2 * n2;
    auto source = models::spinful_fermion();
    std::vector<models::Term> terms;
    std::array<char const*, 2> create = {"c_dag_up", "c_dag_dn"};
    std::array<char const*, 2> destroy = {"c_up", "c_dn"};
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < n; ++q) {
            double h = input.one_body[p * n + q];
            if (h != 0. && std::abs(h) >= input.integral_cutoff)
                for (int a = 0; a < 2; ++a)
                    terms.push_back(expression(source, {{p, create[a]}, {q, destroy[a]}}, h));
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t s = 0; s < n; ++s) {
                    double v = input.two_body[((p * n + q) * n + r) * n + s];
                    if (v == 0. || std::abs(v) < input.integral_cutoff) continue;
                    for (int a = 0; a < 2; ++a)
                        for (int b = 0; b < 2; ++b)
                            terms.push_back(expression(source, {{p, create[a]}, {r, create[b]},
                                {s, destroy[b]}, {q, destroy[a]}}, .5 * v));
                }
        }
    auto canonical = models::TermSum(std::move(terms)).canonicalise(true, 0.);
    models::validate_hermiticity(canonical, std::vector<models::SiteBasisPtr>(n, source));
    // Input tensors are no longer used after construction; a callback cannot
    // change the running Hamiltonian, charge target or scalar energy offset.
    auto const alpha = input.alpha_electrons, beta = input.beta_electrons;
    auto const core = input.core_energy;
    dmtk::InvocationGuard invocation(controls.scratch_parent, controls.log);
    dmtk::QN::init();
    auto const n_index = dmtk::QN::add_qn_index("N", true);
    auto const sz_index = dmtk::QN::add_qn_index("Sz", false);
    int const mask = (1 << n_index) | (1 << sz_index);
    dmtk::QN::set_qn_mask(mask);
    molecular_detail::Site orbitals(source, n_index, sz_index);
    models::SiteBasis neutral;
    neutral.dim = 1;
    neutral.quantum_numbers = {{"N", {0.}}, {"Sz", {0.}}};
    molecular_detail::Site padding(models::finalize_site(std::move(neutral)), n_index, sz_index);
    padding.intern(padding.source->op("I"));
    std::size_t const length = std::max<std::size_t>(n, 4);
    dmtk::Lattice lattice(static_cast<int>(length), dmtk::OBC);
    dmtk::Hami<double> h(lattice);
    for (std::size_t i = 0; i < length; ++i)
        h.sites(i) = i < n ? &orbitals.block : &padding.block;
    for (auto const& term : canonical.terms()) {
        auto native = lower(term, orbitals);
        if (native.term.coef() != 0.) h += native.term;
    }
    // Even a zero Hamiltonian needs a concrete local operator and term.
    if (h.empty()) {
        auto zero = lower(models::Term{}, orbitals).term;
        zero.coef() = 0.;
        h += zero;
    }

    struct Contribution { std::size_t slot, destination; double weight; };
    std::vector<Contribution> contributions;
    std::map<Word, std::size_t> slots;
    dmtk::Hami<double> measurements(lattice);
    auto request = [&](models::Term const& term, std::size_t destination) {
        auto native = lower(term, orbitals);
        auto const weight = native.term.coef();
        if (weight == 0.) return;
        auto [entry, inserted] = slots.emplace(std::move(native.word), measurements.size());
        if (inserted) {
            native.term.coef() = 1.;
            native.term.set_value(std::numeric_limits<double>::quiet_NaN());
            measurements += native.term;
        }
        contributions.push_back({entry->second, destination, weight});
    };
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < n; ++q) {
            for (int a = 0; a < 2; ++a)
                request(expression(source, {{q, create[a]}, {p, destroy[a]}}), static_cast<std::size_t>(a) * n2 + p * n + q);
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t s = 0; s < n; ++s)
                    for (int a = 0; a < 2; ++a)
                        for (int b = 0; b < 2; ++b)
                            request(expression(source, {{p, create[a]}, {r, create[b]},
                                {s, destroy[b]}, {q, destroy[a]}}), 2 * n2 + ((p * n + q) * n + r) * n + s);
        }
    MolecularResult result;
    result.orbitals = n;
    Observer observer{result, callback};
    dmtk::System<double> system(h, lattice, "ALPS_molecular");
    system.set_store_products(true).set_use_hc(false).set_grow_symmetric(false).set_verbose(0);
    system.signal_connect(Observer::observe, dmtk::SYSTEM_SIGNAL_END_ITER, &observer);
    system.qnt[static_cast<int>(n_index)] = static_cast<double>(alpha + beta);
    system.qnt[static_cast<int>(sz_index)] = .5 * (static_cast<double>(alpha) - static_cast<double>(beta));
    system.set_qn_mask(mask);
    if (controls.lanczos_tolerance) system.set_lanczos_tolerance(*controls.lanczos_tolerance);
    if (controls.truncation_error)
        system.set_error(*controls.truncation_error, controls.truncation_max_states
            ? static_cast<int>(*controls.truncation_max_states) : -1);
    dmtk::Matrix<std::size_t> schedule(2, controls.sweeps.size());
    for (std::size_t i = 0; i < controls.sweeps.size(); ++i) {
        schedule(0, i) = controls.sweeps[i].right_to_left;
        schedule(1, i) = controls.sweeps[i].left_to_right;
    }
    system.corr = measurements;
    system.start(controls.sweeps.size(), schedule);
    system.run(static_cast<int>(controls.warmup_states));
    system.set_full_sweep(true);
    system.final_sweep(controls.sweeps.back().left_to_right, dmtk::RIGHT2LEFT, 1, true);
    // Read every expectation from the same final state and retained products.
    // A missing measurement is an error, never an implicit zero RDM element.
    for (auto& term : system.corr) term.set_value(std::numeric_limits<double>::quiet_NaN());
    system.measure();
    std::vector<double> values(2 * n2 + n4, 0.);
    for (auto const& item : contributions) {
        auto const measured = system.corr[item.slot].value();
        if (!std::isfinite(measured))
            throw std::runtime_error("classic molecular engine did not evaluate a required RDM element");
        values[item.destination] += item.weight * measured;
    }
    result.energy = system.energy[0] + core;
    if (!std::isfinite(result.energy))
        throw std::runtime_error("classic molecular engine returned a nonfinite energy");
    result.final_truncation_error = system.truncation_error();
    result.one_rdm_alpha.assign(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(n2));
    result.one_rdm_beta.assign(values.begin() + static_cast<std::ptrdiff_t>(n2), values.begin() + static_cast<std::ptrdiff_t>(2 * n2));
    result.two_rdm.assign(values.begin() + static_cast<std::ptrdiff_t>(2 * n2), values.end());
    return result;
}

} // namespace alps::dmrg::classic
