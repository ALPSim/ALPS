// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Written for this tree: adaptive time-dependent DMRG on the dmtk engine,
// the method of S. R. White and A. E. Feiguin, Phys. Rev. Lett. 93, 076401
// (2004). The sweep schedule (even bonds dt/2, odd bonds dt, even bonds
// dt/2, a measuring sweep with identity gates; gates on the central bond
// and on an edge bond whose block is a single site; no Lanczos solve) was
// read from Adrian Feiguin's teaching program afeiguin/dmrg-tutorials,
// 06-heischain_time/heischain_time.cc, a dense spin-1/2 Heisenberg code
// without quantum numbers. That repository carries no license, so it served
// as a description of the algorithm only: no code was copied or translated
// from it. The code below is new and works on dmtk's charge-sorted states,
// its wave-function transformation (new_seed) and its truncation.
//
// One step of a half sweep, at the configuration (left block, site s, site
// s + 1, right block):
//  1. the closing state of the previous step is carried to this
//     configuration by White's wave-function transformation (new_seed,
//     which folds the grown block's rotation in and unfolds the shrinking
//     block's one; exact up to the previous truncation);
//  2. the gates of this half sweep act on the bonds whose two sites are
//     explicit: the central bond (s, s + 1), and (0, 1) or (L-2, L-1) while
//     the left or the right block is a single site;
//  3. the state is renormalized and the density matrix of the evolved state
//     selects the kept states of the grown block (dmtk's truncate and
//     rotate, with the evolved state as the only target).
// A half sweep runs over the whole chain (left blocks of 1 to L-3 sites), so
// every bond is visited once in each half sweep, and it ends at the
// configuration the next half sweep starts from: the next one continues
// from the same state without a transformation. The Hamiltonian is not
// applied, so its block operators are neither built nor rotated while the
// system evolves; the operators measured (dmtk's corr) still are.
//
// Time-step targeting, the method of A. E. Feiguin and S. R. White, Phys.
// Rev. B 72, 020404(R) (2005), written for this tree from the paper: the
// final Hamiltonian replaces the System's Hamiltonian, and its block
// operators are built by the same hooks the ground-state sweeps use. Its
// block Hamiltonians carry their own name, so a block left from the
// preparation holds none. The first half sweep rebuilds the blocks on one
// side and the second those on the other; from then on dmtk's superblock
// product is the final Hamiltonian. At every configuration the state psi is
// carried there as above, and four Runge-Kutta vectors of
// k = c H (c = -i dt in real time, -dt in imaginary time) give the
// approximations
//   psi(t + dt/3)  = psi + (31 k1 + 14 k2 + 14 k3 - 5 k4) / 162,
//   psi(t + 2dt/3) = psi + (16 k1 + 20 k2 + 20 k3 - 2 k4) / 81,
//   psi(t + dt)    = psi + (k1 + 2 k2 + 2 k3 + k4) / 6
// of the paper. The density matrix targets psi and these three with weights
// 1/3, 1/6, 1/6 and 1/3. At the last configuration of a step, psi becomes
// exp(c H) psi by a Lanczos (Krylov) exponential, which replaces the
// Runge-Kutta psi(t + dt) among the targets.
//
// This private header is included only by classic.cpp.
#pragma once

#include "dmtk/dmtk.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace alps::dmrg::classic::detail {

// A two-site operator on sites (i, i + 1) in the engine's charge-sorted
// local bases: `dense` is [out][in], row-major, with the pair index
// p * second + q for state p of the first site and q of the second. An empty
// gate (first == 0) is the identity.
template<class T>
struct BondGate {
    std::size_t first = 0, second = 0;
    std::vector<T> dense;
    // The nonzero entries of every input column: (out, value).
    std::vector<std::vector<std::pair<std::size_t, T>>> columns;

    bool identity() const { return first == 0; }

    void index() {
        std::size_t const n = first * second;
        columns.assign(n, {});
        for (std::size_t out = 0; out < n; ++out)
            for (std::size_t in = 0; in < n; ++in)
                if (dense[out * n + in] != T(0)) columns[in].emplace_back(out, dense[out * n + in]);
    }
};

// later * earlier: the gate that applies `earlier` first. Either may be the identity.
template<class T>
BondGate<T> compose(BondGate<T> const& later, BondGate<T> const& earlier) {
    if (earlier.identity()) return later;
    if (later.identity()) return earlier;
    if (later.first != earlier.first || later.second != earlier.second)
        throw std::logic_error("classic DMRG time evolution composes gates of different bonds");
    std::size_t const n = later.first * later.second;
    BondGate<T> out{later.first, later.second, std::vector<T>(n * n, T(0)), {}};
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t k = 0; k < n; ++k) {
            T const a = later.dense[r * n + k];
            if (a == T(0)) continue;
            for (std::size_t c = 0; c < n; ++c) out.dense[r * n + c] += a * earlier.dense[k * n + c];
        }
    out.index();
    return out;
}

// exp(-i t h) in real time, exp(-t h) in imaginary time, of a Hermitian
// two-site matrix h ([out][in], dimension first * second). The exponential
// is taken on each connected block of h's nonzero pattern, so states h does
// not connect (other charges or parities) stay exactly unconnected.
template<class T>
BondGate<T> exponential(std::vector<std::complex<double>> const& h, std::size_t first, std::size_t second,
                        double t, bool imaginary) {
    using Complex = std::complex<double>;
    std::size_t const n = first * second;
    if (h.size() != n * n) throw std::logic_error("classic DMRG bond Hamiltonian has the wrong dimension");
    std::vector<std::size_t> parent(n);
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    auto root = [&](std::size_t x) {
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    };
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c)
            if (h[r * n + c] != Complex{}) parent[root(r)] = root(c);
    std::map<std::size_t, std::vector<std::size_t>> blocks;
    for (std::size_t r = 0; r < n; ++r) blocks[root(r)].push_back(r);
    BondGate<T> gate{first, second, std::vector<T>(n * n, T(0)), {}};
    for (auto const& [key, members] : blocks) {
        auto const k = static_cast<Eigen::Index>(members.size());
        Eigen::MatrixXcd block(k, k);
        for (Eigen::Index i = 0; i < k; ++i)
            for (Eigen::Index j = 0; j < k; ++j) {
                // The Hermitian part; the driver checked Hermiticity on the terms.
                Complex const a = h[members[i] * n + members[j]], b = h[members[j] * n + members[i]];
                block(i, j) = 0.5 * (a + std::conj(b));
            }
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(block);
        if (solver.info() != Eigen::Success)
            throw std::runtime_error("classic DMRG could not diagonalize a bond Hamiltonian");
        Eigen::VectorXcd phases(k);
        for (Eigen::Index i = 0; i < k; ++i) {
            double const e = solver.eigenvalues()[i];
            phases[i] = imaginary ? Complex(std::exp(-t * e)) : std::exp(Complex(0., -t * e));
        }
        Eigen::MatrixXcd const u =
            solver.eigenvectors() * phases.asDiagonal() * solver.eigenvectors().adjoint();
        for (Eigen::Index i = 0; i < k; ++i)
            for (Eigen::Index j = 0; j < k; ++j) {
                Complex const v = u(i, j);
                T& entry = gate.dense[members[i] * n + members[j]];
                if constexpr (std::is_same_v<T, double>) {
                    if (std::abs(v.imag()) > 1e-12 * std::max(1., std::abs(v)))
                        throw std::logic_error("classic DMRG real arithmetic produced a complex gate");
                    entry = v.real();
                } else {
                    entry = v;
                }
            }
    }
    gate.index();
    return gate;
}

// Apply a gate to the dmtk blocks `first` and `first + 1` (numbered 1 to 4)
// of a state, both of which must be single sites.
template<class T>
void apply_gate(dmtk::VectorState<T>& state, int first, BondGate<T> const& gate) {
    if (gate.identity()) return;
    std::array<dmtk::PackedBasis const*, 4> const blocks{&state.b1(), &state.b2(), &state.b3(), &state.b4()};
    int const a = first - 1, b = first; // 0-based positions of the two sites
    if (blocks[a]->dim() != gate.first || blocks[b]->dim() != gate.second)
        throw std::logic_error("classic DMRG gate does not match the sites it acts on");
    // The subspace of each site state.
    auto subspace_of = [](dmtk::PackedBasis const& basis) {
        std::vector<std::size_t> out(basis.dim());
        for (auto const& space : basis)
            for (std::size_t s = space.begin(); s <= static_cast<std::size_t>(space.end()); ++s) out[s] = space.begin();
        return out;
    };
    std::vector<std::size_t> const sub_a = subspace_of(*blocks[a]), sub_b = subspace_of(*blocks[b]);
    // The state's subspaces by the first states of their four block subspaces.
    std::vector<dmtk::StateSpace const*> spaces;
    std::map<std::array<std::size_t, 4>, std::size_t> index;
    for (auto it = state.subspace_begin(); it != state.subspace_end(); ++it) {
        dmtk::StateSpace const& space = *it;
        index.emplace(std::array<std::size_t, 4>{space[1].begin(), space[2].begin(), space[3].begin(), space[4].begin()},
                      spaces.size());
        spaces.push_back(&space);
    }
    // The two block positions the gate leaves alone, the outer one first.
    std::array<int, 2> spectators{};
    for (int k = 0, j = 0; k < 4; ++k)
        if (k != a && k != b) spectators[static_cast<std::size_t>(j++)] = k;
    auto strides = [](dmtk::StateSpace const& space) {
        std::array<std::size_t, 4> out{};
        std::size_t stride = 1;
        for (int k = 3; k >= 0; --k) {
            out[static_cast<std::size_t>(k)] = stride;
            stride *= space[static_cast<std::size_t>(k + 1)].dim();
        }
        return out;
    };
    double largest = 0.;
    for (auto const& column : gate.columns)
        for (auto const& entry : column) largest = std::max(largest, std::abs(entry.second));
    dmtk::VectorState<T> result(state);
    result = T(0);
    for (dmtk::StateSpace const* source : spaces) {
        dmtk::StateSpace const& from = *source;
        auto const from_strides = strides(from);
        std::array<std::size_t, 4> key{from[1].begin(), from[2].begin(), from[3].begin(), from[4].begin()};
        std::size_t const outer = from[static_cast<std::size_t>(spectators[0] + 1)].dim();
        std::size_t const inner = from[static_cast<std::size_t>(spectators[1] + 1)].dim();
        for (std::size_t p = 0; p < from[static_cast<std::size_t>(a + 1)].dim(); ++p)
            for (std::size_t q = 0; q < from[static_cast<std::size_t>(b + 1)].dim(); ++q) {
                std::size_t const in = (static_cast<std::size_t>(from[static_cast<std::size_t>(a + 1)].begin()) + p)
                                     * gate.second
                                     + static_cast<std::size_t>(from[static_cast<std::size_t>(b + 1)].begin()) + q;
                std::size_t const source_base =
                    static_cast<std::size_t>(from.start()) + p * from_strides[static_cast<std::size_t>(a)]
                    + q * from_strides[static_cast<std::size_t>(b)];
                for (auto const& [out, value] : gate.columns[in]) {
                    std::size_t const p_out = out / gate.second, q_out = out % gate.second;
                    key[static_cast<std::size_t>(a)] = sub_a[p_out];
                    key[static_cast<std::size_t>(b)] = sub_b[q_out];
                    auto found = index.find(key);
                    if (found == index.end()) {
                        // A charge the sector fixes would change: only rounding may do that.
                        if (std::abs(value) <= 1e-13 * largest) continue;
                        throw std::logic_error("classic DMRG gate leaves the conserved sector");
                    }
                    dmtk::StateSpace const& to = *spaces[found->second];
                    auto const to_strides = strides(to);
                    std::size_t const target_base =
                        static_cast<std::size_t>(to.start())
                        + (p_out - static_cast<std::size_t>(to[static_cast<std::size_t>(a + 1)].begin()))
                              * to_strides[static_cast<std::size_t>(a)]
                        + (q_out - static_cast<std::size_t>(to[static_cast<std::size_t>(b + 1)].begin()))
                              * to_strides[static_cast<std::size_t>(b)];
                    std::size_t const so = from_strides[static_cast<std::size_t>(spectators[0])],
                                      si = from_strides[static_cast<std::size_t>(spectators[1])],
                                      to_o = to_strides[static_cast<std::size_t>(spectators[0])],
                                      to_i = to_strides[static_cast<std::size_t>(spectators[1])];
                    for (std::size_t x = 0; x < outer; ++x)
                        for (std::size_t y = 0; y < inner; ++y)
                            result[target_base + x * to_o + y * to_i] += value * state[source_base + x * so + y * si];
                    key[static_cast<std::size_t>(a)] = from[static_cast<std::size_t>(a + 1)].begin();
                    key[static_cast<std::size_t>(b)] = from[static_cast<std::size_t>(b + 1)].begin();
                }
            }
    }
    state = result;
}

// The dmtk system of a time evolution: the ground state of its Hamiltonian
// by the usual warm-up and sweeps, then half sweeps that apply bond gates.
template<class T>
class EvolvingSystem : public dmtk::System<T> {
    using Base = dmtk::System<T>;

public:
    EvolvingSystem(dmtk::Hami<T> const& h, dmtk::Lattice const& lattice) : Base(h, lattice, "ALPS") {}

    // The warm-up and the scheduled sweeps of System::run, the Lanczos
    // residuals of the last sweep certified down to the tolerance (as
    // final_sweep certifies the states it measures). A four-site chain fits
    // no sweep step: its single configuration is solved directly. With
    // `sectors`, the warm-up's k-site chain targets sectors[k].
    void prepare(int warmup, std::optional<dmtk::Vector<dmtk::QN>> const& sectors = std::nullopt) {
        int const ls = this->lattice().ls();
        if (warmup > 0 && sectors) {
            this->_custom_qns = true;
            this->warmup_loop(static_cast<std::size_t>(warmup), *sectors);
        } else if (warmup > 0) {
            this->warmup_loop(static_cast<std::size_t>(warmup));
        }
        for (std::size_t s = 1; s <= this->_numsweeps; ++s) {
            this->_lanczos_certify = s == this->_numsweeps;
            dmtk::out() << "NEW SWEEP " << s << " " << this->_sweeps(0, s - 1) << " " << this->_sweeps(1, s - 1) << '\n';
            this->sweep(this->_sweeps(0, s - 1), this->_sweeps(1, s - 1), dmtk::RIGHT2LEFT, 1);
        }
        continued_ = false;
        if (ls == 4) {
            this->_lanczos_certify = true;
            this->_in_warmup = false;
            this->dir = dmtk::RIGHT2LEFT;
            this->iter = 1;
            this->qn = this->qnt;
            this->read_block(this->leftblock, 1, dmtk::LEFT);
            this->read_block(this->rightblock, 1, dmtk::RIGHT);
            this->init_iteration(this->leftblock, this->h.get_site(1), this->h.get_site(2), this->rightblock, false);
            this->diagonalize(false);
            continued_ = true; // the first half sweep starts at this configuration
        }
        this->_lanczos_certify = false;
    }

    // From here on no Lanczos solve runs: the Hamiltonian's block operators
    // are no longer built, and truncations keep `states` states or, with a
    // truncation error, as many as it needs up to `states`.
    void begin_evolution(std::size_t states, std::optional<double> truncation_error) {
        evolving_ = true;
        states_ = static_cast<int>(states);
        this->set_rotate_terms(false); // no local-energy diagnostic of the old Hamiltonian
        this->set_store_products(false);
        this->_use_error = truncation_error.has_value();
        if (truncation_error) {
            this->_error = *truncation_error;
            this->_error_max_size = states_;
        }
    }

    bool evolving() const { return evolving_; }

    // One half sweep over the whole chain in `direction`, applying gates[b]
    // to bond (b, b + 1) at the step where both its sites are explicit.
    // `measure` measures corr at every step, as the unstored path of sweep()
    // does, so every term is measured at some configuration. Returns the
    // discarded weight summed over the half sweep's truncations.
    double half_sweep(std::size_t direction, std::vector<BondGate<T> const*> const& gates, bool measure) {
        int const ls = this->lattice().ls();
        if (gates.size() != static_cast<std::size_t>(ls - 1))
            throw std::logic_error("classic DMRG time evolution needs one gate per bond");
        visits_.assign(gates.size(), 0);
        this->_in_warmup = false;
        this->_custom_qns = false;
        this->qn = this->qnt;
        this->dir = direction;
        bool const right_to_left = direction == dmtk::RIGHT2LEFT;
        double discarded = 0.;
        for (int step = 1; step < ls - 2; ++step) {
            this->iter = step;
            this->signal_emit(dmtk::SYSTEM_SIGNAL_START_ITER);
            int const left_sites = right_to_left ? ls - step - 2 : step;
            int const right_sites = ls - left_sites - 2;
            this->read_block(this->leftblock, left_sites, dmtk::LEFT);
            this->read_block(this->rightblock, right_sites, dmtk::RIGHT);
            dmtk::Block<T> const& site1 = this->h.get_site(left_sites);
            dmtk::Block<T> const& site2 = this->h.get_site(left_sites + 1);
            this->m1 = static_cast<int>(this->leftblock.dim());
            this->m4 = static_cast<int>(this->rightblock.dim());
            this->m2 = static_cast<int>(site1.dim());
            this->m3 = static_cast<int>(site2.dim());
            if (step == 1 && continued_) {
                // The configuration the previous half sweep closed on.
                dmtk::VectorState<T> const previous = this->gs;
                this->init_iteration(this->leftblock, site1, site2, this->rightblock, false);
                this->gs = previous;
            } else {
                this->init_iteration(this->leftblock, site1, site2, this->rightblock, true);
                this->gs = this->seed;
            }
            auto apply = [&](int block, int bond) {
                ++visits_[static_cast<std::size_t>(bond)];
                if (gates[static_cast<std::size_t>(bond)]) apply_gate(this->gs, block, *gates[static_cast<std::size_t>(bond)]);
            };
            if (left_sites == 1) apply(1, 0);
            apply(2, left_sites);
            if (right_sites == 1) apply(3, ls - 2);
            double const norm = std::sqrt(std::real(product(this->gs, this->gs)));
            if (!(norm > 0.) || !std::isfinite(norm))
                throw std::runtime_error("classic DMRG time evolution: the state vanished");
            this->gs /= T(norm);
            this->energy[0] = std::numeric_limits<double>::quiet_NaN();
            this->m = states_;
            if (right_to_left) {
                this->m4 = std::min(this->m4 * this->m3, this->m);
                this->truncate(dmtk::RIGHT, this->m4);
                this->rotate(dmtk::RIGHT, this->newblock);
                this->write_iter(dmtk::RIGHT);
            } else {
                this->m1 = std::min(this->m1 * this->m2, this->m);
                this->truncate(dmtk::LEFT, this->m1);
                this->rotate(dmtk::LEFT, this->newblock);
                this->write_iter(dmtk::LEFT);
            }
            discarded += std::max(0., this->truncation_error());
            if (measure) this->measure();
            this->signal_emit(dmtk::SYSTEM_SIGNAL_END_ITER);
        }
        for (int count : visits_)
            if (count != 1) throw std::logic_error("classic DMRG time evolution visited a bond other than once");
        continued_ = true;
        return discarded;
    }

    // Time-step targeting from here on: `final_hamiltonian` (named apart
    // from the initial one) becomes the System's Hamiltonian, and the hooks
    // build and rotate its block operators. The Trotter method's
    // begin_evolution turns off the interaction terms (rotate_terms) that
    // dmtk's superblock product needs, so they are turned on again here.
    void begin_targeting(dmtk::Hami<T> const& final_hamiltonian, double time_step, bool imaginary,
                         std::size_t states, std::optional<double> truncation_error) {
        begin_evolution(states, truncation_error);
        targeting_ = true;
        this->set_rotate_terms(true);
        this->h = final_hamiltonian;
        imaginary_ = imaginary;
        time_step_ = time_step;
        built_ = 0;
    }

    // One half sweep of time-step targeting in `direction`. Until a half
    // sweep has built one side's blocks of the final Hamiltonian, the
    // truncations target the state alone and the iterations carry no energy
    // (NaN). After that every configuration targets the four states and
    // reports <psi|H|psi> as its energy, and with `advance` the last one
    // takes the Krylov step. `measure` measures corr at every configuration.
    // Returns the squared norm the state lost in the half sweep's
    // truncations.
    double targeting_half_sweep(std::size_t direction, bool advance, bool measure) {
        if (!targeting_) throw std::logic_error("classic DMRG time-step targeting has not begun");
        int const ls = this->lattice().ls();
        bool const built = built_ > 0;
        if (advance && !built)
            throw std::logic_error("classic DMRG time-step targeting advanced before its Hamiltonian was built");
        this->_in_warmup = false;
        this->_custom_qns = false;
        this->qn = this->qnt;
        this->dir = direction;
        bool const right_to_left = direction == dmtk::RIGHT2LEFT;
        double lost = 0.;
        for (int step = 1; step < ls - 2; ++step) {
            this->iter = step;
            this->signal_emit(dmtk::SYSTEM_SIGNAL_START_ITER);
            int const left_sites = right_to_left ? ls - step - 2 : step;
            int const right_sites = ls - left_sites - 2;
            this->read_block(this->leftblock, left_sites, dmtk::LEFT);
            this->read_block(this->rightblock, right_sites, dmtk::RIGHT);
            dmtk::Block<T> const& site1 = this->h.get_site(left_sites);
            dmtk::Block<T> const& site2 = this->h.get_site(left_sites + 1);
            this->m1 = static_cast<int>(this->leftblock.dim());
            this->m4 = static_cast<int>(this->rightblock.dim());
            this->m2 = static_cast<int>(site1.dim());
            this->m3 = static_cast<int>(site2.dim());
            if (step == 1 && continued_) {
                dmtk::VectorState<T> const previous = this->gs;
                this->init_iteration(this->leftblock, site1, site2, this->rightblock, false);
                this->gs = previous;
            } else {
                this->init_iteration(this->leftblock, site1, site2, this->rightblock, true);
                this->gs = this->seed;
                // The weight of the state the previous truncation kept.
                lost += std::max(0., 1. - std::real(product(this->gs, this->gs)));
            }
            double const norm = std::sqrt(std::real(product(this->gs, this->gs)));
            if (!(norm > 0.) || !std::isfinite(norm))
                throw std::runtime_error("classic DMRG time evolution: the state vanished");
            this->gs /= T(norm);
            use_targets_ = built;
            this->energy[0] = std::numeric_limits<double>::quiet_NaN();
            if (built) {
                rk_targets();
                if (advance && step == ls - 3) {
                    this->gs = krylov_step(this->gs);
                    double const evolved = std::sqrt(std::real(product(this->gs, this->gs)));
                    if (!(evolved > 0.) || !std::isfinite(evolved))
                        throw std::runtime_error("classic DMRG time evolution: the state vanished");
                    this->gs /= T(evolved);
                    targets_[3] = this->gs;
                }
            }
            this->m = states_;
            if (right_to_left) {
                this->m4 = std::min(this->m4 * this->m3, this->m);
                this->truncate(dmtk::RIGHT, this->m4);
                this->rotate(dmtk::RIGHT, this->newblock);
                this->write_iter(dmtk::RIGHT);
            } else {
                this->m1 = std::min(this->m1 * this->m2, this->m);
                this->truncate(dmtk::LEFT, this->m1);
                this->rotate(dmtk::LEFT, this->newblock);
                this->write_iter(dmtk::LEFT);
            }
            use_targets_ = false;
            if (measure) this->measure();
            this->signal_emit(dmtk::SYSTEM_SIGNAL_END_ITER);
        }
        ++built_;
        continued_ = true;
        return lost;
    }

    // The truncation's targets: the state alone, or during time-step
    // targeting the four states of rk_targets.
    void build_target_states(int position) override {
        if (!use_targets_) {
            this->_target.resize(1);
            this->_target_weight.resize(1);
            Base::build_target_states(position);
            return;
        }
        static constexpr double weights[4] = {1. / 3., 1. / 6., 1. / 6., 1. / 3.};
        this->_target.resize(4);
        this->_target_weight.resize(4);
        for (std::size_t k = 0; k < 4; ++k) {
            this->_target[k] = targets_[k];
            this->_target_weight[k] = weights[k];
        }
    }

    // The Lanczos product's terms and the block operators of the Hamiltonian
    // are built before the evolution begins, and again during time-step
    // targeting (for the final Hamiltonian).
    void init_hamiltonian(bool create_composite) override {
        if (!evolving_ || targeting_) Base::init_hamiltonian(create_composite);
    }
    void rotate_hamiltonian(int position, dmtk::Block<T>& b, dmtk::Basis& basis, dmtk::Basis& rho_basis) override {
        if (!evolving_ || targeting_) Base::rotate_hamiltonian(position, b, basis, rho_basis);
    }

private:
    // c = -i dt in real time, -dt in imaginary time (real arithmetic then).
    T step_factor() const {
        if constexpr (std::is_same_v<T, double>) {
            return -time_step_;
        } else {
            return imaginary_ ? T(-time_step_) : T(0., -time_step_);
        }
    }

    dmtk::VectorState<T> apply(dmtk::VectorState<T> const& v) { return dmtk::product(static_cast<Base&>(*this), v); }

    // psi and its Runge-Kutta states at t + dt/3, t + 2dt/3 and t + dt (the
    // formulas in the header comment); the energy <psi|H|psi> on the way.
    void rk_targets() {
        auto const& psi = this->gs;
        T const c = step_factor();
        dmtk::VectorState<T> const h_psi = apply(psi);
        this->energy[0] = std::real(product(psi, h_psi));
        dmtk::VectorState<T> const k1 = c * h_psi;
        dmtk::VectorState<T> stage = psi;
        stage += T(0.5) * k1;
        dmtk::VectorState<T> const k2 = c * apply(stage);
        stage = psi;
        stage += T(0.5) * k2;
        dmtk::VectorState<T> const k3 = c * apply(stage);
        stage = psi;
        stage += k3;
        dmtk::VectorState<T> const k4 = c * apply(stage);
        auto combine = [&](double a1, double a2, double a3, double a4) {
            dmtk::VectorState<T> out = psi;
            out += T(a1) * k1;
            out += T(a2) * k2;
            out += T(a3) * k3;
            out += T(a4) * k4;
            return out;
        };
        targets_[0] = psi;
        targets_[1] = combine(31. / 162., 14. / 162., 14. / 162., -5. / 162.);
        targets_[2] = combine(16. / 81., 20. / 81., 20. / 81., -2. / 81.);
        targets_[3] = combine(1. / 6., 2. / 6., 2. / 6., 1. / 6.);
    }

    // exp(c H) psi in the Krylov space of H and psi (Lanczos with full
    // reorthogonalization), grown until the coefficient of its newest vector
    // is below 1e-13 or the space closes.
    dmtk::VectorState<T> krylov_step(dmtk::VectorState<T> const& psi) {
        constexpr std::size_t max_dimension = 40;
        T const c = step_factor();
        double const norm = std::sqrt(std::real(product(psi, psi)));
        std::vector<dmtk::VectorState<T>> basis;
        basis.push_back(psi);
        basis.back() /= T(norm);
        std::vector<double> alpha, beta;
        std::vector<T> coefficients;
        for (std::size_t j = 0;; ++j) {
            dmtk::VectorState<T> w = apply(basis[j]);
            alpha.push_back(std::real(product(basis[j], w)));
            for (int pass = 0; pass < 2; ++pass)
                for (auto const& v : basis) w -= product(v, w) * v;
            double const b = std::sqrt(std::real(product(w, w)));
            std::size_t const n = alpha.size();
            Eigen::MatrixXd tridiagonal = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(n));
            for (std::size_t i = 0; i < n; ++i) {
                tridiagonal(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i)) = alpha[i];
                if (i + 1 < n)
                    tridiagonal(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i + 1)) =
                        tridiagonal(static_cast<Eigen::Index>(i + 1), static_cast<Eigen::Index>(i)) = beta[i];
            }
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> const solver(tridiagonal);
            coefficients.assign(n, T(0));
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t e = 0; e < n; ++e) {
                    auto const ei = static_cast<Eigen::Index>(e);
                    coefficients[i] += solver.eigenvectors()(static_cast<Eigen::Index>(i), ei)
                                       * std::exp(c * solver.eigenvalues()(ei)) * solver.eigenvectors()(0, ei);
                }
            bool const closed = b <= 1e-12 * std::max(1., std::abs(alpha.back()));
            if (closed || (n > 1 && std::abs(coefficients.back()) < 1e-13)) break;
            if (n == max_dimension)
                throw std::runtime_error("classic DMRG time-step targeting: the Krylov step did not converge in "
                                         + std::to_string(max_dimension) + " vectors; use a smaller time_step");
            beta.push_back(b);
            w /= T(b);
            basis.push_back(std::move(w));
        }
        dmtk::VectorState<T> out = basis[0];
        out *= T(norm) * coefficients[0];
        for (std::size_t i = 1; i < coefficients.size(); ++i) out += (T(norm) * coefficients[i]) * basis[i];
        return out;
    }

    bool evolving_ = false;
    bool targeting_ = false;
    bool imaginary_ = false;
    bool use_targets_ = false; // the next truncation targets targets_
    double time_step_ = 0.;
    int built_ = 0;            // targeting half sweeps since begin_targeting
    bool continued_ = false; // the state is at the configuration a half sweep starts from
    int states_ = 0;
    std::vector<int> visits_;
    std::array<dmtk::VectorState<T>, 4> targets_;
};

} // namespace alps::dmrg::classic::detail
