// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Adapted from ALPSim/ALPS master, PR 216 by Marcus Rosales (commits
// 46c76b1d2 and 35cea33e7): DMRGTask::extract_mps in
// src/apps/dmrg/dmrg/dmrg.h and split_two_site in
// src/apps/dmrg/dmrg/mps_export.h, which assembled the dmrg application's
// ground state from the same engine's stored block rotations. Changes: the
// closing configuration is read from the System in memory after
// final_sweep (block sizes from its left and right blocks, every targeted
// state, the wave functions the measurements were taken on) instead of the
// ground-state file of the step at L/2; the SVD is Eigen's BDCSVD instead
// of LAPACK gesvd; the tensors are returned as classic::MatrixProductState.
// The row labels of the rotations (left: old block state i1, site state i2;
// right: site state i1, old block state i2) and the subspace walk of the
// wave function are #216's; the charge-sorted site bases map to
// PhysicalBasis indices through State::i1, as in molecular.ipp's CiWalk.
//
// This private header is included only by classic.cpp; the System's
// scratch files must still exist when an MpsExport is built.
#pragma once

#include "dmtk/dmtk.hpp"

#include <alps/dmrg/classic.hpp>

#include <Eigen/SVD>

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace alps::dmrg::classic::detail {

template<class T>
class MpsExport {
public:
    // Reads the block rotations of the final sweep up to its closing
    // configuration: left block of n sites, sites n and n+1, right block.
    explicit MpsExport(dmtk::System<T>& system) : system_(system) {
        int const L = system.lattice().ls();
        n_ = system.leftblock.lattice().ls();
        int const r = system.rightblock.lattice().ls();
        if (L < 4 || n_ < 1 || r < 1 || n_ + r + 2 != L)
            throw std::logic_error("classic DMRG MPS export: the final sweep did not close on the whole chain");
        auto const length = static_cast<std::size_t>(L);
        kit_index_.resize(length);
        for (std::size_t i = 0; i < length; ++i) {
            auto const& basis = system.h.get_site(static_cast<int>(i)).basis();
            for (std::size_t p = 0; p < basis.size(); ++p) kit_index_[i].push_back(basis[p].i1);
        }
        shapes_.resize(length);
        tensors_.resize(length);
        // The single-site blocks at either end: the bond index is the sorted site state.
        for (std::size_t i : {std::size_t(0), length - 1}) {
            std::size_t const d = kit_index_[i].size();
            shapes_[i] = i == 0 ? Shape{1, d, d} : Shape{d, d, 1};
            tensors_[i].assign(d * d, T(0));
            for (std::size_t p = 0; p < d; ++p)
                tensors_[i][i == 0 ? kit_index_[i][p] * d + p : p * d + kit_index_[i][p]] = T(1);
        }
        // Left block k adds site k-1 to block k-1; right block k adds site L-k.
        for (int k = 2; k <= n_; ++k) read_rotation(k, dmtk::LEFT, static_cast<std::size_t>(k - 1));
        for (int k = 2; k <= r; ++k) read_rotation(k, dmtk::RIGHT, static_cast<std::size_t>(L - k));
    }

    // `state` on the closing configuration (System::gs or a target).
    MatrixProductState state(dmtk::VectorState<T> const& state) const {
        auto const n = static_cast<std::size_t>(n_);
        std::size_t const da = shapes_[n - 1][2], ds = kit_index_[n].size(), dt = kit_index_[n + 1].size(),
                          db = shapes_[n + 2][0];
        if (state.b1().dim() != da || state.b2().dim() != ds || state.b3().dim() != dt || state.b4().dim() != db)
            throw std::logic_error("classic DMRG MPS export: the wave function does not match the final blocks");
        std::vector<T> psi(da * ds * dt * db, T(0));
        // dmtk stores the values by charge subspace, each a-major with b fastest.
        for (auto sp = state.subspace_begin(); sp != state.subspace_end(); ++sp) {
            dmtk::StateSpace const& space = *sp;
            std::size_t index = static_cast<std::size_t>(space.start());
            auto const first = [&](std::size_t k) { return static_cast<std::size_t>(space[k].begin()); };
            auto const last = [&](std::size_t k) { return first(k) + space[k].dim(); };
            for (std::size_t a = first(1); a < last(1); ++a)
                for (std::size_t s = first(2); s < last(2); ++s)
                    for (std::size_t t = first(3); t < last(3); ++t)
                        for (std::size_t b = first(4); b < last(4); ++b)
                            psi[((a * ds + kit_index_[n][s]) * dt + kit_index_[n + 1][t]) * db + b] = state[index++];
        }
        auto shapes = shapes_;
        auto tensors = tensors_;
        split(psi, da, ds, dt, db, tensors[n], shapes[n], tensors[n + 1], shapes[n + 1]);
        MatrixProductState out;
        out.center = n + 1;
        out.sites.resize(tensors.size());
        for (std::size_t i = 0; i < tensors.size(); ++i) {
            auto& site = out.sites[i];
            site.left = shapes[i][0];
            site.local = shapes[i][1];
            site.right = shapes[i][2];
            if constexpr (std::is_same_v<T, double>) site.real = std::move(tensors[i]);
            else site.complex = std::move(tensors[i]);
        }
        return out;
    }

private:
    using Shape = std::array<std::size_t, 3>;

    void read_rotation(int k, int position, std::size_t site) {
        bool const left = position == dmtk::LEFT;
        dmtk::BMatrix<T> rotation;
        dmtk::Basis rows;
        system_.read_rho(rotation, rows, k, position);
        std::size_t const old_dim = left ? shapes_[site - 1][2] : shapes_[site + 1][0];
        std::size_t const d = kit_index_[site].size();
        std::size_t new_dim = 0;
        for (auto it = rotation.begin(); it != rotation.end(); ++it)
            new_dim = std::max<std::size_t>(new_dim, static_cast<std::size_t>(it->col_range().end() + 1));
        if (new_dim == 0) throw std::logic_error("classic DMRG MPS export: an empty block rotation");
        shapes_[site] = left ? Shape{old_dim, d, new_dim} : Shape{new_dim, d, old_dim};
        auto& tensor = tensors_[site];
        tensor.assign(old_dim * d * new_dim, T(0));
        for (auto it = rotation.begin(); it != rotation.end(); ++it) {
            auto const r0 = static_cast<std::size_t>(it->row_range().begin());
            auto const c0 = static_cast<std::size_t>(it->col_range().begin());
            for (std::size_t i = 0; i < static_cast<std::size_t>(it->row_range().size()); ++i) {
                dmtk::State const& row = rows[r0 + i];
                std::size_t const old_state = left ? row.i1 : row.i2;
                std::size_t const s = kit_index_[site][left ? row.i2 : row.i1];
                if (old_state >= old_dim)
                    throw std::logic_error("classic DMRG MPS export: the block rotations do not chain");
                for (std::size_t j = 0; j < static_cast<std::size_t>(it->col_range().size()); ++j) {
                    std::size_t const c = c0 + j;
                    T const value = (*it)(static_cast<int>(j), static_cast<int>(i));
                    if (left) tensor[(old_state * d + s) * new_dim + c] = value;
                    else tensor[(c * d + s) * old_dim + old_state] = value;
                }
            }
        }
    }

    // psi[a][s][t][b] = sum_x left[a][s][x] right[x][t][b], left orthonormal
    // and right = S V^H.
    static void split(std::vector<T> const& psi, std::size_t da, std::size_t ds, std::size_t dt, std::size_t db,
                      std::vector<T>& left, Shape& left_shape, std::vector<T>& right, Shape& right_shape) {
        using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        Eigen::Index const rows = static_cast<Eigen::Index>(da * ds), cols = static_cast<Eigen::Index>(dt * db);
        Matrix const m = Eigen::Map<Matrix const>(psi.data(), rows, cols);
#if EIGEN_VERSION_AT_LEAST(3, 4, 90)
        Eigen::BDCSVD<Matrix, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(m);
#else
        Eigen::BDCSVD<Matrix> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
#endif
        auto const& values = svd.singularValues();
        if (values.size() == 0 || !(values[0] > 0.))
            throw std::runtime_error("classic DMRG MPS export: the wave function vanishes");
        std::size_t keep = 0;
        while (keep < static_cast<std::size_t>(values.size()) && values[static_cast<Eigen::Index>(keep)] > 1e-14 * values[0])
            ++keep;
        auto const& u = svd.matrixU();
        auto const& v = svd.matrixV();
        left_shape = {da, ds, keep};
        left.assign(da * ds * keep, T(0));
        for (Eigen::Index row = 0; row < rows; ++row)
            for (std::size_t x = 0; x < keep; ++x)
                left[static_cast<std::size_t>(row) * keep + x] = u(row, static_cast<Eigen::Index>(x));
        right_shape = {keep, dt, db};
        right.assign(keep * dt * db, T(0));
        for (std::size_t x = 0; x < keep; ++x)
            for (Eigen::Index col = 0; col < cols; ++col) {
                T value = v(col, static_cast<Eigen::Index>(x));
                if constexpr (!std::is_same_v<T, double>) value = std::conj(value);
                right[x * static_cast<std::size_t>(cols) + static_cast<std::size_t>(col)] =
                    values[static_cast<Eigen::Index>(x)] * value;
            }
    }

    dmtk::System<T>& system_;
    int n_ = 0;
    std::vector<std::vector<std::size_t>> kit_index_; // [site][sorted state] -> PhysicalBasis index
    std::vector<Shape> shapes_;
    std::vector<std::vector<T>> tensors_;
};

} // namespace alps::dmrg::classic::detail
