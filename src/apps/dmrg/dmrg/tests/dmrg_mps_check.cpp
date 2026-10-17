// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Checks the /mps group a dmrg run with SAVE_MPS=1 writes (the layout of
// ALPSim/ALPS PR 216): the tensor shapes chain from bond one to bond one,
// the tensors before mps/center are left-orthonormal and those after it
// right-orthonormal, the state is normalized, and every product state with
// weight has the requested total of one quantum number, read from
// mps/site_basis in the local-state order of the tensors.
//
// usage: dmrg_mps_check <result.h5> <quantum number> <total>
#include <alps/hdf5.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Tensor {
    std::size_t left, local, right;
    std::vector<double> values;
    double at(std::size_t l, std::size_t s, std::size_t r) const { return values[(l * local + s) * right + r]; }
};

int fail(std::string const& message) {
    std::cerr << "dmrg_mps_check: " << message << '\n';
    return 1;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 4) return fail("usage: dmrg_mps_check <result.h5> <quantum number> <total>");
    std::string const name = argv[2];
    double const total = std::atof(argv[3]);
    alps::hdf5::archive ar(argv[1], "r");
    if (!ar.is_group("/mps")) return fail("no /mps group");
    int length = 0, center = 0;
    ar["/mps/length"] >> length;
    ar["/mps/center"] >> center;
    std::vector<int> types;
    ar["/mps/site_type"] >> types;
    if (length < 4 || static_cast<int>(types.size()) != length || center < 0 || center >= length)
        return fail("inconsistent length, center or site types");
    std::vector<Tensor> sites;
    std::size_t bond = 1;
    for (int i = 0; i < length; ++i) {
        std::string const path = "/mps/tensors/" + std::to_string(i);
        std::vector<std::size_t> shape;
        Tensor t{};
        ar[path + "/shape"] >> shape;
        ar[path + "/values"] >> t.values;
        if (shape.size() != 3 || shape[0] != bond) return fail("the shapes do not chain at site " + std::to_string(i));
        t.left = shape[0], t.local = shape[1], t.right = shape[2];
        if (t.values.size() != t.left * t.local * t.right) return fail("wrong value count at site " + std::to_string(i));
        bond = t.right;
        sites.push_back(std::move(t));
    }
    if (bond != 1) return fail("the last bond is not one");

    double worst = 0.;
    for (int i = 0; i < length; ++i) {
        if (i == center) continue;
        auto const& t = sites[static_cast<std::size_t>(i)];
        bool const left = i < center;
        std::size_t const k = left ? t.right : t.left, other = left ? t.left : t.right;
        for (std::size_t x = 0; x < k; ++x)
            for (std::size_t y = 0; y < k; ++y) {
                double sum = 0.;
                for (std::size_t a = 0; a < other; ++a)
                    for (std::size_t s = 0; s < t.local; ++s)
                        sum += left ? t.at(a, s, x) * t.at(a, s, y) : t.at(x, s, a) * t.at(y, s, a);
                worst = std::max(worst, std::abs(sum - (x == y ? 1. : 0.)));
            }
    }
    if (worst > 1e-10) return fail("a tensor is not orthonormal: " + std::to_string(worst));

    // The dense state with the charge total of every product state.
    std::vector<double> amplitude{1.}, charge{0.};
    std::size_t width = 1;
    for (int i = 0; i < length; ++i) {
        auto const& t = sites[static_cast<std::size_t>(i)];
        std::vector<double> local;
        ar["/mps/site_basis/" + std::to_string(types[static_cast<std::size_t>(i)]) + "/" + name] >> local;
        if (local.size() != t.local) return fail("site basis and tensor disagree at site " + std::to_string(i));
        std::size_t const count = amplitude.size() / width;
        std::vector<double> next(count * t.local * t.right, 0.), next_charge(count * t.local, 0.);
        for (std::size_t c = 0; c < count; ++c)
            for (std::size_t s = 0; s < t.local; ++s) {
                next_charge[c * t.local + s] = charge[c] + local[s];
                for (std::size_t l = 0; l < t.left; ++l)
                    for (std::size_t r = 0; r < t.right; ++r)
                        next[(c * t.local + s) * t.right + r] += amplitude[c * width + l] * t.at(l, s, r);
            }
        amplitude.swap(next);
        charge.swap(next_charge);
        width = t.right;
    }
    double norm = 0., outside = 0.;
    for (std::size_t c = 0; c < amplitude.size(); ++c) {
        norm += amplitude[c] * amplitude[c];
        if (std::abs(charge[c] - total) > 1e-9) outside += amplitude[c] * amplitude[c];
    }
    if (std::abs(norm - 1.) > 1e-10) return fail("the state is not normalized: " + std::to_string(norm));
    if (outside > 1e-20) return fail("weight outside " + name + " = " + argv[3] + ": " + std::to_string(outside));
    std::cout << "mps ok: " << length << " sites, center " << center << '\n';
    return 0;
}
