// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/src/operators/site.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
//
// Shared implementation of alps::models site bases, consumed by Python;
// parity pinned by
// src/model/tests/models/models_kit_parity.cpp.

#include <alps/models/site.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace alps {
namespace models {

namespace {

Matrix zeros(int dim) {
    return Matrix(static_cast<std::size_t>(dim) * dim, Complex(0.0, 0.0));
}

Matrix identity(int dim) {
    Matrix m = zeros(dim);
    for (int i = 0; i < dim; ++i) m[static_cast<std::size_t>(i) * dim + i] = 1.0;
    return m;
}

Complex& at(Matrix& m, int dim, int r, int c) {
    return m[static_cast<std::size_t>(r) * dim + c];
}

Matrix dagger(Matrix const& m, int dim) {
    Matrix out = zeros(dim);
    for (int r = 0; r < dim; ++r)
        for (int c = 0; c < dim; ++c)
            at(out, dim, c, r) = std::conj(m[static_cast<std::size_t>(r) * dim + c]);
    return out;
}

// "%+g"-formatted delta tag matching Python's f"{v:+g}".
std::string fmt_plus_g(double v) {
    char buf[48];
    std::snprintf(buf, sizeof(buf), "%+g", v);
    return buf;
}

} // namespace

std::map<std::string, double> compute_qn_deltas(
    Matrix const& matrix, int dim,
    std::map<std::string, std::vector<double>> const& quantum_numbers,
    double tol) {
    std::map<std::string, double> out;
    if (quantum_numbers.empty()) return out;
    // Non-zero entry list.
    std::vector<std::pair<int, int>> nz;
    for (int r = 0; r < dim; ++r)
        for (int c = 0; c < dim; ++c)
            if (matrix[static_cast<std::size_t>(r) * dim + c] != Complex{})
                nz.emplace_back(r, c);
    for (auto const& [qn, arr] : quantum_numbers) {
        if (static_cast<int>(arr.size()) != dim)
            throw std::invalid_argument("compute_qn_deltas: QN '" + qn +
                                        "' length != matrix dim");
        if (nz.empty()) { out[qn] = 0.0; continue; }
        double const ref = arr[static_cast<std::size_t>(nz[0].first)] -
                           arr[static_cast<std::size_t>(nz[0].second)];
        bool ok = true;
        for (auto const& [r, c] : nz) {
            double const d = arr[static_cast<std::size_t>(r)] - arr[static_cast<std::size_t>(c)];
            if (!(std::abs(d - ref) <= tol)) { ok = false; break; }
        }
        if (ok) out[qn] = ref;
    }
    return out;
}

bool LocalOp::is_identity(double tol) const {
    for (int r = 0; r < dim; ++r)
        for (int c = 0; c < dim; ++c) {
            Complex const expect = (r == c) ? Complex(1.0) : Complex(0.0);
            if (std::abs(matrix[static_cast<std::size_t>(r) * dim + c] - expect) > tol)
                return false;
        }
    return true;
}

LocalOp LocalOp::dag() const {
    LocalOp out;
    bool const has_suffix =
        name.size() >= 4 && name.compare(name.size() - 4, 4, "_dag") == 0;
    out.name = has_suffix ? name.substr(0, name.size() - 4) : name + "_dag";
    out.dim = dim;
    out.matrix = dagger(matrix, dim);
    out.parity = parity;
    for (auto const& [k, v] : qn_deltas) out.qn_deltas[k] = -v + 0.0;
    return out;
}

std::vector<LocalOp> LocalOp::split_by_qn_delta(
    std::map<std::string, std::vector<double>> const& quantum_numbers) const {
    std::vector<LocalOp> out;
    if (quantum_numbers.empty()) {
        LocalOp copy = *this;
        copy.qn_deltas.clear();
        out.push_back(std::move(copy));
        return out;
    }
    // Bucket non-zero entries by their delta tuple (QNs in map order —
    // sorted; see the SiteBasis ordering note in site.hpp).
    using Exact = boost::multiprecision::cpp_rational;
    std::map<std::vector<Exact>, Matrix> buckets;
    for (int r = 0; r < dim; ++r)
        for (int c = 0; c < dim; ++c) {
            Complex const v = matrix[static_cast<std::size_t>(r) * dim + c];
            if (v == Complex{}) continue;
            std::vector<Exact> key;
            key.reserve(quantum_numbers.size());
            for (auto const& [qn, arr] : quantum_numbers)
                key.push_back(Exact(arr[static_cast<std::size_t>(r)]) -
                              Exact(arr[static_cast<std::size_t>(c)]));
            auto it = buckets.find(key);
            if (it == buckets.end())
                it = buckets.emplace(std::move(key), zeros(dim)).first;
            it->second[static_cast<std::size_t>(r) * dim + c] = v;
        }
    if (buckets.empty()) {
        LocalOp zero;
        zero.name = name;
        zero.dim = dim;
        zero.matrix = zeros(dim);
        zero.parity = parity;
        for (auto const& [qn, arr] : quantum_numbers) { (void)arr; zero.qn_deltas[qn] = 0.0; }
        out.push_back(std::move(zero));
        return out;
    }
    for (auto const& [key, sub] : buckets) {   // sorted by delta tuple
        LocalOp comp;
        std::string tag;
        std::size_t i = 0;
        comp.qn_deltas.clear();
        for (auto const& [qn, arr] : quantum_numbers) {
            (void)arr;
            double const delta = key[i].convert_to<double>();
            comp.qn_deltas[qn] = delta;
            if (i) tag += ",";
            tag += qn + "=" + fmt_plus_g(delta);
            ++i;
        }
        comp.name = buckets.size() > 1 ? name + "[" + tag + "]" : name;
        comp.dim = dim;
        comp.matrix = sub;
        comp.parity = parity;
        out.push_back(std::move(comp));
    }
    return out;
}

LocalOp const& SiteBasis::op(std::string const& name) const {
    auto it = operators.find(name);
    if (it == operators.end())
        throw std::out_of_range("SiteBasis has no operator '" + name + "'");
    return it->second;
}

SiteBasisPtr finalize_site(SiteBasis site) {
    if (site.parity_operator.empty())
        site.parity_operator.assign(static_cast<std::size_t>(site.dim), 1.0);
    if (!site.operators.count("I")) {
        LocalOp I;
        I.name = "I";
        I.dim = site.dim;
        I.matrix = identity(site.dim);
        site.operators.emplace("I", std::move(I));
    }
    if (!site.quantum_numbers.empty())
        for (auto& [name, op] : site.operators)
            if (op.qn_deltas.empty())
                op.qn_deltas = compute_qn_deltas(op.matrix, op.dim, site.quantum_numbers);
    return std::make_shared<SiteBasis const>(std::move(site));
}

namespace {

LocalOp make_op(std::string name, int dim, Matrix m, Parity p = Parity::even) {
    LocalOp op;
    op.name = std::move(name);
    op.dim = dim;
    op.matrix = std::move(m);
    op.parity = p;
    return op;
}

Matrix matmul(Matrix const& a, Matrix const& b, int dim) {
    Matrix out = zeros(dim);
    for (int r = 0; r < dim; ++r)
        for (int k = 0; k < dim; ++k) {
            Complex const ark = a[static_cast<std::size_t>(r) * dim + k];
            if (ark == Complex(0.0)) continue;
            for (int c = 0; c < dim; ++c)
                at(out, dim, r, c) += ark * b[static_cast<std::size_t>(k) * dim + c];
        }
    return out;
}

} // namespace

SiteBasisPtr spin(double S) {
    double const twoS = 2.0 * S;
    if (!(S > 0.0) || std::abs(twoS - std::nearbyint(twoS)) > 1e-12)
        throw std::invalid_argument("spin: S must be a positive (half-)integer");
    int const d = static_cast<int>(std::lround(twoS)) + 1;

    std::vector<double> m_values(static_cast<std::size_t>(d));
    for (int i = 0; i < d; ++i) m_values[static_cast<std::size_t>(i)] = S - i;

    Matrix Sz = zeros(d), Sp = zeros(d);
    for (int i = 0; i < d; ++i) {
        at(Sz, d, i, i) = m_values[static_cast<std::size_t>(i)];
        double const m = m_values[static_cast<std::size_t>(i)];
        if (i - 1 >= 0)
            at(Sp, d, i - 1, i) = std::sqrt(S * (S + 1) - m * (m + 1));
    }
    Matrix Sm = dagger(Sp, d);
    Matrix Sx = zeros(d), Sy = zeros(d);
    for (std::size_t k = 0; k < Sx.size(); ++k) {
        Sx[k] = 0.5 * (Sp[k] + Sm[k]);
        Sy[k] = Complex(0.0, -0.5) * (Sp[k] - Sm[k]);
    }
    Matrix S2 = matmul(Sx, Sx, d);
    {
        Matrix const yy = matmul(Sy, Sy, d), zz = matmul(Sz, Sz, d);
        for (std::size_t k = 0; k < S2.size(); ++k) S2[k] += yy[k] + zz[k];
    }

    SiteBasis site;
    site.dim = d;
    site.statistics = Statistics::bosonic;
    site.quantum_numbers["Sz"] = m_values;
    site.operators.emplace("Sz", make_op("Sz", d, std::move(Sz)));
    site.operators.emplace("Sp", make_op("Sp", d, std::move(Sp)));
    site.operators.emplace("Sm", make_op("Sm", d, std::move(Sm)));
    site.operators.emplace("Sx", make_op("Sx", d, std::move(Sx)));
    site.operators.emplace("Sy", make_op("Sy", d, std::move(Sy)));
    site.operators.emplace("S2", make_op("S2", d, std::move(S2)));
    return finalize_site(std::move(site));
}

SiteBasisPtr boson(int n_max) {
    if (n_max < 1) throw std::invalid_argument("boson: n_max must be >= 1");
    int const d = n_max + 1;
    Matrix a = zeros(d), n_op = zeros(d);
    std::vector<double> n_diag(static_cast<std::size_t>(d));
    for (int n = 0; n < d; ++n) {
        n_diag[static_cast<std::size_t>(n)] = n;
        at(n_op, d, n, n) = n;
        if (n >= 1) at(a, d, n - 1, n) = std::sqrt(static_cast<double>(n));
    }
    Matrix a_dag = dagger(a, d);

    SiteBasis site;
    site.dim = d;
    site.statistics = Statistics::bosonic;
    site.quantum_numbers["N"] = n_diag;
    site.operators.emplace("a", make_op("a", d, a));
    site.operators.emplace("a_dag", make_op("a_dag", d, a_dag));
    site.operators.emplace("n", make_op("n", d, std::move(n_op)));
    site.operators.emplace("b", make_op("b", d, std::move(a)));
    site.operators.emplace("b_dag", make_op("b_dag", d, std::move(a_dag)));
    return finalize_site(std::move(site));
}

SiteBasisPtr hardcore_boson() {
    int const d = 2;
    Matrix a = zeros(d);
    at(a, d, 0, 1) = 1.0;
    Matrix a_dag = dagger(a, d);
    Matrix n_op = zeros(d);
    at(n_op, d, 1, 1) = 1.0;

    SiteBasis site;
    site.dim = d;
    site.statistics = Statistics::bosonic;
    site.quantum_numbers["N"] = {0.0, 1.0};
    site.operators.emplace("a", make_op("a", d, a));
    site.operators.emplace("a_dag", make_op("a_dag", d, a_dag));
    site.operators.emplace("n", make_op("n", d, std::move(n_op)));
    site.operators.emplace("b", make_op("b", d, std::move(a)));
    site.operators.emplace("b_dag", make_op("b_dag", d, std::move(a_dag)));
    return finalize_site(std::move(site));
}

SiteBasisPtr spinless_fermion() {
    int const d = 2;
    Matrix c = zeros(d);
    at(c, d, 0, 1) = 1.0;
    Matrix c_dag = dagger(c, d);
    Matrix n_op = zeros(d);
    at(n_op, d, 1, 1) = 1.0;
    Matrix F = zeros(d);
    at(F, d, 0, 0) = 1.0;
    at(F, d, 1, 1) = -1.0;

    SiteBasis site;
    site.dim = d;
    site.statistics = Statistics::fermionic;
    site.quantum_numbers["N"] = {0.0, 1.0};
    site.parity_operator = {1.0, -1.0};
    site.operators.emplace("c", make_op("c", d, std::move(c), Parity::odd));
    site.operators.emplace("c_dag", make_op("c_dag", d, std::move(c_dag), Parity::odd));
    site.operators.emplace("n", make_op("n", d, std::move(n_op)));
    site.operators.emplace("F", make_op("F", d, std::move(F)));
    return finalize_site(std::move(site));
}

SiteBasisPtr spinful_fermion() {
    int const d = 4;
    // Basis |0>=0, |up>=1, |dn>=2, |up,dn>=3 ; c_dag_up|dn> = -|up,dn>.
    Matrix c_dag_up = zeros(d);
    at(c_dag_up, d, 1, 0) = 1.0;
    at(c_dag_up, d, 3, 2) = -1.0;
    Matrix c_up = dagger(c_dag_up, d);   // real -> transpose
    Matrix c_dag_dn = zeros(d);
    at(c_dag_dn, d, 2, 0) = 1.0;
    at(c_dag_dn, d, 3, 1) = 1.0;
    Matrix c_dn = dagger(c_dag_dn, d);

    Matrix n_up = matmul(c_dag_up, c_up, d);
    Matrix n_dn = matmul(c_dag_dn, c_dn, d);
    Matrix n_tot = zeros(d);
    for (std::size_t k = 0; k < n_tot.size(); ++k) n_tot[k] = n_up[k] + n_dn[k];
    Matrix n_updn = matmul(n_up, n_dn, d);

    Matrix Sz = zeros(d);
    for (std::size_t k = 0; k < Sz.size(); ++k) Sz[k] = 0.5 * (n_up[k] - n_dn[k]);
    Matrix Sp = matmul(c_dag_up, c_dn, d);
    Matrix Sm = matmul(c_dag_dn, c_up, d);
    Matrix Sx = zeros(d), Sy = zeros(d);
    for (std::size_t k = 0; k < Sx.size(); ++k) {
        Sx[k] = 0.5 * (Sp[k] + Sm[k]);
        Sy[k] = Complex(0.0, -0.5) * (Sp[k] - Sm[k]);
    }
    Matrix F = zeros(d);
    at(F, d, 0, 0) = 1.0; at(F, d, 1, 1) = -1.0;
    at(F, d, 2, 2) = -1.0; at(F, d, 3, 3) = 1.0;

    SiteBasis site;
    site.dim = d;
    site.statistics = Statistics::fermionic;
    site.quantum_numbers["N"] = {0.0, 1.0, 1.0, 2.0};
    site.quantum_numbers["Sz"] = {0.0, 0.5, -0.5, 0.0};
    site.parity_operator = {1.0, -1.0, -1.0, 1.0};
    site.operators.emplace("c_up", make_op("c_up", d, std::move(c_up), Parity::odd));
    site.operators.emplace("c_dag_up", make_op("c_dag_up", d, std::move(c_dag_up), Parity::odd));
    site.operators.emplace("c_dn", make_op("c_dn", d, std::move(c_dn), Parity::odd));
    site.operators.emplace("c_dag_dn", make_op("c_dag_dn", d, std::move(c_dag_dn), Parity::odd));
    site.operators.emplace("n_up", make_op("n_up", d, std::move(n_up)));
    site.operators.emplace("n_dn", make_op("n_dn", d, std::move(n_dn)));
    site.operators.emplace("n", make_op("n", d, std::move(n_tot)));
    site.operators.emplace("n_updn", make_op("n_updn", d, std::move(n_updn)));
    site.operators.emplace("Sz", make_op("Sz", d, std::move(Sz)));
    site.operators.emplace("Sp", make_op("Sp", d, std::move(Sp)));
    site.operators.emplace("Sm", make_op("Sm", d, std::move(Sm)));
    site.operators.emplace("Sx", make_op("Sx", d, std::move(Sx)));
    site.operators.emplace("Sy", make_op("Sy", d, std::move(Sy)));
    site.operators.emplace("F", make_op("F", d, std::move(F)));
    return finalize_site(std::move(site));
}

} // namespace models
} // namespace alps
