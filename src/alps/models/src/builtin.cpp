// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/src/operators/builtin.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
//
// The 14 built-in model factories, mirroring builtin.py term-for-term
// (including the Sz-eigenoperator bond forms — Sz⊗Sz + (S+⊗S- + S-⊗S+)/2
// rather than Sx/Sy tensor products — which the QMC matrix-model adapter's
// raising/lowering classifier requires). Parity pinned by
// src/model/tests/models/models_kit_parity.cpp.

#include <alps/models/builtin.hpp>

#include <cmath>
#include <stdexcept>

namespace alps {
namespace models {

namespace {

// site.<op>(i) shorthand.
TermSum L(SiteBasisPtr const& s, char const* op, int i) {
    return TermSum::located(s, op, i);
}

// The Sz-eigenoperator Heisenberg bond: Sz⊗Sz + (Sp⊗Sm + Sm⊗Sp)/2.
TermSum ss_bond(SiteBasisPtr const& s, int i, int j) {
    return L(s, "Sz", i) * L(s, "Sz", j) +
           L(s, "Sp", i) * L(s, "Sm", j) * Complex(0.5) +
           L(s, "Sm", i) * L(s, "Sp", j) * Complex(0.5);
}

} // namespace

Hamiltonian heisenberg(lattice::Lattice const& lat, double J, double h, double S,
                       std::optional<int> bond_type) {
    SiteBasisPtr site = spin(S);
    Hamiltonian H(lat, site, FermionConvention::strict, {"Sz"});
    H.add_bond(J, [site](int i, int j) { return ss_bond(site, i, j); }, bond_type);
    if (h != 0.0)
        H.add_site(-h, [site](int i) { return L(site, "Sz", i); });
    return H;
}

Hamiltonian ising(lattice::Lattice const& lat, double J, double h, double S,
                  std::optional<int> bond_type) {
    SiteBasisPtr site = spin(S);
    Hamiltonian H(lat, site, FermionConvention::strict, {"Sz"});
    H.add_bond(J, [site](int i, int j) { return L(site, "Sz", i) * L(site, "Sz", j); },
               bond_type);
    if (h != 0.0)
        H.add_site(-h, [site](int i) { return L(site, "Sz", i); });
    return H;
}

Hamiltonian xxz(lattice::Lattice const& lat, double Jxy, double Jz, double h, double S,
                std::optional<int> bond_type) {
    SiteBasisPtr site = spin(S);
    Hamiltonian H(lat, site, FermionConvention::strict, {"Sz"});
    H.add_bond(Jxy,
               [site](int i, int j) {
                   return L(site, "Sp", i) * L(site, "Sm", j) * Complex(0.5) +
                          L(site, "Sm", i) * L(site, "Sp", j) * Complex(0.5);
               },
               bond_type);
    H.add_bond(Jz, [site](int i, int j) { return L(site, "Sz", i) * L(site, "Sz", j); },
               bond_type);
    if (h != 0.0)
        H.add_site(-h, [site](int i) { return L(site, "Sz", i); });
    return H;
}

Hamiltonian xy(lattice::Lattice const& lat, double J, double S,
               std::optional<int> bond_type) {
    return xxz(lat, /*Jxy=*/J, /*Jz=*/0.0, /*h=*/0.0, S, bond_type);
}

Hamiltonian transverse_field_ising(lattice::Lattice const& lat, double J, double Gamma,
                                   double h, double S, std::optional<int> bond_type) {
    SiteBasisPtr site = spin(S);
    std::vector<std::string> qns = Gamma == 0.0 ? std::vector<std::string>{"Sz"}
                                                : std::vector<std::string>{};
    Hamiltonian H(lat, site, FermionConvention::strict, std::move(qns));
    H.add_bond(J, [site](int i, int j) { return L(site, "Sz", i) * L(site, "Sz", j); },
               bond_type);
    if (Gamma != 0.0)
        H.add_site(-Gamma, [site](int i) { return L(site, "Sx", i); });
    if (h != 0.0)
        H.add_site(-h, [site](int i) { return L(site, "Sz", i); });
    return H;
}

Hamiltonian heisenberg_d(lattice::Lattice const& lat, double J, double D, double h,
                         double S, std::optional<int> bond_type) {
    SiteBasisPtr site = spin(S);
    Hamiltonian H(lat, site, FermionConvention::strict, {"Sz"});
    H.add_bond(J, [site](int i, int j) { return ss_bond(site, i, j); }, bond_type);
    if (D != 0.0)
        H.add_site(D, [site](int i) { return L(site, "Sz", i) * L(site, "Sz", i); });
    if (h != 0.0)
        H.add_site(-h, [site](int i) { return L(site, "Sz", i); });
    return H;
}

Hamiltonian biquadratic(lattice::Lattice const& lat, double J, double K, double S,
                        std::optional<int> bond_type) {
    SiteBasisPtr site = spin(S);
    Hamiltonian H(lat, site, FermionConvention::strict, {"Sz"});
    H.add_bond(J, [site](int i, int j) { return ss_bond(site, i, j); }, bond_type);
    if (K != 0.0)
        H.add_bond(K,
                   [site](int i, int j) {
                       TermSum const ss = ss_bond(site, i, j);
                       return ss * ss;
                   },
                   bond_type);
    return H;
}

Hamiltonian j1j2(lattice::Lattice const& lat, double J1, double J2, double S) {
    SiteBasisPtr site = spin(S);
    Hamiltonian H(lat, site, FermionConvention::strict, {"Sz"});
    H.add_bond(J1, [site](int i, int j) { return ss_bond(site, i, j); }, 0);
    H.add_bond(J2, [site](int i, int j) { return ss_bond(site, i, j); }, 1);
    return H;
}

Hamiltonian hubbard(lattice::Lattice const& lat, double t, double U, double mu,
                    std::optional<int> bond_type) {
    SiteBasisPtr site = spinful_fermion();
    Hamiltonian H(lat, site, FermionConvention::strict, {"N", "Sz"});
    H.add_bond(-t,
               [site](int i, int j) {
                   return L(site, "c_dag_up", i) * L(site, "c_up", j) +
                          L(site, "c_dag_dn", i) * L(site, "c_dn", j);
               },
               bond_type, /*add_hc=*/true);
    if (U != 0.0)
        H.add_site(U, [site](int i) { return L(site, "n_updn", i); });
    if (mu != 0.0)
        H.add_site(-mu, [site](int i) { return L(site, "n", i); });
    return H;
}

Hamiltonian spinless_fermion_model(lattice::Lattice const& lat, double t, double mu,
                                   double V, std::optional<int> bond_type) {
    SiteBasisPtr site = spinless_fermion();
    Hamiltonian H(lat, site, FermionConvention::strict, {"N"});
    H.add_bond(-t,
               [site](int i, int j) { return L(site, "c_dag", i) * L(site, "c", j); },
               bond_type, /*add_hc=*/true);
    if (V != 0.0)
        H.add_bond(V, [site](int i, int j) { return L(site, "n", i) * L(site, "n", j); },
                   bond_type);
    if (mu != 0.0)
        H.add_site(-mu, [site](int i) { return L(site, "n", i); });
    return H;
}

Hamiltonian t_j(lattice::Lattice const& lat, double t, double J,
                std::optional<int> bond_type) {
    // Projected 3-state space |0>, |up>, |dn> (custom fermionic site).
    int const d = 3;
    auto Z = []() { return Matrix(static_cast<std::size_t>(d) * d, Complex(0.0)); };
    auto set = [](Matrix& m, int r, int c, double v) {
        m[static_cast<std::size_t>(r) * d + c] = v;
    };
    Matrix c_up = Z();    set(c_up, 0, 1, 1.0);
    Matrix c_dag_up = Z(); set(c_dag_up, 1, 0, 1.0);
    Matrix c_dn = Z();    set(c_dn, 0, 2, 1.0);
    Matrix c_dag_dn = Z(); set(c_dag_dn, 2, 0, 1.0);
    Matrix n_up = Z();    set(n_up, 1, 1, 1.0);
    Matrix n_dn = Z();    set(n_dn, 2, 2, 1.0);
    Matrix n_tot = Z();   set(n_tot, 1, 1, 1.0); set(n_tot, 2, 2, 1.0);
    Matrix Sz = Z();      set(Sz, 1, 1, 0.5); set(Sz, 2, 2, -0.5);
    Matrix Sp = Z();      set(Sp, 1, 2, 1.0);   // |up><dn|
    Matrix Sm = Z();      set(Sm, 2, 1, 1.0);   // |dn><up|
    Matrix Sx = Z();      set(Sx, 1, 2, 0.5); set(Sx, 2, 1, 0.5);
    Matrix Sy = Z();
    Sy[1 * 3 + 2] = Complex(0.0, -0.5);
    Sy[2 * 3 + 1] = Complex(0.0, 0.5);

    SiteBasis sb;
    sb.dim = d;
    sb.statistics = Statistics::fermionic;
    sb.quantum_numbers["N"] = {0.0, 1.0, 1.0};
    sb.quantum_numbers["Sz"] = {0.0, 0.5, -0.5};
    sb.parity_operator = {1.0, -1.0, -1.0};
    auto add = [&sb](char const* name, Matrix m, Parity p = Parity::even) {
        LocalOp op;
        op.name = name;
        op.dim = d;
        op.matrix = std::move(m);
        op.parity = p;
        sb.operators.emplace(name, std::move(op));
    };
    add("c_up", std::move(c_up), Parity::odd);
    add("c_dag_up", std::move(c_dag_up), Parity::odd);
    add("c_dn", std::move(c_dn), Parity::odd);
    add("c_dag_dn", std::move(c_dag_dn), Parity::odd);
    add("n_up", std::move(n_up));
    add("n_dn", std::move(n_dn));
    add("n", std::move(n_tot));
    add("Sz", std::move(Sz));
    add("Sp", std::move(Sp));
    add("Sm", std::move(Sm));
    add("Sx", std::move(Sx));
    add("Sy", std::move(Sy));
    SiteBasisPtr site = finalize_site(std::move(sb));

    Hamiltonian H(lat, site, FermionConvention::strict, {"N", "Sz"});
    H.add_bond(-t,
               [site](int i, int j) {
                   return L(site, "c_dag_up", i) * L(site, "c_up", j) +
                          L(site, "c_dag_dn", i) * L(site, "c_dn", j);
               },
               bond_type, /*add_hc=*/true);
    H.add_bond(J, [site](int i, int j) { return ss_bond(site, i, j); }, bond_type);
    H.add_bond(-J / 4.0,
               [site](int i, int j) { return L(site, "n", i) * L(site, "n", j); },
               bond_type);
    return H;
}

Hamiltonian kondo_lattice(lattice::Lattice const& lat, double t, double J, double mu) {
    SiteBasisPtr f_site = spinful_fermion();
    SiteBasisPtr s_site = spin(0.5);
    Hamiltonian H(lat, {{0, f_site}, {1, s_site}}, FermionConvention::strict,
                  {"N", "Sz"});
    H.add_bond(-t,
               [f_site](int i, int j) {
                   return L(f_site, "c_dag_up", i) * L(f_site, "c_up", j) +
                          L(f_site, "c_dag_dn", i) * L(f_site, "c_dn", j);
               },
               0, /*add_hc=*/true);
    H.add_bond(J,
               [f_site, s_site](int i, int j) {
                   return L(f_site, "Sz", i) * L(s_site, "Sz", j) +
                          L(f_site, "Sp", i) * L(s_site, "Sm", j) * Complex(0.5) +
                          L(f_site, "Sm", i) * L(s_site, "Sp", j) * Complex(0.5);
               },
               1);
    if (mu != 0.0)
        H.add_site(-mu, [f_site](int i) { return L(f_site, "n", i); }, 0);
    return H;
}

Hamiltonian hardcore_boson_model(lattice::Lattice const& lat, double t, double mu,
                                 double V, std::optional<int> bond_type) {
    SiteBasisPtr site = hardcore_boson();
    Hamiltonian H(lat, site, FermionConvention::strict, {"N"});
    H.add_bond(-t,
               [site](int i, int j) {
                   return L(site, "a_dag", i) * L(site, "a", j) +
                          L(site, "a", i) * L(site, "a_dag", j);
               },
               bond_type);
    if (V != 0.0)
        H.add_bond(V, [site](int i, int j) { return L(site, "n", i) * L(site, "n", j); },
                   bond_type);
    if (mu != 0.0)
        H.add_site(-mu, [site](int i) { return L(site, "n", i); });
    return H;
}

Hamiltonian bose_hubbard(lattice::Lattice const& lat, double t, double U, double mu,
                         double V, int n_max, std::optional<int> bond_type) {
    SiteBasisPtr site = boson(n_max);
    Hamiltonian H(lat, site, FermionConvention::strict, {"N"});
    H.add_bond(-t,
               [site](int i, int j) {
                   return L(site, "b_dag", i) * L(site, "b", j) +
                          L(site, "b_dag", j) * L(site, "b", i);
               },
               bond_type);
    if (U != 0.0)
        H.add_site(U / 2.0, [site](int i) {
            return L(site, "n", i) * L(site, "n", i) - L(site, "n", i);
        });
    if (V != 0.0)
        H.add_bond(V, [site](int i, int j) { return L(site, "n", i) * L(site, "n", j); },
                   bond_type);
    if (mu != 0.0)
        H.add_site(-mu, [site](int i) { return L(site, "n", i); });
    return H;
}

Hamiltonian ab_initio(std::vector<double> const& h1, std::vector<double> const& eri, int n, double tol) {
    if (n < 1) throw std::invalid_argument("ab_initio: the orbital count must be positive");
    auto const m = static_cast<std::size_t>(n);
    if (h1.size() != m * m) throw std::invalid_argument("ab_initio: h1 must hold n x n values");
    if (eri.size() != m * m * m * m) throw std::invalid_argument("ab_initio: eri must hold n^4 values");
    if (!(tol >= 0.)) throw std::invalid_argument("ab_initio: the integral cutoff must be nonnegative");
    SiteBasisPtr site = spinful_fermion();
    Hamiltonian H(lattice::chain(n, 1.0, lattice::Boundary::open), site, FermionConvention::strict, {"N", "Sz"});
    auto const kept = [tol](double value) { return value != 0.0 && std::abs(value) >= tol; };
    char const* const create[] = {"c_dag_up", "c_dag_dn"};
    char const* const annihilate[] = {"c_up", "c_dn"};
    TermSum total;
    // One body: sum_{pq,s} h_pq c+_{ps} c_{qs}.
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q) {
            double const value = h1[std::size_t(p) * m + std::size_t(q)];
            if (!kept(value)) continue;
            for (int s = 0; s < 2; ++s) total += Complex(value) * (L(site, create[s], p) * L(site, annihilate[s], q));
        }
    // Two body: 1/2 sum_{pqrs,st} (pq|rs) c+_{ps} c+_{rt} c_{st} c_{qs}.
    for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q)
            for (int r = 0; r < n; ++r)
                for (int o = 0; o < n; ++o) {
                    double const value = eri[((std::size_t(p) * m + std::size_t(q)) * m + std::size_t(r)) * m + std::size_t(o)];
                    if (!kept(value)) continue;
                    for (int s1 = 0; s1 < 2; ++s1)
                        for (int s2 = 0; s2 < 2; ++s2)
                            total += Complex(0.5 * value) * (L(site, create[s1], p) * L(site, create[s2], r)
                                                             * L(site, annihilate[s2], o) * L(site, annihilate[s1], q));
                }
    H.add(std::move(total));
    return H;
}

} // namespace models
} // namespace alps
