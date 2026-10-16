// SPDX-License-Identifier: MIT
#include <alps/models.hpp>
#include <cmath>
#include <stdexcept>

// The models component builds and compiles a model without the XML catalog,
// legacy parameters, HDF5 or BLAS/LAPACK, and exports it as Pauli strings.
int main() {
    namespace lat = alps::lattice;
    namespace mdl = alps::models;
    auto const chain = lat::chain(4, 1.0, lat::Boundary::open);
    if (chain.num_sites() != 4 || chain.num_bonds() != 3)
        throw std::runtime_error("models component lattice construction failed");
    auto const model = mdl::heisenberg(chain, 1.0);
    auto const compiled = mdl::compile_hamiltonian(model);
    if (compiled.op.site_dims != std::vector<int>(4, 2) || compiled.edges.size() != 3)
        throw std::runtime_error("models component compilation lost the lattice");
    double weight = 0.;
    for (auto const& term : compiled.op.terms) weight += std::abs(term.coefficient);
    if (compiled.op.terms.empty() || !(weight > 0.))
        throw std::runtime_error("models component compiled no terms");
    mdl::validate_hermiticity(model);
    auto const basis = mdl::physical_basis(model);
    if (basis.sites.size() != 4 || basis.conserved_quantum_numbers != std::vector<std::string>{"Sz"})
        throw std::runtime_error("models component physical basis is incomplete");
    // Nine strings per bond: XX, YY, ZZ each as (1/4) of the Pauli product.
    auto const pauli = mdl::export_pauli(compiled.op, basis);
    if (pauli.num_qubits != 4 || pauli.terms.size() != 9 || pauli.discarded_l1 > 1e-12)
        throw std::runtime_error("models component Pauli export is incomplete");
    for (auto const& [word, coefficient] : pauli.terms)
        if (word.size() != 4 || std::abs(coefficient - .25) > 1e-14)
            throw std::runtime_error("models component Pauli export has a wrong coefficient");
    if (pauli.state_index({0, 1, 0, 1}) != 10)
        throw std::runtime_error("models component Pauli register layout changed");
    return 0;
}
