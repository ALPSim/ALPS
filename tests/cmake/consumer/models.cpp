// SPDX-License-Identifier: MIT
#include <alps/models.hpp>
#include <cmath>
#include <stdexcept>

// The models component builds and compiles a model without the XML catalog,
// legacy parameters, HDF5 or BLAS/LAPACK.
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
    return 0;
}
