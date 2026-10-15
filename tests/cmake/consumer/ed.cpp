// SPDX-License-Identifier: MIT
#include <alps/ed.hpp>
#include <alps/models.hpp>
#include <cmath>
#include <stdexcept>

// The ed component diagonalizes a compiled model without the XML catalog,
// legacy parameters, HDF5 or the simulation runtime: an open 4-site
// Heisenberg chain in the Sz=0 sector (dimension 6) densely, and the sparse
// assembly entry point on the same model.
int main() {
    namespace lat = alps::lattice;
    namespace mdl = alps::models;
    auto const model = mdl::heisenberg(lat::chain(4, 1.0, lat::Boundary::open), 1.0);
    alps::ed::ProductOptions sector;
    sector.sector = {{"Sz", 0.}};
    auto const problem = alps::ed::prepare_product(model, sector);
    if (problem.dimension() != 6)
        throw std::runtime_error("ed component sector enumeration failed");
    alps::ed::EigenstateOptions request;
    request.count = 2;
    auto const result = problem.eigenstates(request);
    if (result.eigenvalues.size() != 2 || result.method != alps::ed::EigenstateMethod::dense)
        throw std::runtime_error("ed component dense eigensolve failed");
    // E0 of the open 4-site chain: -(3 + 2 sqrt(3)) / 4 ... as a sanity floor
    if (!(result.eigenvalues[0] < -1.6 && result.eigenvalues[0] > -1.7))
        throw std::runtime_error("ed component ground state energy is off");
    auto const csr = alps::ed::build_sparse_from_local(model.to_local_decomposition());
    if (csr.dim != 16 || csr.indptr.size() != 17)
        throw std::runtime_error("ed component sparse assembly failed");
    auto const energy = alps::ed::expectation_values(result, model.operator_decomposition().op);
    if (energy.size() != 2 || std::abs(energy[0].real() - result.eigenvalues[0]) > 1e-10)
        throw std::runtime_error("ed component expectation value disagrees with the eigenvalue");
    return 0;
}
