// SPDX-License-Identifier: MIT
#include <alps/ed.hpp>
#include <alps/models.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>

// The ed component diagonalizes a compiled model without the XML catalog,
// legacy parameters, HDF5 or the simulation runtime: an open 4-site
// Heisenberg chain in the Sz=0 sector (dimension 6) densely, the sparse
// assembly entry point on the same model, the momentum-resolved and the
// point-group-resolved spectra of the periodic chain and a thermodynamic
// average over its levels.
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

    // Periodic chain: every momentum of the Sz=0 sector, whose six levels
    // include the exact ground state energy -2 of the 4-site ring.
    auto const ring = mdl::heisenberg(lat::chain(4), 1.0);
    alps::ed::ModelSpectrumOptions spectrum;
    spectrum.sector = {{"Sz", 0.}};
    spectrum.all_momenta = true;
    auto const resolved = alps::ed::model_spectrum(ring, spectrum);
    std::vector<double> levels;
    for (auto const& block : resolved.blocks)
        levels.insert(levels.end(), block.eigenvalues.begin(), block.eigenvalues.end());
    if (resolved.blocks.size() != 4 || levels.size() != 6)
        throw std::runtime_error("ed component momentum resolution failed");
    if (std::abs(*std::min_element(levels.begin(), levels.end()) + 2.0) > 1e-10)
        throw std::runtime_error("ed component momentum-resolved ground state energy is off");
    auto const thermal = alps::ed::thermal_surface(levels, {}, {1.0}, {0.0});
    if (thermal.size() != 1 || !(thermal[0].entropy > 0.) || !(thermal[0].energy > -2.0))
        throw std::runtime_error("ed component thermodynamic average is off");

    // The complete D4 irrep family of the ring resolves the same six levels.
    std::vector<alps::ed::FiniteGroupIrrep> const family{
        alps::ed::dihedral_irrep(4, 0, 1), alps::ed::dihedral_irrep(4, 0, -1), alps::ed::dihedral_irrep(4, 1),
        alps::ed::dihedral_irrep(4, 2, 1), alps::ed::dihedral_irrep(4, 2, -1)};
    alps::ed::FiniteSpectrumOptions finite;
    finite.sector = {{"Sz", 0.}};
    auto const by_irrep = alps::ed::finite_spectrum(ring, family, finite);
    std::size_t resolved_states = 0;
    double lowest = 0.;
    for (auto const& block : by_irrep.blocks) {
        if (block.group_order != 8 || block.irrep.empty())
            throw std::runtime_error("ed component finite-group block is unlabelled");
        resolved_states += block.eigenvalues.size();
        lowest = std::min(lowest, block.eigenvalues.front());
    }
    if (resolved_states != 6 || std::abs(lowest + 2.0) > 1e-10)
        throw std::runtime_error("ed component finite-group resolution failed");
    return 0;
}
