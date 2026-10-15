// SPDX-License-Identifier: MIT
#include <alps/ed.hpp>
#include <alps/model/catalog.hpp>
#include <cmath>
#include <stdexcept>
#include <vector>

// The catalog importer of the legacy model module compiles an XML catalog
// model (models.xml on lattices.xml) into an ALPS::models operator that
// ALPS::ed diagonalizes in the catalog's own conserved sector: the 4-site
// periodic Heisenberg ring at Sz=0 has the ground state energy -2.
int main() {
    alps::Parameters p;
    p["MODEL"] = "spin";
    p["LATTICE"] = "chain lattice";
    p["L"] = 4;
    p["J"] = 1.0;
    p["CONSERVED_QUANTUMNUMBERS"] = "Sz";
    p["Sz_total"] = 0;
    alps::models::CatalogModel<> model(p);
    auto const op = model.compile_hamiltonian();
    if (op.site_dims != std::vector<int>(4, 2) || op.terms.empty())
        throw std::runtime_error("catalog importer compiled the wrong local structure");
    alps::ed::ProjectedSpaces spaces;
    spaces.charges = model.conserved_quantum_numbers();
    auto const parity = model.local_parity();
    std::vector<std::vector<std::vector<double>>> charges;
    for (auto const& name : spaces.charges) charges.push_back(model.local_charge(name));
    for (std::size_t site = 0; site < parity.size(); ++site) {
        alps::ed::SpectrumLocalSpace space;
        space.parity = parity[site];
        for (auto const& charge : charges) space.quantum_numbers.push_back(charge[site]);
        spaces.local_spaces.push_back(space);
    }
    alps::ed::ProjectedOptions options;
    options.sector = model.sector();
    auto const problem = alps::ed::prepare_projected(op, spaces, {}, options);
    if (problem.dimension() != 6)
        throw std::runtime_error("catalog importer sector enumeration is off");
    alps::ed::EigenstateOptions solve;
    solve.count = 1;
    auto const result = problem.eigenstates(solve);
    if (std::abs(result.spectrum.eigenvalues.at(0) + 2.0) > 1e-10)
        throw std::runtime_error("catalog importer ground state energy is off");
    return 0;
}
