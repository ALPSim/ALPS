// SPDX-License-Identifier: MIT
#include <alps/dmrg/classic.hpp>
#include <alps/models.hpp>
#include <cmath>
#include <stdexcept>

// The classic DMRG component solves a compiled model without the XML
// catalog, legacy parameters, HDF5, the simulation runtime, LAPACK or Eigen
// in its public interface: the open 6-site Heisenberg chain in the Sz = 0
// sector with enough kept states for exact sweeps, plus one measured
// observable on the ground state.
int main() {
    namespace lat = alps::lattice;
    namespace mdl = alps::models;
    namespace cd = alps::dmrg::classic;
    auto const model = mdl::heisenberg(lat::chain(6, 1.0, lat::Boundary::open), 1.0);
    cd::Options options;
    options.sweeps = cd::fixed_sweeps(3, 32);
    options.warmup_states = 32;
    options.sector = {{"Sz", 0.}};
    cd::Measurements measurements;
    measurements.terms["Sz"] = cd::site_terms(model, "Sz");
    measurements.expectation["H"] = mdl::compile_hamiltonian(model).op;
    auto const result = cd::solve(model, options, measurements);
    // Exact ground state energy of the open 6-site chain.
    if (result.energies.size() != 1 || std::abs(result.energies[0] + 2.493577133888) > 1e-8)
        throw std::runtime_error("classic_dmrg component ground state energy is off");
    if (result.iterations.empty() || !(result.final_truncation_error < 1e-10))
        throw std::runtime_error("classic_dmrg component did not report exact sweeps");
    auto const& sz = result.measurements.at("Sz");
    if (sz.values.size() != 1 || sz.values[0].size() != 6)
        throw std::runtime_error("classic_dmrg component site measurement layout is off");
    double total = 0.;
    for (double value : sz.values[0]) total += value;
    if (std::abs(total) > 1e-8)
        throw std::runtime_error("classic_dmrg component site measurements do not sum to the sector total");
    auto const& energy = result.measurements.at("H");
    if (std::abs(energy.values.at(0).at(0) - result.energies[0]) > 1e-8)
        throw std::runtime_error("classic_dmrg component expectation value disagrees with the energy");
    return 0;
}
