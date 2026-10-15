// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/runtime/projected_problem.cpp.
// Adaptation: C++17 (no std::span or std::numbers), <alps/ed/...> and
// <alps/models/...> include names, namespace alps::ed.
#include <alps/ed/projected_problem.hpp>
#include <alps/models/physical_basis.hpp>
#include <alps/models/site.hpp>
#include <alps/models/terms.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace alps::ed {

MomentumSelection select_momentum(std::vector<std::size_t> const& periods, std::vector<double> const& wavevector) {
    if (wavevector.size() != periods.size())
        throw std::invalid_argument("momentum needs one component per periodic lattice direction");
    double const turn = 2. * std::acos(-1.);
    MomentumSelection result;
    result.momentum.resize(periods.size());
    std::size_t stride = 1;
    for (std::size_t direction = periods.size(); direction-- > 0;) {
        auto const period = periods[direction];
        if (!period) throw std::invalid_argument("translation periods must be positive");
        // Grid indices stay exact in double arithmetic and in the flattened index.
        if (period > (std::size_t{1} << 52) || period > std::numeric_limits<std::size_t>::max() / stride)
            throw std::invalid_argument("translation group is too large");
        if (!std::isfinite(wavevector[direction]))
            throw std::invalid_argument("momentum values must be finite real numbers");
        // Floor-modulo wrapping, round-half-even snapping and an angular
        // tolerance, evaluated in the order the grid convention defines.
        double wrapped = std::fmod(wavevector[direction], turn);
        if (wrapped < 0.) wrapped += turn;
        auto const index = static_cast<std::size_t>(std::nearbyint(wrapped * static_cast<double>(period) / turn)) % period;
        double const canonical = turn * static_cast<double>(index) / static_cast<double>(period);
        double const difference = wrapped - canonical;
        if (std::abs(std::atan2(std::sin(difference), std::cos(difference))) > 1e-10)
            throw std::invalid_argument("momentum is not on the lattice's allowed translation grid");
        result.momentum[direction] = canonical;
        result.index += index * stride;
        stride *= period;
    }
    return result;
}

ProjectedSpaces projected_spaces(models::PhysicalBasis const& physical) {
    models::validate_basis(physical);
    ProjectedSpaces result;
    result.charges = physical.conserved_quantum_numbers;
    std::sort(result.charges.begin(), result.charges.end());
    bool const strict = physical.fermion_convention == models::FermionConvention::strict;
    for (auto const& site : physical.sites) {
        auto& local = result.local_spaces.emplace_back();
        local.parity = strict ? site.parity : std::vector<double>(site.parity.size(), 1.);
        for (auto const& name : result.charges) {
            auto const row = site.quantum_numbers.find(name);
            if (row == site.quantum_numbers.end())
                throw std::invalid_argument("physical basis site lacks declared conserved quantum number " + name);
            local.quantum_numbers.push_back(row->second);
        }
    }
    return result;
}

namespace {
void check_sector(std::vector<std::string> const& conserved, std::map<std::string, double> const& sector) {
    for (auto const& item : sector) {
        if (std::find(conserved.begin(), conserved.end(), item.first) == conserved.end())
            throw std::invalid_argument("sector contains a charge not declared conserved: " + item.first);
        if (!std::isfinite(item.second)) throw std::invalid_argument("sector totals must be finite");
    }
}
std::vector<std::size_t> periods_of(std::vector<SpectrumTranslation> const& translations) {
    std::vector<std::size_t> periods;
    for (auto const& translation : translations) periods.push_back(translation.period);
    return periods;
}
struct PreparedModel {
    models::CompiledOperator op;
    ProjectedSpaces spaces;
    std::vector<std::optional<double>> targets;
};

// A compiled operator already carries its Jordan-Wigner strings. Materialize
// those as ordinary diagonal factors for the shared model validators, using
// ungraded local spaces so no second exchange sign/string can be introduced.
// This is a local-tensor conversion, never a Hilbert-space enumeration. Keep
// the original compiled operator for assembly and its accumulation order.
void validate_compiled_model(models::CompiledOperator const& op, ProjectedSpaces const& spaces) {
    if (!std::is_sorted(spaces.charges.begin(), spaces.charges.end())
        || std::adjacent_find(spaces.charges.begin(), spaces.charges.end()) != spaces.charges.end())
        throw std::invalid_argument("compiled ED charge names must be sorted and unique");
    if (op.site_dims.size() != spaces.local_spaces.size())
        throw std::invalid_argument("compiled ED local space count mismatch");
    std::vector<models::SiteBasisPtr> sites;
    for (std::size_t i = 0; i < op.site_dims.size(); ++i) {
        auto const& local = spaces.local_spaces[i];
        int const dim = op.site_dims[i];
        if (dim <= 0 || local.parity.size() != static_cast<std::size_t>(dim)
            || local.quantum_numbers.size() != spaces.charges.size())
            throw std::invalid_argument("compiled ED local space dimensions mismatch");
        for (auto p : local.parity)
            if (p != 1. && p != -1.) throw std::invalid_argument("local parity must be +1 or -1");
        models::SiteBasis site;
        site.dim = dim;
        for (std::size_t q = 0; q < spaces.charges.size(); ++q)
            site.quantum_numbers.emplace(spaces.charges[q], local.quantum_numbers[q]);
        sites.push_back(models::finalize_site(std::move(site)));
    }
    // Validate all rows, including unused charges and empty operators.
    models::physical_basis(sites, spaces.charges, models::FermionConvention::lattice);
    std::vector<models::Term> terms;
    for (auto const& source : op.terms) {
        auto& term = terms.emplace_back();
        term.coeff = source.coefficient;
        auto append = [&](int site, models::Matrix matrix) {
            if (site < 0 || static_cast<std::size_t>(site) >= sites.size())
                throw std::invalid_argument("compiled ED factor site is out of range");
            models::LocalOp local;
            local.name = "compiled";
            local.dim = sites[static_cast<std::size_t>(site)]->dim;
            local.matrix = std::move(matrix);
            term.factors.push_back({site, std::move(local), sites[static_cast<std::size_t>(site)]});
        };
        for (auto const& factor : source.factors) append(factor.site, factor.matrix);
        // Strings act on the input, hence are the rightmost factors.
        for (int site : source.jordan_wigner) {
            if (site < 0 || static_cast<std::size_t>(site) >= sites.size())
                throw std::invalid_argument("compiled ED parity site is out of range");
            std::size_t const dim = static_cast<std::size_t>(sites[static_cast<std::size_t>(site)]->dim);
            models::Matrix parity(dim * dim);
            for (std::size_t state = 0; state < dim; ++state)
                parity[state * dim + state] = spaces.local_spaces[static_cast<std::size_t>(site)].parity[state];
            append(site, std::move(parity));
        }
    }
    models::TermSum const input(std::move(terms));
    models::validate_hermiticity(input, sites, models::FermionConvention::lattice);
    (void)models::compile_hamiltonian(input, sites, spaces.charges, models::FermionConvention::lattice);
}

std::vector<std::optional<double>> charge_targets(ProjectedSpaces const& spaces,
                                                  std::map<std::string, double> const& sector) {
    check_sector(spaces.charges, sector);
    std::vector<std::optional<double>> result;
    for (auto const& name : spaces.charges) {
        auto const target = sector.find(name);
        result.push_back(target == sector.end() ? std::nullopt : std::optional<double>(target->second));
    }
    return result;
}
// Full-model validation and sorted local rows shared by eigenstates and spectra.
// Unselected charges remain mixed, so every declared charge must be conserved.
PreparedModel prepare_model(models::TermSum const& terms, std::vector<models::SiteBasisPtr> const& sites,
    std::vector<std::string> const& conserved, models::FermionConvention convention,
    std::map<std::string, double> const& sector) {
    PreparedModel result;
    result.spaces = projected_spaces(models::physical_basis(sites, conserved, convention));
    models::validate_hermiticity(terms, sites, convention);
    result.op = models::compile_hamiltonian(terms, sites, conserved, convention);
    for (auto const& name : result.spaces.charges) {
        auto const target = sector.find(name);
        result.targets.push_back(target == sector.end() ? std::nullopt : std::optional<double>(target->second));
    }
    return result;
}
std::vector<SpectrumTranslation> lattice_translations(models::Hamiltonian const& model) {
    auto const& lattice = model.lat();
    auto const directions = lattice.periodic_directions();
    if (directions.empty())
        throw std::invalid_argument("momentum resolution requires at least one periodic lattice direction");
    auto generators = lattice.translation_generators();
    std::vector<SpectrumTranslation> translations;
    for (std::size_t i = 0; i < directions.size(); ++i)
        translations.push_back({std::move(generators[i]),
                                static_cast<std::size_t>(lattice.extent()[static_cast<std::size_t>(directions[i])])});
    return translations;
}
std::vector<models::SiteBasisPtr> model_sites(models::Hamiltonian const& model) {
    std::vector<models::SiteBasisPtr> sites;
    for (int site = 0; site < model.num_sites(); ++site) sites.push_back(model.site_for(site));
    return sites;
}
struct SpectrumRequest {
    ModelSpectrum result;
    std::vector<std::size_t> selected;
};
SpectrumRequest spectrum_request(std::vector<std::string> const& conserved,
    std::vector<SpectrumTranslation> const& translations, ModelSpectrumOptions const& options) {
    validate_eigensolver_controls(options.which, options.tolerance, options.max_iterations);
    if (options.count && *options.count <= 0) throw std::invalid_argument("eigenvalue count must be positive");
    check_sector(conserved, options.sector);
    if (options.all_momenta && !options.momentum.empty())
        throw std::invalid_argument("request either one momentum or all momenta");
    if (!translations.empty() && !options.all_momenta && options.momentum.empty())
        throw std::invalid_argument("translations require a momentum or all-momenta request");
    SpectrumRequest request;
    request.result.sector = options.sector;
    for (std::size_t i = 0; i < translations.size(); ++i) request.result.directions.push_back(static_cast<int>(i));
    if (!options.all_momenta) {
        auto momentum = select_momentum(periods_of(translations), options.momentum);
        if (!translations.empty()) request.selected.push_back(momentum.index);
        request.result.momentum = std::move(momentum.momentum);
    }
    return request;
}
ModelSpectrum finish_spectrum(models::CompiledOperator const& op, ProjectedSpaces const& spaces,
    std::vector<SpectrumTranslation> const& translations, std::vector<std::optional<double>> const& targets,
    ModelSpectrumOptions const& options, SpectrumRequest request) {
    request.result.blocks = symmetry_spectrum_from_local(op, spaces.local_spaces, translations, targets,
        options.count.value_or(std::numeric_limits<int>::max()), request.selected, options.which,
        options.tolerance, options.max_iterations);
    request.result.charges = spaces.charges;
    return std::move(request.result);
}
} // namespace

ProjectedProblem prepare_projected(models::CompiledOperator const& op, ProjectedSpaces const& spaces,
    std::vector<SpectrumTranslation> const& translations, ProjectedOptions const& options) {
    auto targets = charge_targets(spaces, options.sector);
    auto momentum = select_momentum(periods_of(translations), options.momentum);
    validate_compiled_model(op, spaces);
    auto basis = make_projected_basis(spaces.local_spaces, translations, targets, momentum.index);
    return ProjectedProblem::prepare(op, basis, spaces.charges);
}

ModelSpectrum model_spectrum(models::CompiledOperator const& op, ProjectedSpaces const& spaces,
    std::vector<SpectrumTranslation> const& translations, ModelSpectrumOptions const& options) {
    auto request = spectrum_request(spaces.charges, translations, options);
    auto targets = charge_targets(spaces, options.sector);
    validate_compiled_model(op, spaces);
    return finish_spectrum(op, spaces, translations, targets, options, std::move(request));
}

ProjectedProblem prepare_projected(models::TermSum const& terms, std::vector<models::SiteBasisPtr> const& sites,
    std::vector<std::string> const& conserved, models::FermionConvention convention,
    std::vector<SpectrumTranslation> const& translations, ProjectedOptions const& options) {
    check_sector(conserved, options.sector);
    auto const momentum = select_momentum(periods_of(translations), options.momentum);
    auto model = prepare_model(terms, sites, conserved, convention, options.sector);
    auto const basis = make_projected_basis(model.spaces.local_spaces, translations, model.targets, momentum.index);
    return ProjectedProblem::prepare(model.op, basis, std::move(model.spaces.charges));
}

ProjectedProblem prepare_projected(models::Hamiltonian const& model, ProjectedOptions const& options) {
    std::vector<SpectrumTranslation> translations;
    if (!options.momentum.empty()) translations = lattice_translations(model);
    return prepare_projected(model.expand(), model_sites(model), model.conserved_qns(), model.fermion_convention(),
                             translations, options);
}

ModelSpectrum model_spectrum(models::TermSum const& terms, std::vector<models::SiteBasisPtr> const& sites,
    std::vector<std::string> const& conserved, models::FermionConvention convention,
    std::vector<SpectrumTranslation> const& translations, ModelSpectrumOptions const& options) {
    auto request = spectrum_request(conserved, translations, options);
    auto model = prepare_model(terms, sites, conserved, convention, options.sector);
    return finish_spectrum(model.op, model.spaces, translations, model.targets, options, std::move(request));
}

ModelSpectrum model_spectrum(models::Hamiltonian const& model, ModelSpectrumOptions const& options) {
    std::vector<SpectrumTranslation> translations;
    if (options.all_momenta || !options.momentum.empty()) translations = lattice_translations(model);
    auto result = model_spectrum(model.expand(), model_sites(model), model.conserved_qns(), model.fermion_convention(),
                                 translations, options);
    if (!translations.empty()) result.directions = model.lat().periodic_directions();
    return result;
}

ProjectedProblem prepare_finite(models::CompiledOperator const& op, ProjectedSpaces const& spaces,
    FiniteGroupIrrep const& irrep, std::map<std::string, double> const& sector) {
    auto targets = charge_targets(spaces, sector);
    validate_compiled_model(op, spaces);
    return ProjectedProblem::prepare(op, make_finite_group_basis(spaces.local_spaces, irrep, targets), spaces.charges);
}
ProjectedProblem prepare_finite(models::TermSum const& terms, std::vector<models::SiteBasisPtr> const& sites,
    std::vector<std::string> const& conserved, models::FermionConvention convention, FiniteGroupIrrep const& irrep,
    std::map<std::string, double> const& sector) {
    auto model = prepare_model(terms, sites, conserved, convention, sector);
    return ProjectedProblem::prepare(model.op, make_finite_group_basis(model.spaces.local_spaces, irrep, model.targets),
                                     std::move(model.spaces.charges));
}
ProjectedProblem prepare_finite(models::Hamiltonian const& model, FiniteGroupIrrep const& irrep,
                                std::map<std::string, double> const& sector) {
    return prepare_finite(model.expand(), model_sites(model), model.conserved_qns(), model.fermion_convention(), irrep,
                          sector);
}

namespace {
void finite_controls(FiniteSpectrumOptions const& options) {
    validate_eigensolver_controls(options.which, options.tolerance, options.max_iterations);
    if (options.count && *options.count <= 0) throw std::invalid_argument("eigenvalue count must be positive");
}
ModelSpectrum finite_result(models::CompiledOperator const& op, ProjectedSpaces const& spaces,
    std::vector<FiniteGroupIrrep> const& irreps, std::vector<std::optional<double>> const& targets,
    FiniteSpectrumOptions const& options) {
    ModelSpectrum result;
    result.charges = spaces.charges;
    result.sector = options.sector;
    result.blocks = finite_spectrum_from_local(op, spaces.local_spaces, irreps, targets,
        options.count.value_or(std::numeric_limits<int>::max()), options.which, options.tolerance,
        options.max_iterations);
    return result;
}
} // namespace

ModelSpectrum finite_spectrum(models::CompiledOperator const& op, ProjectedSpaces const& spaces,
    std::vector<FiniteGroupIrrep> const& irreps, FiniteSpectrumOptions const& options) {
    finite_controls(options);
    auto targets = charge_targets(spaces, options.sector);
    validate_compiled_model(op, spaces);
    return finite_result(op, spaces, irreps, targets, options);
}
ModelSpectrum finite_spectrum(models::TermSum const& terms, std::vector<models::SiteBasisPtr> const& sites,
    std::vector<std::string> const& conserved, models::FermionConvention convention,
    std::vector<FiniteGroupIrrep> const& irreps, FiniteSpectrumOptions const& options) {
    finite_controls(options);
    auto model = prepare_model(terms, sites, conserved, convention, options.sector);
    return finite_result(model.op, model.spaces, irreps, model.targets, options);
}
ModelSpectrum finite_spectrum(models::Hamiltonian const& model, std::vector<FiniteGroupIrrep> const& irreps,
                              FiniteSpectrumOptions const& options) {
    return finite_spectrum(model.expand(), model_sites(model), model.conserved_qns(), model.fermion_convention(), irreps,
                           options);
}

} // namespace alps::ed
