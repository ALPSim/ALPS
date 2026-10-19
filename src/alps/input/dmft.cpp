// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/input/dmft.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>

#include <stdexcept>

namespace alps { namespace input {

    namespace {

        std::optional<std::string> optional_string(Parameters const & parms, std::string const & name) {
            if (!parms.defined(name))
                return std::nullopt;
            return static_cast<std::string>(parms[name]);
        }

        std::optional<std::string> nonempty_string(Parameters const & parms, std::string const & name) {
            std::optional<std::string> value = optional_string(parms, name);
            if (value && value->empty())
                return std::nullopt;
            return value;
        }

        std::string flavor_name(std::string const & prefix, int f) {
            return prefix + boost::lexical_cast<std::string>(f);
        }

    }

    dmft_sizes read_dmft_sizes(Parameters const & parms) {
        dmft_sizes result;
        result.flavors = parms.value_or_default("FLAVORS", 2);
        result.sites = parms.value_or_default("SITES", 1);
        return result;
    }

    dmft_run_parameters read_dmft_run_parameters(Parameters const & parms) {
        dmft_run_parameters result;
        result.cluster_loop = parms.defined("CLUSTER_LOOP");
        result.omega_loop = parms.value_or_default("OMEGA_LOOP", false);
        if (!result.cluster_loop)
            result.solver = static_cast<std::string>(parms["SOLVER"]);
        result.sc_write_delta = parms.value_or_default("SC_WRITE_DELTA", false);
        result.antiferromagnet = parms.value_or_default("ANTIFERROMAGNET", false);
        result.symmetrization = parms.value_or_default("SYMMETRIZATION", false);
        result.semicircle_hilbert = parms.value_or_default("SEMICIRCLE_HILBERT", false);
        result.sizes = read_dmft_sizes(parms);
        return result;
    }

    dmft_loop_parameters read_dmft_loop_parameters(Parameters const & parms, dmft_loop loop) {
        bool const matsubara = loop == dmft_loop::matsubara;
        dmft_loop_parameters result;
        result.n_tau = static_cast<int>(parms["N"]);
        result.n_matsubara = matsubara ? static_cast<int>(parms["NMATSUBARA"]) : 0;
        result.sizes.flavors = parms.value_or_default("FLAVORS", 2);
        result.sizes.sites = matsubara ? static_cast<int>(parms.value_or_default("SITES", 1)) : 1;
        result.beta = static_cast<double>(parms["BETA"]);
        result.h = parms.value_or_default("H", 0.);
        result.mu = loop != dmft_loop::dressed_imaginary_time ? static_cast<double>(parms["MU"]) : 0.;
        result.converged = static_cast<double>(parms["CONVERGED"]);
        result.symmetrization = static_cast<bool>(parms["SYMMETRIZATION"]);
        result.relax_rate = matsubara ? static_cast<double>(parms.value_or_default("RELAX_RATE", 1.)) : 1.;
        result.max_it = static_cast<int>(parms.value_or_default("MAX_IT", 1000));
        result.basename = static_cast<std::string>(parms["BASENAME"]);
        result.h_init = optional_string(parms, "H_INIT");
        result.g0_tau_output = static_cast<std::string>(parms.value_or_default("G0TAU_output", "G0tau_output"));
        result.g0_omega_output = static_cast<std::string>(parms.value_or_default("G0OMEGA_output", "G0omega_output"));
        if (matsubara)
            result.initial_g0_tau_output = optional_string(parms, "G0TAU_input");
        return result;
    }

    dmft_initial_g0_parameters read_dmft_initial_g0_tau(Parameters const & parms) {
        dmft_initial_g0_parameters result;
        result.n_tau = boost::lexical_cast<int>(parms["N"]);
        result.n_matsubara = 0;
        result.flavors = parms.value_or_default("FLAVORS", 2);
        result.beta = result.mu = result.h = 0.;
        result.insulating = false;
        result.input = nonempty_string(parms, "G0TAU_INPUT");
        result.output = optional_string(parms, "G0TAU_input");
        return result;
    }

    dmft_initial_g0_parameters read_dmft_initial_g0_omega(Parameters const & parms) {
        dmft_initial_g0_parameters result;
        result.n_matsubara = boost::lexical_cast<unsigned int>(parms["NMATSUBARA"]);
        result.n_tau = boost::lexical_cast<unsigned int>(parms["N"]);
        result.flavors = parms.value_or_default("FLAVORS", 2);
        result.beta = static_cast<double>(parms["BETA"]);
        result.mu = static_cast<double>(parms["MU"]);
        result.h = static_cast<double>(parms.value_or_default("H", 0.));
        result.insulating = parms.defined("INSULATING");
        result.input = nonempty_string(parms, "G0OMEGA_INPUT");
        result.output = optional_string(parms, "G0OMEGA_input");
        return result;
    }

    bool read_dmft_antiferromagnet(Parameters const & parms) {
        return parms.value_or_default("ANTIFERROMAGNET", false);
    }

    dmft_bandstructure_choice read_dmft_bandstructure(Parameters const & parms) {
        dmft_bandstructure_choice result;
        if (parms.defined("DOSFILE")) {
            result.kind = dmft_bandstructure::dos;
            result.dosfile = static_cast<std::string>(parms["DOSFILE"]);
        } else if (parms.defined("TWODBS") && parms["TWODBS"] == "hexagonal")
            result.kind = dmft_bandstructure::hexagonal;
        else if (parms.defined("TWODBS"))
            result.kind = dmft_bandstructure::square;
        else
            result.kind = dmft_bandstructure::semicircle;
        return result;
    }

    std::vector<double> read_semicircle_hoppings(Parameters const & parms) {
        unsigned n_flavor = parms.value_or_default("FLAVORS", 2);
        if (n_flavor % 2 != 0)
            boost::throw_exception(std::logic_error("SemicircleBandstructure: current implementation does not allow for odd n_flavor."));
        std::vector<double> result;
        for (unsigned int f = 0; f < n_flavor / 2; ++f) {
            // flavors (2f) and (2f+1) have the same bandwidth
            std::string const t_f = flavor_name("t", f);
            result.push_back(parms.defined(t_f) ? static_cast<double>(parms[t_f]) : static_cast<double>(parms["t"]));
        }
        return result;
    }

    dmft_lattice_hoppings read_lattice_hoppings(Parameters const & parms) {
        dmft_lattice_hoppings result;
        result.t = static_cast<double>(parms.value_or_default("t", 1.));
        result.tprime = static_cast<double>(parms.value_or_default("tprime", 0.));
        result.L = static_cast<int>(parms.value_or_default("L", 128));
        return result;
    }

    dmft_band_moments read_band_moments(Parameters const & parms, int flavors) {
        dmft_band_moments result;
        for (int f = 0; f < flavors; ++f) {
            std::string const eps = flavor_name("EPS_", f);
            std::string const epssq = flavor_name("EPSSQ_", f);
            result.eps.push_back(parms.defined(eps) ? std::optional<double>(static_cast<double>(parms[eps])) : std::nullopt);
            result.epssq.push_back(parms.defined(epssq) ? std::optional<double>(static_cast<double>(parms[epssq])) : std::nullopt);
        }
        return result;
    }

    dmft_fourier_parameters read_dmft_fourier_parameters(Parameters const & parms, bool with_mu, bool with_U) {
        dmft_fourier_parameters result;
        result.sizes = read_dmft_sizes(parms);
        result.h = static_cast<double>(parms.value_or_default("H", 0.));
        result.U = with_U ? static_cast<double>(parms["U"]) : 0.;
        dmft_band_moments const moments = read_band_moments(parms, result.sizes.flavors);
        for (int f = 0; f < result.sizes.flavors; ++f) {
            result.eps.push_back(moments.eps[f].value_or(0.0));
            result.epssq.push_back(moments.epssq[f].value_or(1.0));
        }
        result.beta = static_cast<double>(parms["BETA"]);
        result.mu = with_mu ? static_cast<double>(parms["MU"]) : 0.;
        return result;
    }

    dmft_interaction_parameters read_dmft_interaction(Parameters const & parms) {
        dmft_interaction_parameters result;
        result.sizes = read_dmft_sizes(parms);
        result.U = result.J = result.Uprime = 0.;
        result.u_matrix = optional_string(parms, "U_MATRIX");
        if (result.u_matrix)
            return result;
        result.U = static_cast<double>(parms["U"]);
        if (result.sizes.flavors != 1 && result.sizes.flavors != 2) {
            result.J = static_cast<double>(parms["J"]);
            result.Uprime = parms.value_or_default("U'", result.U - 2 * result.J);
        }
        return result;
    }

    dmft_external_solver_parameters read_dmft_external_solver(Parameters const & parms, dmft_loop loop) {
        bool const matsubara = loop == dmft_loop::matsubara;
        dmft_external_solver_parameters result;
        result.n_matsubara = matsubara ? static_cast<int>(static_cast<unsigned int>(parms["NMATSUBARA"])) : 0;
        result.n_tau = static_cast<unsigned int>(parms["N"]);
        result.sizes.sites = static_cast<unsigned int>(parms.value_or_default("SITES", 1));
        result.sizes.flavors = static_cast<unsigned int>(parms.value_or_default("FLAVORS", 2));
        result.tmpname = optional_string(parms, "TMPNAME");
        result.sc_write_delta = parms.value_or_default("SC_WRITE_DELTA", false);
        result.antiferromagnet = false;
        result.beta = result.mu = result.h = 0.;
        if (result.sc_write_delta) {
            if (!matsubara)
                result.antiferromagnet = parms.value_or_default("ANTIFERROMAGNET", false);
            result.beta = static_cast<double>(parms["BETA"]);
            if (matsubara) {
                result.mu = static_cast<double>(parms["MU"]);
                result.h = static_cast<double>(parms.value_or_default("H", 0.));
            }
        }
        return result;
    }

    dmft_hybridization_handoff read_dmft_hybridization_handoff(Parameters const & parms) {
        dmft_hybridization_handoff result;
        result.mu = static_cast<double>(parms["MU"]);
        result.h = static_cast<double>(parms.value_or_default("H", 0.));
        result.flavors = static_cast<unsigned>(parms["FLAVORS"]);
        return result;
    }

    unsigned int read_dmft_max_time(Parameters const & parms, unsigned int default_time) {
        return static_cast<unsigned int>(parms.value_or_default("MAX_TIME", default_time));
    }

} }
