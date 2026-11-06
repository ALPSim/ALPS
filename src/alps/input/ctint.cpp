// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/input/ctint.hpp>

#include <vector>

namespace alps { namespace input {

    namespace {

        // Defines a parameter that has an older or a derived alternative only
        // if the input supplies it.
        template <typename T>
        void define_if_supplied(params & p, std::string const & name, std::string const & description) {
            if (p.supplied(name))
                p.define<T>(name, description);
        }

    }

    void define_ctint_parameters(params & p) {
        p.define<long>("SWEEPS", "number of measurement steps");
        p.define<long>("THERMALIZATION", "number of thermalization steps");
        p.define<long>("MAX_TIME", "time limit in seconds");
        p.define<double>("BETA", "inverse temperature");
        p.define<double>("U", "on-site interaction");
        p.define<double>("ALPHA", "shift of the auxiliary field");
        p.define<int>("MAX_ORDER", 2048, "maximum perturbation order");
        p.define<int>("SITES", 1, "number of sites");
        p.define<int>("N_ORBITALS", 2, "number of orbitals; the number of flavors unless FLAVORS is set");
        p.define<int>("N_MEAS", 200, "number of updates between measurements unless MEASUREMENT_PERIOD is set");
        p.define<int>("RECALC_PERIOD", 5000, "number of updates between recomputations of the inverse Green's function matrix");
        p.define<double>("ALMOSTZERO", 1.e-16, "roundoff threshold");
        p.define<bool>("HISTOGRAM_MEASUREMENT", false, "measure the self-energy in imaginary time with histograms instead of on Matsubara frequencies");
        p.define<bool>("ATOMIC", false, "start from the atomic limit instead of the bare Green's function in INFILE");
        p.define<bool>("DMFT_FRAMEWORK", false, "read INFILE in the layout of the DMFT framework");
        p.define<int>("TAU_DISCRETIZATION_FOR_EXP", 0, "number of imaginary-time points at which exp(i omega_n tau) is tabulated; 0 computes it on the fly");
        if (!p["ATOMIC"].as<bool>()) p.define<std::string>("INFILE", "HDF5 file with the bare Green's function G0");
        else define_if_supplied<std::string>(p, "INFILE", "HDF5 file with the bare Green's function G0");
        // sizes with older names (which take precedence) or derived defaults
        define_if_supplied<int>(p, "FLAVORS", "number of flavors; N_ORBITALS by default");
        define_if_supplied<int>(p, "NMATSUBARA", "older name of N_MATSUBARA");
        if (p.supplied("N_MATSUBARA") || !p.supplied("NMATSUBARA")) p.define<int>("N_MATSUBARA", "number of Matsubara frequencies");
        define_if_supplied<int>(p, "NMATSUBARA_MEASUREMENTS", "number of Matsubara frequencies measured; N_MATSUBARA by default");
        define_if_supplied<int>(p, "N", "older name of N_TAU");
        if (p.supplied("N_TAU") || !p.supplied("N")) p.define<int>("N_TAU", "number of imaginary-time slices");
        define_if_supplied<int>(p, "NSELF", "number of imaginary-time points of the self-energy; 10*N_TAU by default");
        define_if_supplied<int>(p, "MEASUREMENT_PERIOD", "number of updates between measurements; N_MEAS by default");
        define_if_supplied<int>(p, "CONVERGENCE_CHECK_PERIOD", "number of updates between convergence checks; RECALC_PERIOD by default");
        // read by the Fourier transformer and the interaction matrix through alps::Parameters
        define_if_supplied<double>(p, "MU", "chemical potential");
        define_if_supplied<double>(p, "H", "magnetic field");
        define_if_supplied<double>(p, "J", "Hund's coupling");
        define_if_supplied<double>(p, "U'", "interaction between different orbitals");
        define_if_supplied<std::string>(p, "U_MATRIX", "interaction matrix");
        std::vector<std::string> flavor_energies;
        for (auto const & supplied : p.supplied_text())
            if (supplied.first.compare(0, 4, "EPS_") == 0 || supplied.first.compare(0, 6, "EPSSQ_") == 0)
                flavor_energies.push_back(supplied.first);
        for (auto const & name : flavor_energies)
            p.define<double>(name, "EPS_f: energy of flavor f; EPSSQ_f: its second moment");
    }

    ctint_parameters read_ctint_parameters(params const & p) {
        ctint_parameters result;
        result.max_order = static_cast<int>(p["MAX_ORDER"]);
        result.flavors = p.defined("FLAVORS") ? static_cast<int>(p["FLAVORS"]) : static_cast<int>(p["N_ORBITALS"]);
        result.sites = static_cast<int>(p["SITES"]);
        result.n_matsubara = p.defined("NMATSUBARA") ? static_cast<int>(p["NMATSUBARA"]) : static_cast<int>(p["N_MATSUBARA"]);
        result.n_matsubara_measurements = p.defined("NMATSUBARA_MEASUREMENTS")
            ? static_cast<int>(p["NMATSUBARA_MEASUREMENTS"]) : static_cast<int>(result.n_matsubara);
        result.n_tau = p.defined("N") ? static_cast<int>(p["N"]) : static_cast<int>(p["N_TAU"]);
        result.n_self = p.defined("NSELF") ? static_cast<int>(p["NSELF"]) : static_cast<int>(10 * result.n_tau);
        result.sweeps = static_cast<std::uint64_t>(p["SWEEPS"].as<long>());
        result.thermalization = static_cast<unsigned int>(p["THERMALIZATION"].as<long>());
        result.max_time = p["MAX_TIME"].as<long>();
        result.beta = p["BETA"];
        result.U = p["U"];
        result.alpha = p["ALPHA"];
        result.recalc_period = static_cast<int>(p["RECALC_PERIOD"]);
        result.measurement_period = p.defined("MEASUREMENT_PERIOD")
            ? static_cast<int>(p["MEASUREMENT_PERIOD"]) : static_cast<int>(p["N_MEAS"]);
        result.convergence_check_period = p.defined("CONVERGENCE_CHECK_PERIOD")
            ? static_cast<int>(p["CONVERGENCE_CHECK_PERIOD"]) : static_cast<int>(p["RECALC_PERIOD"]);
        result.almost_zero = p["ALMOSTZERO"];
        result.histogram_measurement = p["HISTOGRAM_MEASUREMENT"].as<bool>();
        result.atomic = p["ATOMIC"].as<bool>();
        result.dmft_framework = p["DMFT_FRAMEWORK"].as<bool>();
        result.tau_discretization_for_exp = static_cast<int>(p["TAU_DISCRETIZATION_FOR_EXP"]);
        result.infile = p.defined("INFILE") ? p["INFILE"].as<std::string>() : std::string();
        return result;
    }

} }
