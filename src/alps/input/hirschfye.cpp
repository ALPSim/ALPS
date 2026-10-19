// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/input/hirschfye.hpp>

namespace alps { namespace input {

    namespace {

        std::optional<std::string> optional_string(Parameters const & parms, std::string const & name) {
            if (!parms.defined(name))
                return std::nullopt;
            return static_cast<std::string>(parms[name]);
        }

    }

    hirschfye_parameters read_hirschfye_parameters(Parameters const & parms) {
        hirschfye_parameters result;
        result.thermalization = static_cast<int>(parms["THERMALIZATION"]);
        result.sweeps = static_cast<int>(parms["SWEEPS"]);
        result.beta = static_cast<double>(parms["BETA"]);
        result.n_tau = static_cast<int>(parms["N"]);
        result.sites = static_cast<int>(parms.value_or_default("SITES", 1));
        result.U = static_cast<double>(parms["U"]);
        result.fourpoint_interval = static_cast<int>(parms.value_or_default("FOURPOINT_INTERVAL", 1000));
        result.measure_fourpoint = static_cast<bool>(parms.value_or_default("MEASURE_FOURPOINT_FUNCTION", false));
        result.tolerance = static_cast<double>(parms["TOLERANCE"]);
        result.max_time = static_cast<int>(parms.value_or_default("MAX_TIME", 86400));
        result.n_matsubara = static_cast<unsigned int>(parms["NMATSUBARA"]);
        return result;
    }

    hirschfye_result_parameters read_hirschfye_result_parameters(Parameters const & parms) {
        hirschfye_result_parameters result;
        result.n_tau = parms["N"];
        result.n_matsubara = parms["NMATSUBARA"];
        result.sites = parms["SITES"];
        result.checkpoint = optional_string(parms, "CHECKPOINT");
        result.outfile = optional_string(parms, "OUTFILE");
        return result;
    }

    hirschfye_input_sizes read_hirschfye_input_sizes(Parameters const & parms) {
        hirschfye_input_sizes result;
        result.n_matsubara = static_cast<int>(parms["NMATSUBARA"]);
        result.sites = parms.value_or_default("SITES", 1);
        result.flavors = parms.value_or_default("FLAVORS", 2);
        return result;
    }

} }
