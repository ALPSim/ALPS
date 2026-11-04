// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/input/looper.hpp>

#include <alps/expression.h>
#include <boost/lexical_cast.hpp>

namespace alps { namespace input {

    looper_temperature_parameters read_looper_temperature_parameters(Parameters const & p) {
        looper_temperature_parameters result;
        if (p.defined("T")) {
            result.temperature = alps::evaluate("T", p);
        } else if (p.defined("TEMPERATURE")) {
            result.temperature = alps::evaluate("TEMPERATURE", p);
        } else if (p.defined("BETA")) {
            result.temperature = 1 / alps::evaluate("BETA", p);
        } else if (p.defined("INVERSE_TEMPERATURE")) {
            result.temperature = 1 / alps::evaluate("INVERSE_TEMPERATURE", p);
        } else {
            return result;
        }
        for (int i = 0;; ++i) {
            int td = -1;
            double ts = -1;
            std::string ns = boost::lexical_cast<std::string>(i);
            if (p.defined("T_DURATION_" + ns) && p.defined("T_START_" + ns)) {
                td = static_cast<int>(alps::evaluate("T_DURATION_" + ns, p));
                ts = alps::evaluate("T_START_" + ns, p);
            } else if (p.defined("TEMPERATURE_DURATION_" + ns) && p.defined("TEMPERATURE_START_" + ns)) {
                td = static_cast<int>(alps::evaluate("TEMPERATURE_DURATION_" + ns, p));
                ts = alps::evaluate("TEMPERATURE_START_" + ns, p);
            } else if (p.defined("BETA_DURATION_" + ns) && p.defined("BETA_START_" + ns)) {
                td = static_cast<int>(alps::evaluate("BETA_DURATION_" + ns, p));
                ts = 1 / alps::evaluate("BETA_START_" + ns, p);
            } else if (p.defined("INVERSE_TEMPERATURE_DURATION_" + ns) &&
                       p.defined("INVERSE_TEMPERATURE_START_" + ns)) {
                td = static_cast<int>(alps::evaluate("INVERSE_TEMPERATURE_DURATION_" + ns, p));
                ts = 1 / alps::evaluate("INVERSE_TEMPERATURE_START_" + ns, p);
            }
            if (td <= 0 || ts <= 0) break;
            result.annealing.push_back(std::make_pair(td, ts));
        }
        return result;
    }

    looper_worker_parameters read_looper_worker_parameters(Parameters const & p) {
        looper_worker_parameters result;
        result.disable_improved_estimator = p.defined("DISABLE_IMPROVED_ESTIMATOR");
        return result;
    }

    looper_correlation_parameters read_looper_correlation_parameters(Parameters const & p) {
        looper_correlation_parameters result;
        result.measure_correlations = p.value_or_default("MEASURE[Correlations]", false);
        result.measure_green_function = p.value_or_default("MEASURE[Green Function]", false);
        result.measure_structure_factor = p.value_or_default("MEASURE[Structure Factor]", false);
        return result;
    }

    std::optional<int> read_looper_initial_site(Parameters const & p) {
        if (!p.defined("INITIAL_SITE"))
            return std::nullopt;
        return static_cast<int>(p["INITIAL_SITE"]);
    }

    looper_weight_parameters read_looper_weight_parameters(Parameters const & p) {
        looper_weight_parameters result;
        if (p.defined("FORCE_SCATTER"))
            result.force_scatter = static_cast<double>(p["FORCE_SCATTER"]);
        result.scattering_ratio = p.value_or_default("SCATTERING_RATIO", 0.5);
        return result;
    }

    looper_model_parameters read_looper_model_parameters(Parameters const & p) {
        looper_model_parameters result;
        result.disorder_seed = p.value_or_default("DISORDER_SEED", 0);
        result.use_site_indices_as_types = p.value_or_default("USE_SITE_INDICES_AS_TYPES", false);
        result.use_bond_indices_as_types = p.value_or_default("USE_BOND_INDICES_AS_TYPES", false);
        if (p.defined("LOOPER_DEBUG[MODEL OUTPUT]"))
            result.debug_model_output = static_cast<std::string>(p["LOOPER_DEBUG[MODEL OUTPUT]"]);
        return result;
    }

    looper_wanglandau_parameters read_looper_wanglandau_parameters(Parameters const & p) {
        looper_wanglandau_parameters result;
        result.multicanonical = !p.defined("DISABLE_MULTICANONICAL_MEASUREMENT");
        result.zhou_bhatt = !p.defined("DISABLE_ZHOU_BHATT");
        result.block_sweeps = result.zhou_bhatt ? 0 : static_cast<int>(p.value_or_default("BLOCK_SWEEPS", 65536));
        result.sweeps = p.value_or_default("SWEEPS", 65536);
        result.iterations = p.value_or_default("WANG_LANDAU_ITERATIONS", 16);
        return result;
    }

} }
