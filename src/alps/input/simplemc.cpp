// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/input/simplemc.hpp>

#include <alps/expression.h>
#include <boost/lexical_cast.hpp>

namespace alps { namespace input {

    simplemc_parameters read_simplemc_parameters(Parameters const & p, int num_bond_types) {
        simplemc_parameters result;
        result.beta = (p.defined("T") ? 1 / evaluate("T", p) : 0);
        result.field = (p.defined("H") ? evaluate("H", p) : 0);
        result.couplings.resize(num_bond_types, p.defined("J") ? evaluate("J", p) : 1);
        for (int t = 0; t < num_bond_types; ++t) {
            std::string name = std::string("J") + boost::lexical_cast<std::string>(t);
            if (p.defined(name)) result.couplings[t] = evaluate(name, p);
        }
        result.snapshot_interval =
            static_cast<int>(p.defined("SNAPSHOT_INTERVAL") ? evaluate("SNAPSHOT_INTERVAL", p) : 0);
        result.base_name = static_cast<std::string>(p["BASE_NAME"]);
        result.clone_id = static_cast<int>(p["CLONE_ID"]);
        return result;
    }

} }
