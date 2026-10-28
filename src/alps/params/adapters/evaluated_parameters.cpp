// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include "evaluated_parameters.hpp"

#include <alps/expression.h>
#include <alps/params/text_values.hpp>

#include <cstdlib>
#include <exception>
#include <string>

namespace alps {
    namespace detail {

        namespace {
            // Text that strtod reads to its end, such as `10`, `2.25` or `1e-3`
            bool is_number(std::string const & text) {
                char const * begin = text.c_str();
                char * end = nullptr;
                std::strtod(begin, &end);
                return end != begin && *end == '\0';
            }
        }

        params make_parameters_evaluated(Parameters const & parameters) {
            params res;
            for (Parameters::const_iterator it = parameters.begin(); it != parameters.end(); ++it) {
                std::string text = it->value();
                if (!text.empty() && !is_number(text)) {
                    try {
                        if (can_evaluate(text, parameters))
                            text = params_ns::to_text(evaluate<double>(text, parameters));
                    } catch (std::exception const &) {
                        // The expression parser rejects it, as it does `chain lattice`
                    }
                }
                res.supply(it->key(), text);
            }
            return res;
        }

    }
}
