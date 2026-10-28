/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/make_parameters_from_file.hpp>
#include <alps/parameter.h>
#include <alps/parameter/parameters_p.h>

#include "evaluated_parameters.hpp"

#include <boost/filesystem/fstream.hpp>
#include <boost/spirit/include/classic_core.hpp>

#include <sstream>
#include <stdexcept>
#include <string>

namespace alps {

    namespace {
        // `KEY = value` lines, then one { } block and nothing after it, with the
        // grammar and comments of ParameterList; the block's keys win
        bool parse_one_block(std::string const & text, Parameters & result) {
            namespace bs = boost::spirit::classic;
            Parameters global, local;
            ParametersParser global_p(global), local_p(local);
            bs::parse_info<std::string::const_iterator> info = bs::parse(
                text.begin(), text.end(),
                *bs::eol_p >> !global_p
                    >> bs::ch_p('{') >> *bs::eol_p >> !local_p >> bs::ch_p('}') >> *bs::eol_p
                    >> bs::end_p,
                bs::blank_p | bs::comment_p("//") | bs::comment_p("/*", "*/"));
            if (!info.full)
                return false;
            result = global;
            for (Parameters::const_iterator it = local.begin(); it != local.end(); ++it)
                result[it->key()] = it->value();
            result.replace_envvar();
            return true;
        }
    }

    params make_parameters_from_file(boost::filesystem::path const & path) {
        boost::filesystem::ifstream ifs(path);
        std::ostringstream contents;
        contents << ifs.rdbuf();

        std::istringstream text(contents.str());
        try {
            return detail::make_parameters_evaluated(Parameters(text));
        } catch (std::runtime_error const & error) {
            Parameters block;
            if (parse_one_block(contents.str(), block))
                return detail::make_parameters_evaluated(block);

            // ParameterList holds one set per { } block, each with the keys set before it
            ParameterList sets;
            try {
                std::istringstream blocks(contents.str());
                sets.parse(blocks);
            } catch (std::runtime_error const &) {
                throw error;  // not blocks either
            }
            if (sets.size() > 1)
                throw std::runtime_error(path.string() + " has " + std::to_string(sets.size())
                    + " parameter sets in { } blocks, and a program reads one. "
                      "parameter2xml writes one XML job file per set.");
            if (sets.size() == 1)
                throw std::runtime_error(path.string() + " has text after its { } block "
                    "or a #clear or #stop, which a program does not read. "
                    "Set keys before the block or inside it.");
            throw error;
        }
    }

}
