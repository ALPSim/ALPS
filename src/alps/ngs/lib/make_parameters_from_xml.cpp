/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>

#include <string>

namespace alps {

    ALPS_DECL params make_parameters_from_xml(boost::filesystem::path const & arg) {
        Parameters par;
        boost::filesystem::ifstream infile(arg.string());

        // read outermost tag (e.g. <SIMULATION>)
        XMLTag tag = parse_tag(infile, true);
        std::string closingtag = "/" + tag.name;

        // scan for <PARAMETERS> and read them
        tag = parse_tag(infile, true);
        while (tag.name != "PARAMETERS" && tag.name != closingtag) {
            std::cerr << "skipping tag with name " << tag.name << "\n";
            skip_element(infile, tag);
            tag = parse_tag(infile, true);
        }

        par.read_xml(tag, infile, true);
        if (!par.defined("SEED"))
            par["SEED"] = 0;
        
        params res;
        for (Parameters::const_iterator it = par.begin(); it != par.end(); ++it)
            res[it->key()] = it->value();
        return res;
    }
}
