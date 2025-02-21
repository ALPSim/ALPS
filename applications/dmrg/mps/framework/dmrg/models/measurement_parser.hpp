/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2016 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef MODELS_MEASUREMENT_PARSER_HPP
#define MODELS_MEASUREMENT_PARSER_HPP

#include "dmrg/utils/BaseParameters.h"

#include <string>
#include <vector>
#include <boost/regex.hpp>

struct MeasurementSpecs {
    std::string meas_class;
    std::string name;
    std::string args;

    MeasurementSpecs() { }
    MeasurementSpecs(std::string const& class_, std::string const& name_, std::string const& args_)
    : meas_class(class_)
    , name(name_)
    , args(args_)
    { }
};


inline std::vector<MeasurementSpecs> parse_measurements(BaseParameters const& parms)
{
    
    std::vector<MeasurementSpecs> meas_specs;
    
    boost::regex measurement_expr("MEASURE_([A-Z_]+)[[:space:]]*\[[[:space:]]*(.*)[[:space:]]*]");
    
    for( BaseParameters::const_iterator it_p = parms.begin() ; it_p != parms.end() ; ++it_p)
    {
        std::string param_name = it_p->key();
        
        if( param_name.find("MEASURE_") == std::string::npos )
            continue;
        
        boost::smatch what;
        if (boost::regex_search(param_name, what, measurement_expr)) {
            std::string measurement_type( what.str(1) );
            std::string measurement_name( what.str(2) );
            
            meas_specs.push_back( MeasurementSpecs( measurement_type , measurement_name , it_p->value() ) );
        }
    }
    
    return meas_specs;
}

#endif