/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#if !defined(BASEPARAMETERS_H) && !defined(DMRGPARAMETERS_H)
#define BASEPARAMETERS_H

#include <string>
#include <fstream>
#include <iostream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>

#include <alps/parameter.h>
#include "utils/io.hpp"

#include "dmrg/utils/parameter_proxy.h"

                
namespace parameters {
    class value {
    public:
        value () : val_(""), empty_(true) { }
        
        template <class T>
        value (const T & val)
        : val_(boost::lexical_cast<std::string>(val))
        , empty_(false)
        { }
        
        std::string get() const {return val_;}
        bool empty() const {return empty_;}
    private:
        std::string val_;
        bool empty_;
        
    };
}
                
class BaseParameters : public alps::Parameters
{
public:
    
    BaseParameters ()
    : alps::Parameters()
    { }
    
    BaseParameters (std::istream& param_file)
    : alps::Parameters()
    {
        try {
            alps::Parameters temp(param_file);
            *static_cast<alps::Parameters*>(this) = temp;
        } catch (std::exception & e) {
            maquis::cerr << "Exception thrown when parsing parameters:" << std::endl;
            maquis::cerr << e.what() << std::endl;
            exit(1);
        }
    }

    BaseParameters(alps::Parameters const& p)
    : alps::Parameters(p)
    { }

    bool is_set (std::string const & key) const
    {
        return defined(key);
    }

    parameters::proxy operator[](std::string const& key)
    {
        if (!defined(key)) {
            std::map<std::string, std::string>::const_iterator match = defaults.find(key);
            if (match != defaults.end())
                alps::Parameters::operator[](key) = match->second;
            else
                boost::throw_exception(std::runtime_error("parameter " + key + " not defined"));
        }
        return parameters::proxy(alps::Parameters::operator[](key));
    }
    
    template<class T> T get(std::string const & key)
    {
        parameters::proxy const& p = (*this)[key];
        return p.as<T>();
    }
    
    template<class T> void set(std::string const & key, T const & value)
    {
        alps::Parameters::operator[](key) = boost::lexical_cast<std::string>(value);
    }
    
    BaseParameters iteration_params(std::string const & var, std::size_t val)
    {
        BaseParameters p;
        
        boost::regex expression("^(.*)\\[" + var + "\\]$");
        boost::smatch what;
        for (alps::Parameters::const_iterator it=this->begin();it != this->end();++it) {
            std::string key = it->key();
            if (boost::regex_match(key, what, expression)) {
                std::vector<value_type> v = (*this)[key];
                if (val < v.size())
                    p.set(what.str(1), v[val]);
                else
                    p.set(what.str(1), *(v.rbegin()));
            }
        }
        p.set(var, val);
        return p;
    }
    
    BaseParameters & operator<<(BaseParameters const& p)
    {
        for (alps::Parameters::const_iterator it=p.begin(); it!=p.end(); ++it)
            alps::Parameters::operator[](it->key()) = it->value();
        defaults.insert(p.defaults.begin(), p.defaults.end());

        return *this;
    }
    
protected:
    void add_option(std::string const & name,
                    std::string const & desc,
                    parameters::value const & val = parameters::value())
    {
        if (!val.empty())
            defaults[name] = val.get();
        descriptions[name] = desc;        
    }
    
    std::map<std::string, std::string> defaults;
    std::map<std::string, std::string> descriptions;
    
};

#endif
