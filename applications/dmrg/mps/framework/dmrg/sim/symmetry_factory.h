/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/block_matrix/symmetry.h"

#include <boost/function.hpp>
#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>
#include <map>
#include <string>

#include "utils/io.hpp"

std::string guess_alps_symmetry(BaseParameters & parms)
{
    std::map<int, std::string> symm_names;
    symm_names[0] = "none";
    symm_names[1] = "u1";
    symm_names[2] = "2u1";
    
    int n=0;
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    if (parms.defined("CONSERVED_QUANTUMNUMBERS")) {
        boost::char_separator<char> sep(" ,");
        std::string qn_string = parms["CONSERVED_QUANTUMNUMBERS"].str();
        tokenizer qn_tokens(qn_string, sep);
        for (tokenizer::iterator it=qn_tokens.begin(); it != qn_tokens.end(); it++) {
            if (parms.defined(*it + "_total"))
                n += 1;
        }
    }
    return symm_names[n];
}

namespace dmrg {
    
    template <class TR>
    typename TR::shared_ptr symmetry_factory(DmrgParameters & parms)
    {
        typedef typename TR::shared_ptr ptr_type;
        std::map<std::string, ptr_type> factory_map;
        
        maquis::cout << "This binary contains symmetries: ";
#ifdef HAVE_NU1
        factory_map["nu1"] = ptr_type(new typename TR::template F<NU1>::type());
        maquis::cout << "nu1 ";
#endif
#ifdef HAVE_TrivialGroup
        factory_map["none"] = ptr_type(new typename TR::template F<TrivialGroup>::type());
        maquis::cout << "none ";
#endif
#ifdef HAVE_U1
        factory_map["u1"] = ptr_type(new typename TR::template F<U1>::type());
        maquis::cout << "u1 ";
#endif
#ifdef HAVE_TwoU1
        factory_map["2u1"] = ptr_type(new typename TR::template F<TwoU1>::type());
        maquis::cout << "2u1 ";
#endif
#ifdef HAVE_TwoU1PG
        factory_map["2u1pg"] = ptr_type(new typename TR::template F<TwoU1PG>::type());
        maquis::cout << "2u1pg ";
#endif
#ifdef HAVE_Ztwo
        factory_map["Z2"] = ptr_type(new typename TR::template F<Ztwo>::type());
        maquis::cout << "Z2 ";
#endif
        maquis::cout << std::endl;
        
        
        std::string symm_name;
        if (!parms.is_set("symmetry")) {
#ifdef HAVE_NU1
            symm_name = "nu1";
#else
            if (parms["model_library"] == "alps")
                symm_name = guess_alps_symmetry(parms);
#endif
        } else {
            symm_name = parms["symmetry"].str();
        }
        
        if (factory_map.find(symm_name) != factory_map.end())
            return factory_map[symm_name];
        else
            throw std::runtime_error("Don't know this symmetry group. Please, check your compilation flags.");
#ifdef AMBIENT
        ambient::sync();
#endif
        return ptr_type();
    }

}
