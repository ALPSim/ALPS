/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include "dmrg/models/alps/symm_handler.hpp"
#include "dmrg/block_matrix/symmetry/z2.h"
	
// Symmetry dependent implementation

// Z2 Symmetry
template <>
Ztwo::charge init_charge<Ztwo> (const alps::Parameters& parms, std::map<std::string, int> const& all_conserved_qn)
{
    typedef std::map<std::string, int> qn_map_type;
    assert(all_conserved_qn.size() == 1);
    qn_map_type::const_iterator it = all_conserved_qn.begin();
    int tmp = alps::evaluate<double>(static_cast<std::string>(parms[it->first+"_total"]), parms);
    if (!(tmp == 0 || tmp == 1))
        throw std::runtime_error("Invalid value for " + it->first + "_total");
    return (tmp == 1) ? Ztwo::Minus : Ztwo::Plus;
}

template <>
Ztwo::charge state_to_charge<Ztwo>(alps::site_state<short> const & state, alps::SiteBasisDescriptor<short> const& b,
                               std::map<std::string, int> const& all_conserved_qn)
{
    typedef std::map<std::string, int> qn_map_type;
    int tmp = 0;
    for (typename alps::SiteBasisDescriptor<short>::const_iterator it = b.begin(); it != b.end(); ++it) {
        qn_map_type::const_iterator match = all_conserved_qn.find(it->name());
        if (match != all_conserved_qn.end())
            tmp = detail::to_integer( get_quantumnumber(state, it->name(), b) );
    }
    return (tmp % 4 == 0 ? Ztwo::Plus : Ztwo::Minus);
}
