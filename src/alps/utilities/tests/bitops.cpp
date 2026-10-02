/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/utility/bitops.hpp>
#include <boost/cstdint.hpp>

int main()
{
    using boost::uint32_t;
    using boost::uint64_t;
    
    uint32_t ui  = 5;
    uint64_t uli = 5;
    ui   <<= sizeof(uint32_t)*8-4;
    uli  <<= sizeof(uint64_t)*8-4;

    bool succ = alps::popcnt(ui)  == 2;
    succ     &= alps::popcnt(uli) == 2;
    return succ ? 0 : -1;
}
