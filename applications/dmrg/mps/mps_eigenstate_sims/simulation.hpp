/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef ALPS_MPS_SIM_RUN_H
#define ALPS_MPS_SIM_RUN_H

#include <boost/shared_ptr.hpp>
#include "dmrg/utils/DmrgParameters.h"

enum run_type {optim_and_measure, optim_only, measure_only};

struct simulation_base {
    virtual ~simulation_base() {}
    virtual void run(DmrgParameters & parms, bool write_xml, run_type rt) =0;
};

template <class SymmGroup>
struct simulation : public simulation_base {
    void run(DmrgParameters & parms, bool write_xml, run_type rt);
};

struct simulation_traits {
    typedef boost::shared_ptr<simulation_base> shared_ptr;
    template <class SymmGroup> struct F {
        typedef simulation<SymmGroup> type;
    };
};

#endif
