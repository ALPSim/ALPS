/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Institute for Theoretical Physics, ETH Zurich
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

#ifndef MAQUIS_DMRG__UTILS_STOP_CALLBACKS_REL_ENERGY_CALLBACK_H
#define MAQUIS_DMRG__UTILS_STOP_CALLBACKS_REL_ENERGY_CALLBACK_H

#include "dmrg/utils/stop_callbacks/stop_callback_base.h"
#include "dmrg/utils/simulation_terminated_exception.h"

namespace dmrg {
    class rel_energy_callback : public stop_callback_base {
    public:
        rel_energy_callback(int N, double rel_en_thresh_, std::string const& en_thresh_at)
        : valid(false)
        , rel_en_thresh(rel_en_thresh_)
        , at_site(-1)
        {
            if (en_thresh_at == "half") {
                at_site = N / 2 - 1;
            } else if (en_thresh_at == "end") {
                at_site = N - 1;
            }
        }
        virtual bool operator()(int site, double en_new)
        {
            if (at_site < 0) return false; // never check
            
            if (!valid && at_site == site) { // set en_prev the first time
                en_prev = en_new;
                valid = true;
                return false;
            }
            
            if (valid && at_site == site) {
                bool stop = std::abs((en_prev-en_new) / en_new) < rel_en_thresh;
                en_prev = en_new;
                return stop;
            } else {
                return false;
            }
        }
        virtual void throw_exception(int sw, int st) const
        {
            throw simulation_terminated(sw, st, "Rel Energy converged.");
        }
        virtual stop_callback_base* clone() const
        {
            return new rel_energy_callback(*this);
        }
        virtual ~rel_energy_callback() {}
    private:
        bool valid;
        double rel_en_thresh;
        int at_site;
        double en_prev;
    };
}

#endif
