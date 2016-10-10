/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Institute for Theoretical Physics, ETH Zurich
 *               2013-2016 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
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
