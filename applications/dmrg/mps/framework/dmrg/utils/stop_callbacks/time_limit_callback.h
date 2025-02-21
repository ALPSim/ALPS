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

#ifndef MAQUIS_DMRG__UTILS_STOP_CALLBACKS_TIME_LIMIT_CALLBACK_H
#define MAQUIS_DMRG__UTILS_STOP_CALLBACKS_TIME_LIMIT_CALLBACK_H

#include "dmrg/utils/stop_callbacks/stop_callback_base.h"
#include "dmrg/utils/time_stopper.h"
#include "dmrg/utils/time_limit_exception.h"

namespace dmrg {
    class time_limit_callback : public stop_callback_base {
    public:
        time_limit_callback(time_stopper const& stopper)
        : stopper_(stopper)
        { }
        virtual bool operator()(int site, double en)
        {
            return stopper_();
        }
        virtual void throw_exception(int sw, int st) const
        {
            throw time_limit(sw, st);
        }
        virtual stop_callback_base* clone() const
        {
            return new time_limit_callback(*this);
        }
        virtual ~time_limit_callback() {}
    private:
        time_stopper stopper_;
    };
}

#endif
