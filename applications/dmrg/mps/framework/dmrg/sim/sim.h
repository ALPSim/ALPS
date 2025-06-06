/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_SIM_H
#define APP_SIM_H

#include "dmrg/version.h"

#include <cmath>
#include <iterator>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>

#include "dmrg/utils/DmrgParameters.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/twositetensor.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"

#include "dmrg/utils/random.hpp"
#include "dmrg/utils/time_stopper.h"
#include "utils/timings.h"

#include "dmrg/models/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"


class abstract_sim {
public:
    virtual ~abstract_sim() {}
    virtual void run() =0;
};


template <class Matrix, class SymmGroup>
class sim : public abstract_sim {
public:
    sim(DmrgParameters const &);
    virtual ~sim();
    
    virtual void run() =0;
    
protected:
    typedef boost::ptr_vector<measurement<Matrix, SymmGroup> > measurements_type;
    typedef std::map<std::string, int> status_type;
    
    virtual std::string results_archive_path(status_type const&) const;
    
    measurements_type iteration_measurements(int sweep);
    virtual void measure(std::string archive_path, measurements_type & meas);
    // TODO: can be made const, now only problem are parameters
    
    virtual void checkpoint_simulation(MPS<Matrix, SymmGroup> const& state, status_type const&);
    
protected:
    DmrgParameters parms;
    
    int init_sweep, init_site;
    bool restore;
    bool dns;
    std::string chkpfile;
    std::string rfile;
    
    time_stopper stop_callback;
    
    Lattice lat;
    Model<Matrix, SymmGroup> model;
    MPS<Matrix, SymmGroup> mps;
    MPO<Matrix, SymmGroup> mpo, mpoc;
    measurements_type all_measurements, sweep_measurements;
};

#include "sim.hpp"
#endif
