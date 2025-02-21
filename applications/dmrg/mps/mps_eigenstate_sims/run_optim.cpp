/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include <boost/filesystem/fstream.hpp>
#include <string>

#include "libpscan/run_sim.hpp"

#include "mps_eigenstate_sims/simulation.hpp"
#include "dmrg/sim/symmetry_factory.h"


void run_sim(const boost::filesystem::path& infile, const boost::filesystem::path& outfile,
             bool write_xml, double time_limit)
{
    maquis::cout.precision(10);
    
    /// Load parameters
    DmrgParameters parms;
    {
        boost::filesystem::ifstream param_file(infile);
        if (!param_file)
            throw std::runtime_error(std::string("Could not open parameter file ") + infile.string() +".");
        alps::Parameters p; p.extract_from_xml(param_file);
        parms = DmrgParameters(p);
    }
    
    /// Match parameters of ALPS DMRG
    if (!parms.defined("chkpfile"))
        parms.set("chkpfile",   (outfile.parent_path() / outfile.stem()).string() + ".chkp");
    parms.set("resultfile", (outfile.parent_path() / outfile.stem()).string() + ".h5");
    parms.set("run_seconds", time_limit);
    if (parms.defined("SWEEPS"))           parms.set("nsweeps", int(parms["SWEEPS"]));
    if (parms.defined("TRUNCATION_ERROR")) parms.set("truncation_final", parms["TRUNCATION_ERROR"]);
    
    if (parms.defined("NUMSTATES")) {
        parms.set("max_bond_dimension", int(parms["NUMSTATES"]));
    } else if (parms.defined("MAXSTATES")) {
        std::stringstream ss;
        
        int init_bond = 20;
        double step = (double(parms["MAXSTATES"]) - init_bond) / (parms["nsweeps"]-1);
        
        ss << init_bond;
        for (int sweep = 1; sweep < parms["nsweeps"]; ++sweep) {
            ss << "," << int(init_bond + step*sweep);
        }
        
        parms.set("sweep_bond_dimensions", ss.str());
    } else if (parms.defined("STATES")) {
        parms.set("sweep_bond_dimensions", parms["NUMSTATES"].str());
    }
    if (!parms.defined("max_bond_dimension") && !parms.defined("sweep_bond_dimensions"))
        throw std::runtime_error("Number of renormalized states not set.");
    
    /// Start simulation
    simulation_traits::shared_ptr sim = dmrg::symmetry_factory<simulation_traits>(parms);
    sim->run(parms, write_xml, optim_and_measure);
}
