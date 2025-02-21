/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Michele Dolfi <dolfim@phys.ethz.ch>
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

// MPS APPLY OPERATOR
// ------------------
// This example program reads the path to an MPS and applies a site operator
// to all sites in the lattice. Result MPSs are stored back to disk as
// individual files.
// The model, lattice and the name of the operator to apply are specified in
// the code below.


#include <alps/utility/copyright.hpp>
#include <iostream>

#include "dmrg/version.h"

#include <complex>
#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double> matrix;
typedef alps::numeric::matrix<std::complex<double> > cmatrix;

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/models/model.h"
#include "dmrg/models/lattice.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>


/// build with NU1 symmetry
typedef NU1 grp;

template <class Matrix, class SymmGroup>
void run(std::string const& chkp, std::string const& bname)
{
    typedef block_matrix<Matrix, SymmGroup> operator_type;
    
    // Load MPS
    MPS<Matrix, grp> mps;
    load(chkp, mps);
    
    // Define the minimal set of parameters needed to retrieve the operators
    DmrgParameters parms;
    parms.set("LATTICE", "open chain lattice");
    parms.set("L",       10                  );
    parms.set("MODEL",   "spin"              );
    parms.set("CONSERVED_QUANTUMNUMBERS", "Sz");
    parms.set("Sz_total", 0);
    
    // Build lattice and model
    Lattice lattice(parms);
    Model<Matrix, SymmGroup> model(lattice, parms);
    
    const std::string op_name = "Sminus"; // Name of the operator to apply
    
    // Apply the operator at every site and save one MPS it
    for (Lattice::pos_t p=0; p<lattice.size(); ++p) {
        std::cout << boost::format("Applying %s at site %d.") % op_name % p << std::flush;
        operator_type op = model.get_operator(op_name, lattice.get_prop<int>("type", p));
        
        MPS<Matrix, grp> mps_out = mps;
        mps_out.apply(op, p);
        
        const std::string oname = boost::str(boost::format("%s.%s_%d.mps") % bname % op_name % p);
        boost::filesystem::create_directory(oname);
        save(oname, mps_out);
        
        std::cout << " Saved in " << oname << std::endl;
    }
}


int main(int argc, char ** argv)
{
    try {
        std::cout << "ALPS/MPS version " DMRG_VERSION_STRING " (2013-2015)\n"
        << "  Density Matrix Renormalization Group algorithm\n"
        << "  available from http://alps.comp-phys.org/\n"
        << "  copyright (c) 2015 Institute for Theoretical Physics, ETH Zurich\n"
        << "  copyright (c) 2010-2011 by Bela Bauer\n"
        << "  copyright (c) 2011-2015 by Michele Dolfi\n"
        << "  for details see the publication: \n"
        << "  M. Dolfi et al., Computer Physics Communications 185, 3430 (2014).\n"
        << "                   doi: 10.1016/j.cpc.2014.08.019\n"
        << std::endl;
        alps::print_copyright(std::cout);
        
        /// parse options
        std::string chkp, oname;
        
        namespace po = boost::program_options;
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version", "print program version")
        ("complex", "use complex numbers")
        ("mps", po::value<std::string>(&chkp)->required(), "path to chkp of on which operators are applied")
        ("oname", po::value<std::string>(&oname)->required(), "basename for the output");
        po::positional_options_description p;
        p.add("mps", 1);
        p.add("oname", 1);
        
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }
        if (vm.count("version")) {
            std::cout << alps::version_string() << std::endl;
            std::cout << DMRG_VERSION_STRING << std::endl;
            return 0;
        }
        
        po::notify(vm);
        
        /// compute
        if (vm.count("complex")) run<cmatrix, grp>(chkp, oname);
        else                     run<matrix,  grp>(chkp, oname);
        
    } catch (std::exception & e) {
        std::cerr << "Exception thrown:" << std::endl;
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}
