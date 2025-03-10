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

#include <alps/utility/copyright.hpp>
#include <iostream>

#include "dmrg/version.h"

#include <complex>
#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double> matrix;
typedef alps::numeric::matrix<std::complex<double> > cmatrix;

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"

#include <boost/program_options.hpp>


/// build with NU1 symmetry
typedef NU1 grp;


template <class Matrix, class SymmGroup>
void run (std::string const& chkp1, std::string const& chkp2)
{
    MPS<Matrix, grp> mps1, mps2;
    load(chkp1, mps1);
    load(chkp2, mps2);
    
    std::cout << "<mps1 | mps2> = " << overlap(mps1, mps2) << std::endl;
}

int main(int argc, char ** argv)
{
    try {
        std::cout << "ALPS/MPS version " DMRG_VERSION_STRING " (2013-2014)\n"
                  << "  Density Matrix Renormalization Group algorithm\n"
                  << "  available from http://alps.comp-phys.org/\n"
                  << "  copyright (c) 2013 Institute for Theoretical Physics, ETH Zurich\n"
                  << "  copyright (c) 2010-2011 by Bela Bauer\n"
                  << "  copyright (c) 2011-2013 by Michele Dolfi\n"
                  << "  for details see the publication: \n"
                  << "  M. Dolfi et al., Computer Physics Communications 185, 3430 (2014).\n"
                  << "                   doi: 10.1016/j.cpc.2014.08.019\n"
                  << std::endl;
        alps::print_copyright(std::cout);
        
        /// parse options
        std::string chkp1, chkp2;
        
        namespace po = boost::program_options;
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version", "print program version")
        ("complex", "use complex numbers")
        ("mps1", po::value<std::string>(&chkp1)->required(), "path to chkp of mps1")
        ("mps2", po::value<std::string>(&chkp2)->required(), "path to chkp of mps2");
        po::positional_options_description p;
        p.add("mps1", 1);
        p.add("mps2", 1);
        
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
        if (vm.count("complex")) run<cmatrix, grp>(chkp1, chkp2);
        else                     run<matrix,  grp>(chkp1, chkp2);
        
    } catch (std::exception & e) {
        std::cerr << "Exception thrown:" << std::endl;
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}
