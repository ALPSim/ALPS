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

#include "libpscan/scheduler.hpp"

#include <alps/utility/copyright.hpp>
#include <iostream>

#include "dmrg/version.h"

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
        
        Options opt(argc,argv);
        if (opt.valid) {
            Scheduler pscan(opt);
            pscan.run();
        }
    } catch (std::exception & e) {
        std::cerr << "Exception thrown:" << std::endl;
        std::cerr << e.what() << std::endl;
        exit(1);
    }
}

