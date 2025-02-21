/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Institute for Theoretical Physics, ETH Zurich
 *               2016-2016 by Michele Dolfi <dolfim@phys.ethz.ch>
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


#include <iostream>
#include <string>

#include "dmrg/evolve/trotter_decomposer.h"

void print_decomp(unsigned nterms, std::string const& te_order)
{
    std::cout << "====================================\n"
              << " nterms=" << nterms << ", order=" << te_order << "\n"
              << "====================================\n"
              << std::endl;
    
    trotter_decomposer decomp(nterms, te_order, true);
    
    std::cout << "## TERMS:\n";
    for (unsigned i=0; i<decomp.size(); ++i)
        std::cout << " " << i << " : " << "x=" << decomp.trotter_term(i).first << ", alpha=" << decomp.trotter_term(i).second << std::endl;
    
    {
        std::cout << "## SIMPLE SEQUENCE:" << std::endl;
        trotter_decomposer::sequence_type Useq = decomp.simple_terms_sequence();
        for (unsigned i=0; i<Useq.size(); ++i)
            std::cout << "  " << Useq[i] << std::endl;
    }
    {
        std::cout << "## INITIAL SEQUENCE:" << std::endl;
        trotter_decomposer::sequence_type Useq = decomp.initial_terms_sequence();
        for (unsigned i=0; i<Useq.size(); ++i)
            std::cout << "  " << Useq[i] << std::endl;
    }
    {
        std::cout << "## DOUBLE SEQUENCE:" << std::endl;
        trotter_decomposer::sequence_type Useq = decomp.double_terms_sequence();
        for (unsigned i=0; i<Useq.size(); ++i)
            std::cout << "  " << Useq[i] << std::endl;
    }
    {
        std::cout << "## FINAL SEQUENCE:" << std::endl;
        trotter_decomposer::sequence_type Useq = decomp.final_terms_sequence();
        for (unsigned i=0; i<Useq.size(); ++i)
            std::cout << "  " << Useq[i] << std::endl;
    }
}

int main(int argc, char ** argv)
{
    print_decomp(2, "first");
    print_decomp(2, "second");
    print_decomp(2, "fourth");

    print_decomp(4, "first");
    print_decomp(4, "second");
    print_decomp(4, "fourth");
}
