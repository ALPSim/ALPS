/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Institute for Theoretical Physics, ETH Zurich
 *               2016-2016 by Michele Dolfi <dolfim@phys.ethz.ch>
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
