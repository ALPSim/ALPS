/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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

#include <iostream>

#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double> matrix;

#include "dmrg/models/model.h"
#include "dmrg/models/lattice.h"

#include "dmrg/utils/DmrgParameters.h"

#include <iomanip>

/// use NU1 symmetry
typedef NU1 grp;


std::ostream& operator<<(std::ostream& os, std::vector<std::string> const& v)
{
    std::copy(v.begin(), v.end()-1, std::ostream_iterator<std::string>(os, ","));
    os << v.back();
    return os;
}

template <class Matrix, class SymmGroup>
void print_test(Model<Matrix, SymmGroup> const& model, std::vector<std::string> const& op_name)
{
    std::cout << op_name << ":\n" << model.get_operator(op_name) << std::endl;
    std::string prod_name = op_name[0]+"(i)";
    for(int i=1;i<op_name.size(); ++i) {
        prod_name += "*" + op_name[i]+"(i)";
    }
    std::cout << prod_name << ":\n" << model.get_operator(prod_name) << std::endl;
}

template <class Matrix, class SymmGroup>
void run1()
{
    // Define the minimal set of parameters
    DmrgParameters parms;
    parms.set("LATTICE", "open chain lattice");
    parms.set("L",       10                  );
    parms.set("MODEL",   "fermion Hubbard"   );
    parms.set("CONSERVED_QUANTUMNUMBERS", "Nup,Ndown");
    parms.set("Nup_total",                5);
    parms.set("Ndown_total",              5);
    
    /// Build lattice and model
    Lattice lattice(parms);
    Model<Matrix, SymmGroup> model(lattice, parms);
    
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("c_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        prod_names.push_back("c_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("c_up");
        prod_names.push_back("cdag_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        prod_names.push_back("c_down");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        prod_names.push_back("cdag_down");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Sz");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Sz");
        prod_names.push_back("Sz");
        print_test(model, prod_names);
    }
}

template <class Matrix, class SymmGroup>
void run2()
{
    // Define the minimal set of parameters
    DmrgParameters parms;
    parms.set("LATTICE", "open chain lattice");
    parms.set("L",       10                  );
    parms.set("MODEL",   "fermion Hubbard"   );
    
    /// Build lattice and model
    Lattice lattice(parms);
    Model<Matrix, SymmGroup> model(lattice, parms);
    
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("c_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        prod_names.push_back("c_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("c_up");
        prod_names.push_back("cdag_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        prod_names.push_back("c_down");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        prod_names.push_back("cdag_down");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Sz");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Sz");
        prod_names.push_back("Sz");
        print_test(model, prod_names);
    }
}


template <class Matrix, class SymmGroup>
void run3()
{
    // Define the minimal set of parameters
    DmrgParameters parms;
    parms.set("LATTICE", "open chain lattice");
    parms.set("L",       10                  );
    parms.set("MODEL",   "alternative fermion Hubbard"   );
    parms.set("CONSERVED_QUANTUMNUMBERS", "N");
    parms.set("N_total",                5);
    
    /// Build lattice and model
    Lattice lattice(parms);
    Model<Matrix, SymmGroup> model(lattice, parms);
    
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("c_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        prod_names.push_back("c_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("c_up");
        prod_names.push_back("cdag_up");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        prod_names.push_back("c_down");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("cdag_up");
        prod_names.push_back("cdag_down");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Sz");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Sz");
        prod_names.push_back("Sz");
        print_test(model, prod_names);
    }
}


template <class Matrix, class SymmGroup>
void run4()
{
    // Define the minimal set of parameters
    DmrgParameters parms;
    parms.set("LATTICE", "open chain lattice");
    parms.set("L",       10                  );
    parms.set("MODEL",   "spin"   );
    parms.set("CONSERVED_QUANTUMNUMBERS", "Z");
    parms.set("Sz_total",                0);
    
    /// Build lattice and model
    Lattice lattice(parms);
    Model<Matrix, SymmGroup> model(lattice, parms);
    
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Splus");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Sminus");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Splus");
        prod_names.push_back("Sminus");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Sminus");
        prod_names.push_back("Splus");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Sz");
        print_test(model, prod_names);
    }
    {
        std::vector<std::string> prod_names;
        prod_names.push_back("Sz");
        prod_names.push_back("Sz");
        print_test(model, prod_names);
    }
}


int main(int argc, char ** argv)
{
    std::cout << "----------------------------------------" << std::endl;
    run1<matrix, grp>();
    std::cout << "----------------------------------------" << std::endl;
    run2<matrix, grp>();
    std::cout << "----------------------------------------" << std::endl;
    run3<matrix, grp>();
    std::cout << "----------------------------------------" << std::endl;
    run4<matrix, grp>();
}
