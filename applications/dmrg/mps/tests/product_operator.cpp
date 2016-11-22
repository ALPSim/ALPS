/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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
