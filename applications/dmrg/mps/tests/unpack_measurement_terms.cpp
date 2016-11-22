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

/// use NU1 symmetry
typedef NU1 grp;

template <class Matrix, class SymmGroup>
void print_test(Model<Matrix, SymmGroup> const& model, std::string const& op_name)
{
    std::cout << op_name << ":\n" << model.unpack_measurement_terms(op_name) << std::endl;
}

template <class Matrix, class SymmGroup>
void run1()
{
    // Define the minimal set of parameters
    DmrgParameters parms;
    parms.set("LATTICE", "open chain lattice");
    parms.set("L",       10                  );
    parms.set("MODEL",   "fermion Hubbard"   );
    
    /// Build lattice and model
    Lattice lattice(parms);
    Model<Matrix, SymmGroup> model(lattice, parms);
    
    print_test(model, "n_up");
    print_test(model, "Splus");
    print_test(model, "fermion_hop");
    print_test(model, "double_occupancy");
    print_test(model, "exchange_xy");
    print_test(model, "Sx");
    print_test(model, "biquadratic");
}

template <class Matrix, class SymmGroup>
void run2()
{
    // Define the minimal set of parameters
    DmrgParameters parms;
    parms.set("LATTICE", "open chain lattice");
    parms.set("L",       10                  );
    parms.set("MODEL",   "t-J"   );
    
    /// Build lattice and model
    Lattice lattice(parms);
    Model<Matrix, SymmGroup> model(lattice, parms);
    
    print_test(model, "n_up");
    print_test(model, "Splus");
    print_test(model, "fermion_hop");
    print_test(model, "double_occupancy");
    print_test(model, "exchange_xy");
    print_test(model, "Sx");
}

// TODO: test when one site op is not defined in some site type basis
template <class Matrix, class SymmGroup>
void run3()
{
}


int main(int argc, char ** argv)
{
    run1<matrix, grp>();
    run2<matrix, grp>();
}
