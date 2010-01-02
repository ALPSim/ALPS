/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2006 by Fabien Alet <alet@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
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

/* $Id$ */

#include "SSE.hpp"
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

// determine sse bond types and calculate their bond matrices
void SSE::initialize_hamiltonian()
{
  typedef boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int,unsigned int> bond_tuple_type;
  number_of_bond_types=0;
  std::map<bond_tuple_type,unsigned int> sse_bond_type;
  std::set<std::string> allops;

  for (bond_iterator it=bonds().first; it!=bonds().second;++it) {
    bond_tuple_type this_bond(inhomogeneous_bond_type(*it),inhomogeneous_site_type(source(*it)),inhomogeneous_site_type(target(*it)),num_neighbors(source(*it)),num_neighbors(target(*it)));

   // TODO : Save here the original bond type if needed
  // original_bond_type[*it] = bond_type[*it];

    std::map<bond_tuple_type,unsigned int>::const_iterator found = sse_bond_type.find(this_bond);
    if (found!=sse_bond_type.end()) {
      bond_type[*it]=found->second;
      number_of_bonds_for_bond_type[found->second]++;
    }
    else {
      // calculate matrix
      // Here uses the "original" bond_type of the lattice (still available)
      matrix_element.insert(std::make_pair(number_of_bond_types,bond_hamiltonian(*it)));      

      site_type_for_bond_type_.resize(number_of_bond_types+1);
      site_type_for_bond_type_[number_of_bond_types].first=site_type(source(*it));
      site_type_for_bond_type_[number_of_bond_types].second=site_type(target(*it));
      
      // Here uses the "original" bond_type of the lattice (still available)
      std::set<std::string> ops = model().bond_term(bond_type(*it)).operator_names(parms);
      allops.insert(ops.begin(),ops.end());
      

      // Now bond_type is given the sse_bond_type
      // save a map to the original bond type
      original_bond_type.push_back(bond_type[*it]);
      // "original" bond type is erased and will not be available
      bond_type[*it]=number_of_bond_types; // new number
      number_of_bonds_for_bond_type.push_back(1); // first bond
      sse_bond_type[this_bond]=number_of_bond_types;
      number_of_bond_types++;

    }
  }

      
  // now look through all terms for raising/lowering operators
  // TODO : For biquadratic terms, one should be able to have two different raising/lowering operators
  
  std::set<int> sitetypes;
  for (site_iterator it=sites().first; it!=sites().second;++it)
    sitetypes.insert(site_type(*it));
  int numsitetypes = *std::max_element(sitetypes.begin(),sitetypes.end());
  matrix_element_lower_.resize(numsitetypes+1);
  matrix_element_raise_.resize(numsitetypes+1);
  for (std::set<int>::iterator it=sitetypes.begin();it!=sitetypes.end();++it) {
    // get matrix 
    bool found_raise=false;
    bool found_lower=false;
    for (std::set<std::string>::iterator t=allops.begin();t!=allops.end();++t) {
      int num_states=model().basis().site_basis(*it).num_states();
      boost::multi_array<double,2> opmatrix = alps::get_matrix(double(),SiteOperator(*t),
        model().basis().site_basis(*it),parms);
      bool is_upper=false;
      bool is_lower=false;
      bool is_other=false;
      for (int i=0;i!=num_states;++i) 
        for (int j=0;j!=num_states;++j) 
          if (opmatrix[i][j]!=0.) {
//            std::cerr << "MAT: " << i << " " << j << " " << opmatrix[i][j] << "\n";
            if (i==j+1) is_lower=true;
            else if (i==j-1) is_upper=true;
            else is_other=true;
          }
      if (is_lower && !is_upper && !is_other) {
        if (found_lower)
          boost::throw_exception(std::runtime_error("Found two lowering operators"));
        found_lower=true;
        std::vector<double> vec (num_states);
        for (int i=0;i<num_states-1;++i)
          vec[i+1]=opmatrix[i+1][i];
        matrix_element_lower_[*it]=vec;
      }
      else if (is_upper && !is_lower && !is_other) {
        if (found_raise)
          boost::throw_exception(std::runtime_error("Found two raising operators"));
        found_raise=true;
        std::vector<double> vec (num_states);
        for (int i=0;i<num_states-1;++i)
          vec[i]=opmatrix[i][i+1];
        matrix_element_raise_[*it]=vec;
      }
    }
    if (!( found_raise && found_lower))
      boost::throw_exception(std::runtime_error("Did not find raising and lowering operators"));
  }
}

// build the bond Hamiltonian matrix
boost::multi_array<double,4> SSE::bond_hamiltonian(const bond_descriptor& b)
{
  unsigned int bond_t        = bond_type(b);
  unsigned int site1_t       = site_type(source(b));
  unsigned int site2_t       = site_type(target(b));
  unsigned int num_neighbor1 = num_neighbors(source(b));
  unsigned int num_neighbor2 = num_neighbors(target(b));
  unsigned int dim1          = model().basis().site_basis(site1_t).num_states();
  unsigned int dim2                         = model().basis().site_basis(site2_t).num_states();
  
  // get site and bond terms
  alps::Parameters p(parms);
  if(inhomogeneous())
    throw_if_xyz_defined(parms,b); // check whether x, y, or z is set

  if (inhomogeneous_sites())
    p << coordinate_as_parameter(source(b)); // set x, y and z
  boost::multi_array<double,2> sitematrix1 = 
    alps::get_matrix(double(),model().site_term(site1_t),model().basis().site_basis(site1_t),p);

  if (inhomogeneous_sites())
    p << coordinate_as_parameter(target(b)); // set x, y and z
  boost::multi_array<double,2> sitematrix2 = 
    alps::get_matrix(double(),model().site_term(site2_t),model().basis().site_basis(site2_t),p);

  if (inhomogeneous_bonds())
    p << coordinate_as_parameter(b); // set x, y and z
  boost::multi_array<double,4> bondtensor = 
    alps::get_matrix(double(),model().bond_term(bond_t),model().basis().site_basis(site1_t),model().basis().site_basis(site2_t),p);
    
  // check bond terms
  for (int i=0;i<dim1;++i)
    for (int j=0;j<dim2;++j)
      for (int k=0;k<dim1;++k)
        for (int l=0;l<dim2;++l)
          if (bondtensor[i][j][k][l]!=0. && (i+j!=k+l ))
            boost::throw_exception(std::runtime_error("Cannot simulate this bond term with the dirloop_sse code"));
  
  
  // add site terms to bond terms
  for (int j=0;j<dim1;++j)
    for (int k=0;k<dim1;++k) {
      // check that site Hamiltonian is diagonal
      if (j!=k && sitematrix1[j][k]!=0.)
        boost::throw_exception(std::runtime_error("Offdiagonal site terms are not implemented in this SSE code"));
      else
        for (int i=0;i<dim2;++i)
          bondtensor[j][i][k][i]+=sitematrix1[j][k]/num_neighbor1;
    }

  for (int j=0;j<dim2;++j)
    for (int k=0;k<dim2;++k) {
      // check that site Hamiltonian is diagonal
      if (j!=k && sitematrix2[j][k]!=0.)
        boost::throw_exception(std::runtime_error("Offdiagonal site terms are not implemented in this SSE code"));
      else
        for (int i=0;i<dim1;++i)
          bondtensor[i][j][i][k]+=sitematrix2[j][k]/num_neighbor2;
    }
  
  return bondtensor;
}

// Calculate the bond offset
void SSE::determine_bonds_offset()
{ 
  energy_offset.resize(number_of_bond_types);
  matrix_sign.resize(boost::extents[number_of_bond_types][maximum_number_of_states][maximum_number_of_states][maximum_number_of_states][maximum_number_of_states]);

  for (alps::uint32_t i=0;i<number_of_bond_types;++i) {
    /********* Old part. Please keep it like this at the moment **************

    double max=-std::numeric_limits<double>::max();
    bool diag=0;
    for (state_type l1=0;l1<number_states_for_site_type_[site_type_for_bond_type_[i].first];l1++) 
      for (state_type l2=0;l2<number_states_for_site_type_[site_type_for_bond_type_[i].second];l2++) { 
        double me=std::fabs(matrix_element[i][l1][l2][l1][l2]);
        if (me!=0) { diag=1;}
        if (me>max)
          max=me; 
      } 
    // Offset is the maximum diagonal matrix element + epsilon
       if ((!(diag)) && (epsilon==0))
         {  boost::throw_exception(std::runtime_error("Hamiltonian looks purely off-diagonal. Please put the parameter EPSILON non zero for sse to work."));
           //cout << "Hello\n"; 
         }
       else {
         energy_offset[i]=max+epsilon;
       }
   
   ***************************************************/

    /******************* Replacement *******************/

    double min=std::numeric_limits<double>::max BOOST_PREVENT_MACRO_SUBSTITUTION ();
    bool diag=0;
    for (state_type l1=0;l1<number_states_for_site_type_[site_type_for_bond_type_[i].first];l1++) 
      for (state_type l2=0;l2<number_states_for_site_type_[site_type_for_bond_type_[i].second];l2++) { 
        double me=-(matrix_element[i][l1][l2][l1][l2]);
        if (me!=0) 
          diag=1;
        if (me<min)
          min=me; 
      }

    if (min>0) { min=0;}

    bool found_non_zero=0;

    if ((!(diag)) && (epsilon==0)) {

       for (state_type l1=0;l1<number_states_for_site_type_[site_type_for_bond_type_[i].first];l1++)
     for (state_type l2=0;l2<number_states_for_site_type_[site_type_for_bond_type_[i].second];l2++) 
       for (state_type l3=0;l3<number_states_for_site_type_[site_type_for_bond_type_[i].first];l3++)
         for (state_type l4=0;l4<number_states_for_site_type_[site_type_for_bond_type_[i].second];l4++) {
           if ((matrix_element[i][l1][l2][l3][l4])!=0) { found_non_zero=1; break;}
         }

    }
    else {found_non_zero=1;}

    if (!(found_non_zero))
      { cout << "WARNING : One bond Hamiltonian seems to be zero. Please modify your model for sse to be optimal.\n";}
    else {
    if ((!(diag)) && (epsilon==0)) 
      boost::throw_exception(std::runtime_error("Hamiltonian looks purely off-diagonal. Please put the parameter EPSILON non zero for sse to work."));

    }
  

    // Offset is the maximum diagonal matrix element + epsilon
    if (is_signed_)
            { cout << "WARNING : Hamiltonian has a sign problem...\n";
              if (epsilon!=0)
                cout << " Please make sure that the value of EPSILON that you use makes the algorithm ergodic\n";
          else
            // random value of epsilon, just to make sure we scan all sign sectors ...
            epsilon=0.0314159;
            }
  


    energy_offset[i]=std::fabs(min)+epsilon;
  
  
   /**************************************************************/

   // invert matrix sign and add diagonal offset
   for (state_type l1=0;l1<number_states_for_site_type_[site_type_for_bond_type_[i].first];l1++)
     for (state_type l2=0;l2<number_states_for_site_type_[site_type_for_bond_type_[i].second];l2++) 
       for (state_type l3=0;l3<number_states_for_site_type_[site_type_for_bond_type_[i].first];l3++)
         for (state_type l4=0;l4<number_states_for_site_type_[site_type_for_bond_type_[i].second];l4++) {
           /************ For sign problem support the full value must be kept or the sign stores separately ****/
           matrix_element[i][l1][l2][l3][l4] = -matrix_element[i][l1][l2][l3][l4];
           if ((l1==l3) && (l2==l4))
             matrix_element[i][l1][l2][l1][l2]+=energy_offset[i];

           /********* For sign problem *****************/
           // Just to be sure there is no sign problem
           matrix_sign[i][l1][l2][l3][l4] = (matrix_element[i][l1][l2][l3][l4]>0 ? 0 : 1);
           matrix_element[i][l1][l2][l3][l4]=std::fabs(matrix_element[i][l1][l2][l3][l4]);
         }
  }
}

// Initialise part of the diagonal update probabilities
void SSE::initialize_diagonal_update_probabilities()
{ 
  // resize the diagonal probability array
  proba_diagonal.resize(boost::extents[number_of_bond_types][maximum_number_of_states][maximum_number_of_states]);

  for (alps::uint32_t i=0;i<number_of_bond_types;++i)
    for (state_type l1=0;l1<number_states_for_site_type_[site_type_for_bond_type_[i].first];l1++)
      for (state_type l2=0;l2<number_states_for_site_type_[site_type_for_bond_type_[i].second];l2++)
        proba_diagonal[i][l1][l2]=std::fabs(matrix_element[i][l1][l2][l1][l2])*num_bonds()*beta;
}

// Initalisation functions
void SSE::initialize_simulation()
{
  initialize_site_states();
  initialize_hamiltonian();

  determine_bonds_offset();
  initialize_diagonal_update_probabilities();
  calculate_proba_worm(); 
  
  print_arrays();
  total_shift=0.;
  for (alps::uint32_t r=0;r<number_of_bond_types;++r)
    total_shift+=energy_offset[r]*number_of_bonds_for_bond_type[r];
}


// Print function
void SSE::print_arrays()
{
if (parms.defined("PRINT_MATRIX")) { 
  cout << "----------- Matrix elements ----------\n";

  for (alps::uint32_t i1=0;i1<number_of_bond_types;++i1) { 
    cout << "**** Bond type " << i1 << endl;
    for (alps::uint32_t i2=0;i2<number_states_for_site_type_[site_type_for_bond_type_[i1].first];++i2) 
      for (alps::uint32_t i3=0;i3<number_states_for_site_type_[site_type_for_bond_type_[i1].second];++i3)
        for (alps::uint32_t i4=0;i4<number_states_for_site_type_[site_type_for_bond_type_[i1].first];++i4)
          for (alps::uint32_t i5=0;i5<number_states_for_site_type_[site_type_for_bond_type_[i1].second];++i5) { 
            cout << "-----" << i2 << " " << i3 << " " << i4 << " " << i5 << endl;
            if (matrix_element[i1][i2][i3][i4][i5]!=0.) {
              cout << i2 << " " << i3 << " " << i4 << " " << i5 << " : " << matrix_element[i1][i2][i3][i4][i5]; 
              if ((i2==i4) && (i3==i5)) 
                cout << " --> " << proba_diagonal[i1][i2][i3];
              cout << endl;
            }
          }
      }
    }

  if (parms.defined("PRINT_PROBA"))  {
    cout << "----------- Worm Proba ----------\n";

  for (alps::uint32_t i1=0;i1<number_of_bond_types;++i1) { 
    cout << "**** Bond type " << i1 << endl;
    for (alps::uint32_t i2=0;i2<number_states_for_site_type_[site_type_for_bond_type_[i1].first];++i2) 
      for (alps::uint32_t i3=0;i3<number_states_for_site_type_[site_type_for_bond_type_[i1].second];++i3)
        for (alps::uint32_t i4=0;i4<number_states_for_site_type_[site_type_for_bond_type_[i1].first];++i4)
          for (alps::uint32_t i5=0;i5<number_states_for_site_type_[site_type_for_bond_type_[i1].second];++i5) { 
            cout << "-----" << i2 << " " << i3 << " " << i4 << " " << i5 << endl;
            
            for (alps::uint32_t r1=0;r1<4;++r1) { 
              for (alps::uint32_t r2=0;r2<4;++r2)
                cout << proba_worm[i1][i2][i3][i4][i5][r1][r2] << " ";
                cout << endl;
            }
          }
      }
    }
}
