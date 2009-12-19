/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
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

#include "WRun.h"
#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <iomanip>

using namespace alps;

//---------------------------------------------------------------------------------------------------

void WRun::initialize_hamiltonian() {

  // Fills the following data structures
  //
  // number_of_bond_types .. number of bond types
  //
  // bond_tye[s] ... the new bond type of a site (taking into account different
  //                 Hilbert spaces on neighboring sites, 
  //                 values 0 .. number_of_bond_types-1
  //
  // site_type_for_bond_type_[bt] ... pair of site types of the two sites
  //
  // site_matrix[st] .. a vector of the diagonal matrix element, 
  //                    indexed by the local state (mu and U in Lode's code)
  //
  // diagonal_matrix_element[bt] .. diagonal interaction between two sites, 
  //                                indexed by the two local states
  //                                (Vnn in Lode's code)
  //
  // hopping_matrix[bt] ... hopping matrix element 
  //
  // matrix_element_raise_[st] ... matrix element of raising operator, 
  // matrix_element_lower_[st] ... matrix element of lowering operator, 
  //                               indexed by local state
  
  typedef boost::tuple<int,int,int> bond_tuple_type;
  std::map<bond_tuple_type,int> worm_bond_type;
  std::set<std::string> allops;
  std::map<int,boost::multi_array<double,4> > matrix_element;
  
  number_of_bond_types = 0;

  // iterate all bonds
  for(bond_iterator it=bonds().first; it!=bonds().second; ++it) {
    bond_tuple_type this_bond(inhomogeneous_bond_type(*it),site_type(source(*it)),site_type(target(*it)));

    std::map<bond_tuple_type,int>::const_iterator found = worm_bond_type.find(this_bond);          
    if (found!=worm_bond_type.end())
      bond_type[*it] = found->second;
    else {
      // calculate matrix
      // Here uses the "original" bond_type of the lattice (still available)
      matrix_element.insert(std::make_pair(number_of_bond_types,bond_hamiltonian(*it)));      

      site_type_for_bond_type_.resize(number_of_bond_types+1);
      site_type_for_bond_type_[number_of_bond_types].first  = site_type(source(*it));
      site_type_for_bond_type_[number_of_bond_types].second = site_type(target(*it));
      
      // Here uses the "original" bond_type of the lattice (still available)
      std::set<std::string> ops = model().bond_term(bond_type(*it)).operator_names(parms);
      allops.insert(ops.begin(),ops.end());

      // lattice bond_type is relabeled as worm_bond_type, "original" lattice bond type is erased
      bond_type[*it]=number_of_bond_types; // new label (worm_bond_type)
      worm_bond_type[this_bond] = number_of_bond_types;
      number_of_bond_types++;
    }
  }

  // iterate all sites and build site matrix
  std::set<unsigned int> site_types;
  for(site_iterator it=sites().first; it!=sites().second; ++it) {
    int this_site_type = inhomogeneous_site_type(*it);
    if (site_types.find(this_site_type)==site_types.end()) {
      if(this_site_type >= site_matrix.size())
        site_matrix.resize(this_site_type+1);
      site_matrix[this_site_type]=site_hamiltonian(*it);
    }
  }
  
  // now look through all terms for raising/lowering operators
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

  // find hopping matrix elements
  for(bond_iterator it=bonds().first; it!=bonds().second; ++it) {
    int bt = bond_type[*it];
    int source_site_type = site_type(source(*it));
    int target_site_type = site_type(target(*it));
    
    double t = matrix_element[bt][0][1][1][0] / ( matrix_element_raise_[source_site_type][0] * matrix_element_lower_[target_site_type][1] );
    
    // consistency check
    int num_source_states = model().basis().site_basis(source_site_type).num_states();
    int num_target_states = model().basis().site_basis(target_site_type).num_states();
    const double epsilon = 1e-8;
    double tp;

    diagonal_matrix_element.insert(std::make_pair(bt, boost::multi_array<double,2>(boost::extents[num_source_states][num_target_states])));
        
    for(int i=0; i!=num_source_states; ++i) {
      for(int j=0; j!=num_target_states; ++j) {
        // set diagonal matrix element
    diagonal_matrix_element[bt][i][j] = matrix_element[bt][i][j][i][j];

    // check hopping 
    if( (i+1<num_source_states) && j ) {
      tp = matrix_element[bt][i][j][i+1][j-1] / ( matrix_element_raise_[source_site_type][i] * matrix_element_lower_[target_site_type][j] );
      
      if(fabs(t-tp)>epsilon) 
        boost::throw_exception(std::runtime_error("Found inconsistent hopping matrix elements."));
    }
      }
    }

    // if everything is fine, add hopping matrix element
    if(bt >= hopping_matrix.size()) {
      hopping_matrix.resize(bt+1);
    }
    hopping_matrix[bt] = t;
  }
}   // WRun::initialize_hamiltonian

//---------------------------------------------------------------------------------------------------

boost::multi_array<double,4> WRun::bond_hamiltonian(const bond_descriptor& b) {
  // get bond terms
  unsigned int bond_t        = bond_type(b);
  unsigned int site1_t       = site_type(source(b));
  unsigned int site2_t       = site_type(target(b));

  alps::Parameters p(parms);

  if (inhomogeneous_bonds()) {
    throw_if_xyz_defined(parms,b); // check whether x, y, or z is set
    p << coordinate_as_parameter(b); // set x, y and z
  }
  
  boost::multi_array<double,4> bondtensor = 
    alps::get_matrix(double(),model().bond_term(bond_t),model().basis().site_basis(site1_t),model().basis().site_basis(site2_t),p);

  unsigned int dim1 = model().basis().site_basis(site1_t).num_states();
  unsigned int dim2 = model().basis().site_basis(site2_t).num_states();

  for (int i=0;i<dim1;++i)
    for (int j=0;j<dim2;++j)
      for (int k=0;k<dim1;++k)
        for (int l=0;l<dim2;++l)
          if (bondtensor[i][j][k][l]!=0. && (i+j!=k+l ))
            boost::throw_exception(std::runtime_error("Cannot simulate this bond term with the worm code"));

  return bondtensor;
}   // WRun::bond_hamiltonian

std::vector<double> WRun::site_hamiltonian(const site_descriptor& s) {
  unsigned int site_t = site_type(s);
  alps::Parameters p(parms);

  if(inhomogeneous_sites()) {
    throw_if_xyz_defined(parms,s);   // check whether x, y, or z is set
    p << coordinate_as_parameter(s); // set x, y and z
  }

  boost::multi_array<double,2> this_site_matrix = 
    alps::get_matrix(double(), model().site_term(site_t),
             model().basis().site_basis(site_t), p);
             
  std::vector<double> diag;
  for (int i=0;i<this_site_matrix.shape()[0];++i)
    diag.push_back(this_site_matrix[i][i]);
    
  return diag;
}   // WRun::site_hamiltonian

//---------------------------------------------------------------------------------------------------

void WRun::print_hamiltonian() {
  std::cout << "Hamiltonian generated by model library" << std::endl
            << " (*) hopping matrix elements " << std::endl;
  for(bond_iterator it=bonds().first; it!=bonds().second; ++it) {
    std::cout << "     t = " << hopping_matrix_element(*it) 
          << "  for bond " << *it << " bond type = " << bond_type(*it) << std::endl; 
  }
  std::cout << " (*) site matrix elements " << std::endl;
  for(site_iterator it=sites().first; it!=sites().second; ++it) {
    std::cout << "     site " << *it 
          << " site type = " << site_type(*it) 
          << " dis. site type = " << inhomogeneous_site_type(*it) << std::endl; 
    int this_site_type = inhomogeneous_site_type(*it);
    int nstates = site_matrix[this_site_type].size();
    for(int i=0; i<nstates; ++i) {
      std::cout << "       state = " << i << "    onsite energy = " << onsite_energy(i, *it) 
        << "    creation = " << creation_matrix_element(i, true, *it) 
        << "    annihilation = " << creation_matrix_element(i, false, *it)<< std::endl;
    }
  }
  /*
  std::cout << " (*) diagonal matrix elements " << std::endl;
  for(int i=0; i<number_states_for_site_type_[0]; ++i) {
    for(int j=0; j<number_states_for_site_type_[0]; ++j) {
      std::cout << "     state = (" << i << "," << j << ")    neighbor energy = " << neighbor_energy(i, j, 0) << std::endl; 
    }
    }
  */
}   // WRun::print_hamiltonian
