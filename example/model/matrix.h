/***************************************************************************
* ALPS application
*
* example/model/matrix.h an matrix class application for full diagonalization
*
* $Id$
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@comp-phys.org>,
**
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
**************************************************************************/

//#include <alps/parameters.h>
#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/osiris/os.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>

template <class T>
struct symbolic_traits {
 BOOST_STATIC_CONSTANT(bool, is_symbolic=false);
};

template <class T>
const bool symbolic_traits<T>::is_symbolic;
template<>
struct symbolic_traits<alps::Expression> {
 BOOST_STATIC_CONSTANT(bool, is_symbolic=true);
};
const bool symbolic_traits<alps::Expression>::is_symbolic;

template <class T, class M = boost::numeric::ublas::matrix<T> >
class HamiltonianMatrix // : public alps::scheduler::Worker
{
public:
  typedef T value_type;
  typedef M matrix_type;

  HamiltonianMatrix (const alps::Parameters&, std::size_t s=1);
  void output(std::ostream& o, bool is_single=false) const 
  { 
    if (!built_) build(); 
    o << matrix_;
  }
  void build() const;

private:
  typedef alps::graph_factory<>::graph_type graph_type;
  
  const graph_type& lattice() const { return graph_factory_.graph();}
  
  alps::graph_factory<graph_type> graph_factory_;
  alps::ModelLibrary models_;
  alps::Parameters parms_;
  mutable bool built_;
  mutable matrix_type matrix_;
};

template <class T, class M>
std::ostream& operator<<(std::ostream& os, const HamiltonianMatrix<T,M>& mat)
{
  mat.output(os);
  return os;
}

template <class T, class M>
HamiltonianMatrix<T,M>::HamiltonianMatrix(const alps::Parameters& p, std::size_t s)
  : graph_factory_(p),
    models_(p),
    parms_(p),
    built_(false),
    matrix_(s,s)
{
}

template <class T, class M>
void HamiltonianMatrix<T,M>::build() const
{
  // using namespace alps;
  
  // get Hamilton operator
  alps::HamiltonianDescriptor<short> ham(models_.hamiltonian(parms_["MODEL"]));
  alps::Parameters p(parms_);
  if (!symbolic_traits<T>::is_symbolic)
    p.copy_undefined(ham.default_parameters());
  ham.set_parameters(p);
  
  // get all site matrices
  alps::property_map<alps::site_type_t, graph_type, int>::const_type
  site_type(alps::get_or_default(alps::site_type_t(), lattice(), 0));

  std::map<int,boost::multi_array<T,2> > site_matrix;
  std::map<int,bool> site_visited;
  double t=alps::dclock();
  for (boost::graph_traits<graph_type>::vertex_iterator it=boost::vertices(lattice()).first; 
       it!=boost::vertices(lattice()).second ; ++it)
    if (!site_visited[site_type[*it]]) {
      int type=site_type[*it];
      std::cout << "Creating site matrix for type " << type << "\n";
      site_visited[type]=true;
      site_matrix.insert(std::make_pair(type,ham.site_term(type).template matrix<T>(
                ham.basis().site_basis(type),models_.simple_operators(),p)));
    }

  std::cout << "Took " << alps::dclock()-t << " seconds\n";
  t=alps::dclock();

  // get all bond matrices
  alps::property_map<alps::bond_type_t,  graph_type, int>::const_type
  bond_type(alps::get_or_default(alps::bond_type_t(),lattice(),0));

  std::map<int,boost::multi_array<T,4> > bond_matrix;
  std::map<boost::tuple<int,int,int>,bool> bond_visited;
  for (boost::graph_traits<graph_type>::edge_iterator it=boost::edges(lattice()).first; 
       it!=boost::edges(lattice()).second ; ++it) {
    int btype  = bond_type[*it];
    int stype1 = site_type[boost::source(*it,lattice())];
    int stype2 = site_type[boost::target(*it,lattice())];
    if (!bond_visited[boost::make_tuple(btype,stype1,stype2)]) {
      std::cout << "Creating bond matrix for type " << btype << "\n";
      bond_visited[boost::make_tuple(btype,stype1,stype2)]=true;
      bond_matrix.insert(std::make_pair(btype,ham.bond_term(btype).template matrix<T>(
                              ham.basis().site_basis(stype1),ham.basis().site_basis(stype2),
                              models_.simple_operators(),p)));
    }
  }

  std::cout << "Took " << alps::dclock()-t << " seconds\n";
  t=alps::dclock();
  
  // create basis set
  std::cout << "Creating basis set\n";
  alps::BasisStatesDescriptor<short> basis(ham.basis(),lattice());
  typedef alps::LookupBasisStates<unsigned int> basis_states_type;
  typedef basis_states_type::value_type state_type;
  basis_states_type states(basis);

  std::cout << "Took " << alps::dclock()-t << " seconds" << std::endl;
  t=alps::dclock();
  
  // build matrix
  
  std::cerr << "Creating matrix\n";
  matrix_.resize(states.size(),states.size());

  // loop over sites
    for (int i=0;i<states.size();++i) {        // loop basis states
      state_type state=states[i];      // get source state
      state_type newstate=state;       // prepare target state
  int s=0;
  for (boost::graph_traits<graph_type>::vertex_iterator it=boost::vertices(lattice()).first; 
      it!=boost::vertices(lattice()).second ; ++it,++s) {
    boost::multi_array<T,2>& mat = site_matrix[site_type[*it]];
      int is=state[s];                         // get site basis index
      for (int js=0;js<basis[s].size();++js) { // loop over target site states
        T val=mat[is][js];                     // get matrix element
        if (val) {                             // if matrix element is nonzero
          newstate[s]=js;                      // build target state
          int j = states.index(newstate);      // lookup target state
          if (j<states.size())
            matrix_(i,j)+=val;                 // set matrix element
        }
      }
    }
  }

  std::cerr << "Took " << alps::dclock()-t << " seconds\n";
  double tot=0.;
  // loop over bonds
    for (int i=0;i<states.size();++i) {      // loop over source states
      state_type state=states[i];           // get source state
      state_type newstate=state;            // prepare target state
  for (boost::graph_traits<graph_type>::edge_iterator it=boost::edges(lattice()).first; 
      it!=boost::edges(lattice()).second ; ++it) {
    boost::multi_array<T,4>& mat = bond_matrix[bond_type[*it]];
    int s1=boost::source(*it,lattice());
    int s2=boost::target(*it,lattice());
      int is1=state[s1];                            // get source site states
      int is2=state[s2];
      for (int js1=0;js1<basis[s1].size();++js1) {  // loop over target site states
        for (int js2=0;js2<basis[s2].size();++js2) {
          T val=mat[is1][is2][js1][js2];            // get matrix element
          if (val) {                                // if nonzero matrix element
            newstate[s1]=js1;                       // build target state
            newstate[s2]=js2;
            int j = states.index(newstate);         // lookup target state
            if (j<states.size())
              matrix_(i,j)+=val;                    // set matrix element
          }
        }
      }
    }
  }  
  built_ = true;
   
  std::cerr << "Took " << alps::dclock()-t << " seconds\n";
}
