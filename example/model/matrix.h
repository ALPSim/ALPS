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
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>

template <class T>
struct symbolic_traits {
 static const bool is_symbolic=false;
};

template<>
struct symbolic_traits<alps::Expression> {
 static const bool is_symbolic=true;
};

template <class T, class M = boost::numeric::ublas::matrix<T> >
class HamiltonianMatrix // : public alps::scheduler::Worker
{
public:
  typedef T value_type;
  typedef M matrix_type;
  //FullDiag(const alps::ProcessList&,const alps::Parameters&,int);
  HamiltonianMatrix (const alps::Parameters&);
  //void dostep();
  //void run();
  void output(std::ostream& o)const  { if (!built_) build(); o << matrix_;}
private:
  typedef alps::graph_factory<>::graph_type graph_type;
  
  const graph_type& lattice() const { return graph_factory_.graph();}
  void build() const;
  
  alps::graph_factory<graph_type> graph_factory_;
  alps::ModelLibrary models_;
  alps::Parameters parms_;
  mutable bool built_;
  //mutable boost::multi_array<T,2> matrix_;
  mutable matrix_type matrix_;
};

template <class T, class M>
std::ostream& operator<<(std::ostream& os, const HamiltonianMatrix<T,M>& mat)
{
  mat.output(os);
  return os;
}

template <class T, class M>
HamiltonianMatrix<T,M>::HamiltonianMatrix(const alps::Parameters& p)
  : graph_factory_(p),
    models_(p),
    parms_(p),
    built_(false)
{
}

template <class T, class M>
void HamiltonianMatrix<T,M>::build() const
{
  using namespace alps;
  
  // get Hamilton operator
  HamiltonianDescriptor<short> ham(models_.hamiltonian(parms_["MODEL"]));
  alps::Parameters p(parms_);
  if (!symbolic_traits<T>::is_symbolic)
    p.copy_undefined(ham.default_parameters());
  ham.set_parameters(p);
  
  // get all site matrices
  property_map<site_type_t,const graph_type,int>::type site_type(
            get_or_default(site_type_t(),lattice(),0));

  std::map<int,boost::multi_array<T,2> > site_matrix;
  std::map<int,bool> site_visited;
  for (boost::graph_traits<graph_type>::vertex_iterator it=boost::vertices(lattice()).first; 
       it!=boost::vertices(lattice()).second ; ++it)
    if (!site_visited[site_type[*it]]) {
      int type=site_type[*it];
      std::cout << "Creating site matrix for type " << type << "\n";
      site_visited[type]=true;
      site_matrix.insert(std::make_pair(type,ham.site_term(type).template matrix<T>(
                ham.basis().site_basis(type),models_.simple_operators(),p)));
    }
  
  // get all bond matrices
  property_map<bond_type_t,const graph_type,int>::type bond_type(
            get_or_default(bond_type_t(),lattice(),0));

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
  
  // create basis set
  std::cout << "Creating basis set\n";
  BasisStatesDescriptor<short> basis(ham.basis(),lattice());
  BasisStates<short> states(basis);
  
  // build matrix
  
  std::cout << "Creating matrix\n";
  //matrix_.resize(boost::extents[states.size()][states.size()]);
  matrix_.resize(states.size(),states.size());

  // loop basis states
  for (int i=0;i<states.size();++i) {
    // loop over sites
    int s=0;
    for (boost::graph_traits<graph_type>::vertex_iterator it=boost::vertices(lattice()).first; 
      it!=boost::vertices(lattice()).second ; ++it,++s) {
      // get site basis index
      int is=basis[s].index(states[i][s]);
      // loop over target index
      std::vector<StateDescriptor<short> > state=states[i];
      for (int js=0;js<basis[s].size();++js)
        if (site_matrix[site_type[*it]][is][js]) {
        // build target state
          state[s]=basis[s][js];
	// lookup target state
	  int j = states.index(state);
	// set matrix element
	  //matrix_[i][j]+=site_matrix[site_type[*it]][is][js];
	  matrix_(i,j)+=site_matrix[site_type[*it]][is][js];
	}
    }
    
    // loop over bonds
    for (boost::graph_traits<graph_type>::edge_iterator it=boost::edges(lattice()).first; 
      it!=boost::edges(lattice()).second ; ++it,++s) {
      // get site basis index
      int s1=boost::source(*it,lattice());
      int s2=boost::target(*it,lattice());
      int is1=basis[s1].index(states[i][s1]);
      int is2=basis[s2].index(states[i][s2]);
      // loop over target indices
      std::vector<StateDescriptor<short> > state=states[i];
      for (int js1=0;js1<basis[s1].size();++js1)
        for (int js2=0;js2<basis[s2].size();++js2)
        if (bond_matrix[bond_type[*it]][is1][is2][js1][js2]) {
        // build target state
          state[s1]=basis[s1][js1];
          state[s2]=basis[s2][js2];
	// lookup target state
	  int j = states.index(state);
	// set matrix element
	  //matrix_[i][j]+=bond_matrix[bond_type[*it]][is1][is2][js1][js2];
	  matrix_(i,j)+=bond_matrix[bond_type[*it]][is1][is2][js1][js2];
	}
    }
  }  
  built_ = true;
}
