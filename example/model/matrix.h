/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2005 by Matthias Troyer <troyer@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

#include <alps/expression.h>
#include <alps/model/model_helper.h>
#include <alps/lattice/graph_helper.h>

#include <alps/numeric/is_nonzero.hpp>
#include <alps/type_traits/is_symbolic.hpp>
#include <alps/multi_array.hpp>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>


template <class T, class M = boost::numeric::ublas::matrix<T> >
class HamiltonianMatrix : public alps::graph_helper<>
{
public:
  typedef T value_type;
  typedef M matrix_type;

  HamiltonianMatrix (const alps::Parameters&);
  void output(std::ostream& o, bool /* is_single */ = false) const 
  { 
    if (!built_) build(); 
    o << matrix_;
  }

  int rows() { return matrix_.size1(); }
  int cols() { return matrix_.size2(); }
  const T& operator()(int i, int j) const { return matrix_(i,j); }


  void build() const;

private:
  typedef alps::graph_helper<>::graph_type graph_type;
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
HamiltonianMatrix<T,M>::HamiltonianMatrix(const alps::Parameters& p)
  : alps::graph_helper<>(p),
    models_(p),
    parms_(p),
    built_(false)
{
}

template <class T, class M>
void HamiltonianMatrix<T,M>::build() const
{
#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
  using namespace alps;
#endif

  // get Hamilton operator
  alps::HamiltonianDescriptor<short> ham(models_.get_hamiltonian(*this,parms_,alps::is_symbolic<T>::type::value));
    
  // get all site matrices
  std::map<unsigned int,alps::multi_array<T,2> > site_matrix;
  std::map<unsigned int,bool> site_visited;
  
  alps::Disorder::seed(parms_.value_or_default("DISORDER_SEED",0));

  
  alps::Parameters parms(parms_);
  for (site_iterator it=sites().first; it!=sites().second ; ++it)
    if (!site_visited[inhomogeneous_site_type(*it)]) {
      unsigned int inhomogeneous_type=inhomogeneous_site_type(*it);
      unsigned int type=site_type(*it);
      site_visited[inhomogeneous_type]=true;
      // set coordinate in case of site disorder
      if (inhomogeneous_sites()) {
        throw_if_xyz_defined(parms_,*it); // check whether x, y, or z is set
        parms << coordinate_as_parameter(*it); // set x, y and z
      }
      site_matrix.insert(std::make_pair(inhomogeneous_type,get_matrix(T(),ham.site_term(type),
        ham.basis().site_basis(type),parms)));
    }

  // get all bond matrices
  std::map<boost::tuple<unsigned int,unsigned int,unsigned int>,alps::multi_array<std::pair<T,bool>,4> > bond_matrix;
  std::map<boost::tuple<unsigned int,unsigned int,unsigned int>,bool> bond_visited;
  for (bond_iterator it=bonds().first; it!=bonds().second ; ++it) {
    unsigned int inhomogeneous_btype  = inhomogeneous_bond_type(*it);
    //unsigned int inhomogeneous_stype1 = inhomogeneous_site_type(source(*it));
    //unsigned int inhomogeneous_stype2 = inhomogeneous_site_type(target(*it));
    unsigned int btype  = bond_type(*it);
    unsigned int stype1 = site_type(source(*it));
    unsigned int stype2 = site_type(target(*it));
    boost::tuple<unsigned int,unsigned int,unsigned int> type(inhomogeneous_btype,stype1,stype2);
    if (!bond_visited[type]) {
      if (inhomogeneous_bonds()) {
        throw_if_xyz_defined(parms_,*it); // check whether x, y, or z is set
        parms << coordinate_as_parameter(*it); // set x, y and z
      }
      bond_visited[type]=true;
      bond_matrix.insert(std::make_pair(type,get_fermionic_matrix(T(),ham.bond_term(btype),
        ham.basis().site_basis(stype1),ham.basis().site_basis(stype2),parms)));
    }
  }

  // create basis set
  std::cerr << "Creating basis set\n";
  alps::basis_states_descriptor<short> basis(ham.basis(),graph());
  //typedef alps::LookupBasisStates<unsigned int> basis_states_type;
  typedef alps::basis_states<short> basis_states_type;
  typedef basis_states_type::value_type state_type;
  basis_states_type states(basis);

  // build matrix
  
  matrix_.resize(states.size(),states.size(),false);
  matrix_.clear();
  std::cout << states << "\n";

  // loop over sites
    for (int i=0;i<states.size();++i) {        // loop basis states
      state_type state=states[i];              // get source state
  int s=0;
  for (site_iterator it=sites().first; it!=sites().second ; ++it,++s) {
    alps::multi_array<T,2>& mat = site_matrix[inhomogeneous_site_type(*it)];
      int is=state[s];                         // get site basis index
      for (int js=0;js<basis[s].size();++js) { // loop over target site states
        T val=mat[is][js];                     // get matrix element
        if (alps::numeric::is_nonzero(val)) {           // if matrix element is nonzero
          state_type newstate=state;
          newstate[s]=js;                      // build target state
          int j = states.index(newstate);      // lookup target state
          if (j<states.size()) {
            matrix_(i,j)+=val;                 // set matrix element
            alps::simplify(matrix_(i,j));
          }
        }
      }
    }
  }

  // loop over bonds
  for (int i=0;i<states.size();++i) {     // loop over source states
    state_type state=states[i];           // get source state
    for (bond_iterator it=bonds().first; it!=bonds().second ; ++it) {
      int s1=source(*it);
      int s2=target(*it);
      alps::multi_array<std::pair<T,bool>,4>& mat = bond_matrix[boost::make_tuple(inhomogeneous_bond_type(*it),site_type(s1),site_type(s2))];
      int is1=state[s1];                           // get source site states
      int is2=state[s2];
      for (int js1=0;js1<basis[s1].size();++js1) { // loop over target site states
        for (int js2=0;js2<basis[s2].size();++js2) {
          T val=mat[is1][is2][js1][js2].first;     // get matrix element
          if (alps::numeric::is_nonzero(val)) {             // if nonzero matrix element
            state_type newstate=state;             // prepare target state
            newstate[s1]=js1;                      // build target state
            newstate[s2]=js2;
            int j = states.index(newstate);        // lookup target state
            if (j<states.size()) {
              if (mat[is1][is2][js1][js2].second) {
                // calculate fermionic sign
                int start = std::min(s1,s2);
                int end = std::max(s1,s2);
                bool f=(s2>=s1);
                for (int i=start;i<end;++i)
                  if (is_fermionic(ham.basis().site_basis(site_type(i)),basis[i][state[i]]))
                    f=!f;
                if (f)
                  val=-val;
              }
              matrix_(i,j)+=val;                    // set matrix element
              alps::simplify(matrix_(i,j));
            }
          }
        }
      }
    }
  }  
  built_ = true;
}
