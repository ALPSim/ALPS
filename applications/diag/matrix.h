/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>
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

#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/scheduler.h>
#include <alps/matrix_as_vector.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/multi_array.hpp>
#include <boost/tokenizer.hpp>
#include <boost/timer.hpp>
#include <cmath>
#include <cstddef>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>


template <class T>
T conj (T x) { return x;}

template <class T>
std::complex<T> conj (std::complex<T> x) { return std::conj(x);}



template <class T>
struct optional_real {
  template <class U>
  static U real(const U& x) { return x;}
  template <class U>
  static U real(const std::complex<U>& x) { return std::real(x);}
};

template <class T>
struct optional_real<std::complex<T> > {
  template <class U>
  static U real(const U& x) { return x;}
};


template <class T, class M>
class HamiltonianMatrix : public alps::scheduler::Task, public alps::graph_helper<>, public alps::model_helper<>
{
public:
  typedef T value_type;
  typedef typename alps::type_traits<value_type>::norm_t magnitude_type;
  typedef M matrix_type;
  typedef alps::model_helper<>::half_integer_type half_integer_type;
  typedef typename alps::graph_helper<>::graph_type graph_type;
  typedef typename alps::graph_helper<>::site_descriptor site_descriptor;
  typedef typename alps::graph_helper<>::bond_descriptor bond_descriptor;
  typedef typename alps::graph_helper<>::vector_type vector_type;
  typedef alps::basis_states<short> basis_states_type;
  typedef alps::bloch_basis_states<short> bloch_basis_states_type;
  typedef basis_states_type::value_type state_type;
  
  HamiltonianMatrix (const alps::ProcessList& , const boost::filesystem::path& );
  matrix_type& matrix() {if (!built_matrix_) build(); return matrix_;}
  const matrix_type& matrix() const {if (!built_matrix_) build(); return matrix_;}
  std::size_t dimension() const { if (!built_basis_) build_basis(); return use_bloch_ ? bloch_states.size() : states.size();}
  void dostep();


  template <class STATES, class V, class W>
  void apply_operator(const STATES&, const alps::SiteOperator& op, site_descriptor s, const V&, W&) const;

  template <class V, class W> 
  void apply_operator(const alps::SiteOperator& op, site_descriptor s, const V& x, W& y) const
  {
    if (use_bloch_)
      apply_operator(bloch_states,op,s,x,y);
    else
      apply_operator(states,op,s,x,y);
  }
  
  template <class V, class W> 
  void apply_operator(const alps::BondOperator& op, bond_descriptor b, const V&, W&) const;

  template <class V, class W> 
  void apply_operator(const alps::BondOperator& op, site_descriptor s1, site_descriptor s2, const V&, W&) const;

  template <class V, class W> 
  void apply_operator(const boost::multi_array<std::pair<T,bool>,4>& mat, site_descriptor s1, site_descriptor s2, const V& x, W& y) const
  {
    if (use_bloch_)
      apply_operator(bloch_states,mat,s1,s2,x,y);
    else
      apply_operator(states,mat,s1,s2,x,y);
  }

  template <class STATES, class V, class W>
  void apply_operator(const STATES&, const boost::multi_array<std::pair<T,bool>,4>& mat, site_descriptor s1, site_descriptor s2, const V& x, W& y) const;

  template <class V, class W> 
  void apply_operator(const alps::SiteOperator& op, const V&, W&) const;
  
  template <class V, class W> 
  void apply_operator(const alps::BondOperator& op, const V&, W&) const;

  template <class V, class W> 
  void apply_operator(const alps::GlobalOperator& op, const V&, W&) const;

  template <class MM, class OP> 
  MM operator_matrix(const OP& op) const 
  {
    MM m(dimension(),dimension());
    m.clear();
    for (int i=0;i<dimension();++i) {
      m(i,i) +=1.;
      m(i,i) -=1.;
    } // Boost bug workaround
    add_operator_matrix(m,op);
    return m;
  }

 
  template <class MM, class OP> 
  void add_operator_matrix(MM& m, const OP& op) const 
  {
    alps::matrix_as_vector<MM> v(m);
    apply_operator(op,v,v);
  }

 
  template <class MM, class OP, class D> 
  MM operator_matrix(const OP& op, D d) const 
  {
    MM m(dimension(),dimension());
    m.clear();
    for (int i=0;i<dimension();++i) {
      m(i,i) +=1.;
      m(i,i) -=1.;
    } // Boost bug workaround
    add_operator_matrix(m,op,d);
    return m;
  }

  template <class MM, class OP, class D> 
  void add_operator_matrix(MM& m, const OP& op, D d) const 
  {
    alps::matrix_as_vector<MM> v(m);
    apply_operator(op,d,v,v);
  }

  template <class MM, class OP> 
  MM operator_matrix(const OP& op, site_descriptor s1, site_descriptor s2) const 
  {
    MM m(dimension(),dimension());
    m.clear();
    for (int i=0;i<dimension();++i) {
      m(i,i) +=1.;
      m(i,i) -=1.;
    } // Boost bug workaround
    add_operator_matrix(m,op,s1,s2);
    return m;
  }

  template <class MM, class OP> 
  void add_operator_matrix(MM& m, const OP& op, site_descriptor s1, site_descriptor s2) const 
  {
    alps::matrix_as_vector<MM> v(m);
    apply_operator(op,s1,s2,v,v);
  }

  boost::multi_array<T,2> local_matrix(const alps::SiteOperator& op, site_descriptor s) const;
  boost::multi_array<std::pair<T,bool>,4> local_matrix(const alps::BondOperator& op, const bond_descriptor& b) const;
  boost::multi_array<std::pair<T,bool>,4> local_matrix(const alps::BondOperator& op, const site_descriptor& s1, const site_descriptor& s2) const;

  template <class V, class W> 
  void apply_operator(const std::string& name, bond_descriptor b, const V& x, W& y) const
  { apply_operator(get_bond_operator(name),b,x,y); }

  template <class V, class W>
  void apply_operator(const std::string& op, site_descriptor s, const V& x, W& y) const;

  template <class V, class W> 
  void apply_operator(const std::string& name, const V& x, W& y) const;


protected:
  void build() const;
  void build_basis() const;
  bool uses_translation_invariance() const { return use_bloch_;}

  std::vector<std::vector<std::pair<std::string,std::string> > > quantumnumbervalues_;
  mutable basis_states_type states;
  mutable bloch_basis_states_type bloch_states;

private:
  typedef std::pair<std::string,std::string> string_pair;
  typedef std::pair<half_integer_type,half_integer_type> half_integer_pair;
  typedef boost::tuple<half_integer_type,half_integer_type,half_integer_type> half_integer_tuple;
  typedef std::vector<std::pair<string_pair,half_integer_tuple> > QNRangeType;

  void build_subspaces(const std::string&);
  virtual void do_subspace()=0;
  void invalidate() { built_matrix_=built_basis_=false;}
  
  QNRangeType ranges_;
  mutable bool built_matrix_;
  mutable bool built_basis_;
  mutable matrix_type matrix_;
  mutable alps::basis_states_descriptor<short> basis_;
  mutable vector_type total_momentum;
  mutable bool use_bloch_;
  
};

   
template <class T, class M>
void HamiltonianMatrix<T,M>::dostep() 
{
  if (finished()) 
    return;
  build_subspaces(parms["CONSERVED_QUANTUMNUMBERS"]);
  std::vector<half_integer_type> indices(ranges_.size());
  std::vector<std::string> momenta;
  if (parms.value_or_default("TRANSLATION_SYMMETRY",true)) {
    std::vector<vector_type> k = translation_momenta();
    for (typename std::vector<vector_type>::const_iterator it=k.begin();it!=k.end();++it)
      momenta.push_back(alps::write_vector(*it));
  }
  int ik=0;
  bool loop_momenta = ! parms.defined("TOTAL_MOMENTUM");	// Loop over momenta ?
  bool done;
  do { 
    // set QN
    std::vector<std::pair<std::string,std::string> > qns;
    for (int i=0;i<indices.size();++i) {
      parms[ranges_[i].first.second]=ranges_[i].second.get<0>()+indices[i];
      qns.push_back(std::make_pair(
        ranges_[i].first.first,
        boost::lexical_cast<std::string>(ranges_[i].second.get<0>()+indices[i])));
    }

    if (loop_momenta && ik<momenta.size())
      parms["TOTAL_MOMENTUM"]=momenta[ik];
    invalidate();
    if (parms.defined("TOTAL_MOMENTUM")) {
      if(loop_momenta)
        qns.push_back(std::make_pair(std::string("TOTAL_MOMENTUM"),momenta[ik]));
      std::vector<alps::Expression> k;
      alps::read_vector_resize(parms["TOTAL_MOMENTUM"],k);
      alps::ParameterEvaluator eval(parms);
      total_momentum.clear();
      for (int i=0;i<k.size();++i)
        total_momentum.push_back(std::real(k[i].value(eval)));
    }
    use_bloch_=alps::dimension(total_momentum)!=0;
    basis().set_parameters(parms);
    if (dimension()) {
      quantumnumbervalues_.push_back(qns);
      // get spectrum
      do_subspace();
    }
    
    // increment indices
    int j=0;
    ++ik;
    if (ik>=momenta.size() || (! loop_momenta)) {
      ik=0;
      while (j!=indices.size()) {
        indices[j] += ranges_[j].second.get<2>();
        if (ranges_[j].second.get<0>()+indices[j]<=ranges_[j].second.get<1>())
          break;
        indices[j]=0;
        ++j;
      }
    }
    done = (indices.size()==0 ? ik==0 : j==indices.size());
  } while (!done);
  finish();
}


template <class T, class M>
HamiltonianMatrix<T,M>::HamiltonianMatrix(const alps::ProcessList& w,const boost::filesystem::path& fn)
  : alps::scheduler::Task(w,fn), 
    alps::graph_helper<>(parms), 
    alps::model_helper<>(*this,parms),
    built_matrix_(false),
    built_basis_(false),
    use_bloch_(false)
{}

template <class T, class M>
void HamiltonianMatrix<T,M>::build_subspaces(const std::string& quantumnumbers)
{
  // split the string into tokens for the quantum numbers
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" ,");
  tokenizer tokens(quantumnumbers, sep);
  std::vector<std::string> quantumnumber_names;
  std::copy(tokens.begin(),tokens.end(),std::back_inserter(quantumnumber_names));
  // check for unevaluated constraints on these quantum numbers
  std::vector<string_pair> constraints;
  for (basis_descriptor_type::unevaluated_constraints_type::const_iterator  
         it =  basis().unevaluated_constraints().begin();
         it != basis().unevaluated_constraints().end(); ++it) {
    if (std::find(quantumnumber_names.begin(),quantumnumber_names.end(),it->first)!=quantumnumber_names.end())
      constraints.push_back(std::make_pair(it->first,boost::lexical_cast<std::string>(it->second)));
  }
  // get the range for each unevaluated quantum number      
  alps::basis_states_descriptor<short> basis_(basis(),graph());
  for (unsigned int i=0;i<constraints.size();++i) {
    half_integer_tuple qn= boost::make_tuple(half_integer_type(0),half_integer_type(0),half_integer_type(1));
    for (unsigned int s=0;s<basis_.size();++s) {
      unsigned int k=alps::get_quantumnumber_index(constraints[i].first,basis_[s].basis());
      qn.get<0>() += basis_[s].basis()[k].global_min();
      qn.get<1>() += basis_[s].basis()[k].global_max();
      if (basis_[s].basis()[k].global_increment().is_odd())
        qn.get<2>() = 0.5;
    }
    ranges_.push_back(std::make_pair(constraints[i],qn));
    std::cerr << "Quantumnumber " << constraints[i].first << " going from " 
              << qn.get<0>() << " to " << qn.get<1>() << " with increment "
              << qn.get<2>() << "\n";
  }
}


template <class T, class M>
boost::multi_array<T,2> HamiltonianMatrix<T,M>::local_matrix(const alps::SiteOperator& op, site_descriptor s) const
{
  alps::Parameters p(parms);
  if (inhomogeneous_sites()) {
    throw_if_xyz_defined(parms,s); // check whether x, y, or z is set
    p << coordinate_as_parameter(s); // set x, y and z
  }
  return alps::get_matrix(T(),op,model().basis().site_basis(site_type(s)),p);
}

template <class T, class M>
boost::multi_array<std::pair<T,bool>,4> HamiltonianMatrix<T,M>::local_matrix(const alps::BondOperator& op, const bond_descriptor& b) const
{
  unsigned int stype1 = site_type(source(b));
  unsigned int stype2 = site_type(target(b));
  alps::Parameters p(parms);
  if (inhomogeneous_bonds()) {
    throw_if_xyz_defined(parms,b); // check whether x, y, or z is set
    p << coordinate_as_parameter(b); // set x, y and z
  }
  return alps::get_fermionic_matrix(T(),op,model().basis().site_basis(stype1),
                                      model().basis().site_basis(stype2),p);  
}

template <class T, class M>
boost::multi_array<std::pair<T,bool>,4> HamiltonianMatrix<T,M>::local_matrix(const alps::BondOperator& op, const site_descriptor& s1, const site_descriptor& s2) const
{
  unsigned int stype1 = site_type(s1);
  unsigned int stype2 = site_type(s2);
  return alps::get_fermionic_matrix(T(),op,model().basis().site_basis(stype1),
                                      model().basis().site_basis(stype2),parms);  
}

template <class T, class M> template <class STATES, class V, class W>
void HamiltonianMatrix<T,M>::apply_operator(const STATES& states, const alps::SiteOperator& op, site_descriptor s, const V& x, W& y) const
{
  boost::multi_array<T,2> mat = local_matrix(op,s);
  for (int i=0;i<dimension();++i) {           // loop basis states
    state_type state=states[i];               // get source state
    int is=state[s];                          // get site basis index
    for (int js=0;js<basis_[s].size();++js) { // loop over target site states
      T val=mat[is][js];                      // get matrix element
      if (alps::is_nonzero(val)) {            // if matrix element is nonzero
        state_type newstate=state;            // prepare target state
        newstate[s]=js;                       // build target state
        std::complex<double> phase;
        std::size_t j;       
        boost::tie(j,phase) = states.index_and_phase(newstate);       // lookup target state
        if (j<dimension()) {
          val *= optional_real<value_type>::real(phase)  * states.normalization(j)/states.normalization(i);
          y[i] += val*x[j];                   // set matrix element
        }
      }
    }
  }
}

template <class T, class M> template <class V, class W>
void HamiltonianMatrix<T,M>::apply_operator(const alps::BondOperator& op, bond_descriptor b, const V& x, W& y) const
{
  apply_operator(local_matrix(op,b),source(b),target(b),x,y);
}  


template <class T, class M> template <class V, class W>
void HamiltonianMatrix<T,M>::apply_operator(const alps::BondOperator& op, site_descriptor s1, site_descriptor s2, const V& x, W& y) const
{
  apply_operator(local_matrix(op,s1,s2),s1,s2,x,y);
}


template <class T, class M> template <class STATES, class V, class W>
void HamiltonianMatrix<T,M>::apply_operator(const STATES& states, const boost::multi_array<std::pair<T,bool>,4>& mat, site_descriptor s1, site_descriptor s2, const V& x, W& y) const
{
  for (std::size_t i=0;i<dimension();++i) {             // loop over source states
    state_type state=states[i];                   // get source state
    int is1=state[s1];                            // get source site states
    int is2=state[s2];
    for (int js1=0;js1<basis_[s1].size();++js1) { // loop over target site states
      for (int js2=0;js2<basis_[s2].size();++js2) {
        T val=mat[is1][is2][js1][js2].first;      // get matrix element
        if (alps::is_nonzero(val)) {              // if nonzero matrix element
          state_type newstate=state;              // prepare target state
          newstate[s1]=js1;                       // build target state
          newstate[s2]=js2;
          std::complex<double> phase;
          std::size_t j;       
          boost::tie(j,phase) = states.index_and_phase(newstate);       // lookup target state
          if (j<dimension()) {
            if (mat[is1][is2][js1][js2].second) {
              // calculate fermionic sign
              bool f=(s2>=s1);
              int start = std::min(s1,s2);
              int end = std::max(s1,s2);

              for (int i=start;i<end;++i)
                if (is_fermionic(model().basis().site_basis(site_type(i)),basis_[i][state[i]]))
                  f=!f;
              if (f)
                val=-val;
            }
            val *= optional_real<value_type>::real(phase) * states.normalization(j)/states.normalization(i);
            y[i] += val*x[j];                   // set matrix element
          }
        }
      }
    }
  }
}



template <class T, class M> template <class V, class W>
void HamiltonianMatrix<T,M>::apply_operator(const std::string& op, site_descriptor s, const V& x, W& y) const
{
  if (has_site_operator(op))
    apply_operator(get_site_operator(op),s,x,y);
  else
    apply_operator(alps::SiteOperator(op+"(i)","i"),s,x,y);
}


template <class T, class M> template <class V, class W>
void HamiltonianMatrix<T,M>::apply_operator(const std::string& name, const V& x, W& y) const
{ 
  if (has_site_operator(name))
    apply_operator(get_site_operator(name),x,y); 
  else if (has_bond_operator(name))
    apply_operator(get_bond_operator(name),x,y);
  else if (has_global_operator(name))
    apply_operator(get_global_operator(name),x,y);
  else // assume site operator
    apply_operator(alps::SiteOperator(name+"(i)","i"),x,y);
}


template <class T, class M> template <class V, class W>
void HamiltonianMatrix<T,M>::apply_operator(const alps::SiteOperator& op, const V& x, W& y) const
{
  for (site_iterator it=sites().first; it!=sites().second ; ++it) 
    apply_operator(op,*it,x,y);
}
  
template <class T, class M> template <class V, class W>
void HamiltonianMatrix<T,M>::apply_operator(const alps::BondOperator& op, const V& x, W& y) const
{
  for (bond_iterator it=bonds().first; it!=bonds().second ; ++it) 
    apply_operator(op,*it,x,y);
}

template <class T, class M> template <class V, class W>
void HamiltonianMatrix<T,M>::apply_operator(const alps::GlobalOperator& op, const V& x, W& y) const
{
  // loop over sites
  for (site_iterator it=sites().first; it!=sites().second ; ++it)
    apply_operator(op.site_term(site_type(*it)),*it,x,y);

  // loop over bonds
  for (bond_iterator it=bonds().first; it!=bonds().second ; ++it) 
      apply_operator(op.bond_term(bond_type(*it)),*it,x,y);
}




template <class T, class M>
void HamiltonianMatrix<T,M>::build_basis() const
{
  boost::timer t;
  // create basis set
  //std::cerr << "Building basis\n";
  basis_ = alps::basis_states_descriptor<short>(basis(),graph());
  if (use_bloch_)
    bloch_states = bloch_basis_states_type(basis_,translations(total_momentum));
  else
    states = basis_states_type(basis_);

  built_basis_ = true;
  //std::cerr << "Time to build basis: " << t.elapsed() << "\n";
}    

template <class T, class M>
void HamiltonianMatrix<T,M>::build() const
{
  // build matrix
  alps::Disorder::seed(parms.value_or_default("DISORDER_SEED",0));
  boost::timer t;
  //std::cerr << "Building matrix\n";
  matrix_ = matrix_type(dimension(),dimension());
  matrix_.clear();
  built_matrix_ = true;
  add_operator_matrix(matrix_,model());
  // counter ublas bug
  for (int i=0;i<dimension();++i)
    matrix_(i,i) += 1.;
  for (int i=0;i<dimension();++i)
    matrix_(i,i) -= 1.;
    
  //std::cerr << "Time to build matrix: " << t.elapsed() << "\n";
  if (this->parms.value_or_default("CHECK_SYMMETRY",false)) {
    std::cerr << "Checking symmetry\n";
    for (int i=0;i<dimension();++i)
      for (int j=0;j<i;++j)
        if (std::abs(static_cast<value_type>(matrix_(i,j))-::conj(static_cast<value_type>(matrix_(j,i)))) > 1e-10) {
          std::cerr << "Symmetry problem: " << i << " " << j << " " << static_cast<value_type>(matrix_(i,j)) << " " << ::conj(static_cast<value_type>(matrix_(j,i))) << "\n";
          std::cout << basis() << "\n";
          std::cout << states << "\n";
          std::cout << bloch_states << "\n";
          std::cout << matrix() << "\n";
          std::abort();
        }
    std::cerr << "Checked symmetry\n";
  }
}
