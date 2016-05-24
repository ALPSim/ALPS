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

#include <alps/model/hamiltonian_matrix.hpp>
#include <alps/lattice.h>
#include <alps/scheduler/diag.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tokenizer.hpp>
#include <cmath>
#include <cstddef>


template <class T, class M>
class DiagMatrix : public alps::scheduler::DiagTask<T>, public alps::hamiltonian_matrix<M>
{
public:
  typedef T value_type;
  typedef alps::model_helper<>::half_integer_type half_integer_type;
  typedef typename alps::hamiltonian_matrix<M>::vector_type vector_type;
  typedef typename alps::hamiltonian_matrix<M>::basis_descriptor_type basis_descriptor_type;
  typedef typename alps::hamiltonian_matrix<M>::site_iterator site_iterator;
  typedef typename alps::hamiltonian_matrix<M>::site_descriptor site_descriptor;
  typedef typename alps::hamiltonian_matrix<M>::bond_iterator bond_iterator;

  typedef boost::numeric::ublas::mapped_vector_of_mapped_vector<T, boost::numeric::ublas::row_major>  operator_matrix_type;

  DiagMatrix (const alps::ProcessList& , const boost::filesystem::path&,bool delay_construct=false);

  void dostep();

  void perform_measurements();
  
  std::size_t dimension() const { return this->alps::hamiltonian_matrix<M>::dimension();}
  
private:
  typedef std::pair<std::string,std::string> string_pair;
  typedef std::pair<half_integer_type,half_integer_type> half_integer_pair;
  typedef boost::tuple<half_integer_type,half_integer_type,half_integer_type> half_integer_tuple;
  typedef std::vector<std::pair<string_pair,half_integer_tuple> > QNRangeType;

  void build_subspaces(const std::string&);
  
  virtual void do_subspace() =0;
  
  virtual std::vector<value_type> calculate(operator_matrix_type const&) const =0;
  
  template <class Op> 
  std::vector<value_type> calculate(Op const& op) const;
  
  template <class Op, class D> 
  std::vector<value_type> calculate(Op const& op, D const&) const;
  
  template <class Op, class D> 
  std::vector<value_type> calculate(Op const& op, std::pair<D,D>  const&) const;
  
  void print() const;
  virtual void print_eigenvectors(std::ostream& os) const=0;

  std::vector<unsigned int> multiplicities_;    
  QNRangeType ranges_;
};


template <class T, class M>
DiagMatrix<T,M>::DiagMatrix(const alps::ProcessList& where , const boost::filesystem::path& p, bool delay_construct) 
    : alps::scheduler::DiagTask<T>(where,p,delay_construct)
    , alps::hamiltonian_matrix<M>(this->get_parameters())
{ 
  if (this->calc_averages())
    multiplicities_ = this->distance_multiplicities();
}

template <class T, class M>
void DiagMatrix<T,M>::dostep() 
{
  if (this->finished()) 
    return;
  build_subspaces(this->alps::scheduler::Task::parms["CONSERVED_QUANTUMNUMBERS"]);
  std::vector<half_integer_type> indices(ranges_.size());
  std::vector<std::string> momenta;
  if (this->alps::scheduler::Task::parms.value_or_default("TRANSLATION_SYMMETRY",true)) {
    std::vector<vector_type> k = this->translation_momenta();
    for (typename std::vector<vector_type>::const_iterator it=k.begin();it!=k.end();++it)
      momenta.push_back(alps::write_vector(*it));
  }
  unsigned ik=0;
  bool loop_momenta = ! this->alps::scheduler::Task::parms.defined("TOTAL_MOMENTUM");    // Loop over momenta ?
  bool done;
  do { 
    // set QN
    std::vector<std::pair<std::string,std::string> > qns;
    for (unsigned i=0;i<indices.size();++i) {
      this->alps::scheduler::Task::parms[ranges_[i].first.second]=boost::get<0>(ranges_[i].second)+indices[i];
      qns.push_back(std::make_pair(
        ranges_[i].first.first,
        boost::lexical_cast<std::string>(boost::get<0>(ranges_[i].second)+indices[i])));
    }

    if (loop_momenta && ik<momenta.size())
      this->alps::scheduler::Task::parms["TOTAL_MOMENTUM"]=momenta[ik];
    if (this->alps::scheduler::Task::parms.defined("TOTAL_MOMENTUM") && loop_momenta)
      qns.push_back(std::make_pair(std::string("TOTAL_MOMENTUM"),momenta[ik]));
    this->set_parameters(this->alps::scheduler::Task::parms);
    if (this->dimension()) {
      this->quantumnumbervalues_.push_back(qns);
      // get spectrum
      do_subspace();
    }
    
    // increment indices
    int j=0;
    ++ik;
    if (ik>=momenta.size() || (! loop_momenta)) {
      ik=0;
      while (j!=indices.size()) {
        indices[j] += boost::get<2>(ranges_[j].second);
        if (boost::get<0>(ranges_[j].second)+indices[j]<=boost::get<1>(ranges_[j].second))
          break;
        indices[j]=0;
        ++j;
      }
    }
    done = (indices.size()==0 ? ik==0 : j==indices.size());
  } while (!done);
  this->finish();
}



template <class T, class M>
void DiagMatrix<T,M>::build_subspaces(const std::string& quantumnumbers)
{
  // split the string into tokens for the quantum numbers
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" ,");
  tokenizer tokens(quantumnumbers, sep);
  std::vector<std::string> quantumnumber_names;
  std::copy(tokens.begin(),tokens.end(),std::back_inserter(quantumnumber_names));
  // check for unevaluated constraints on these quantum numbers
  std::vector<string_pair> constraints;
  for (typename basis_descriptor_type::unevaluated_constraints_type::const_iterator  
         it =  this->basis().unevaluated_constraints().begin();
         it != this->basis().unevaluated_constraints().end(); ++it) {
    if (std::find(quantumnumber_names.begin(),quantumnumber_names.end(),it->first)!=quantumnumber_names.end())
      constraints.push_back(std::make_pair(it->first,boost::lexical_cast<std::string>(it->second)));
  }
  // get the range for each unevaluated quantum number      
  alps::basis_states_descriptor<short> basis_(this->basis(),this->graph());
  for (unsigned int i=0;i<constraints.size();++i) {
    half_integer_tuple qn= boost::make_tuple(half_integer_type(0),half_integer_type(0),half_integer_type(1));
    for (unsigned int s=0;s<basis_.size();++s) {
      unsigned int k=alps::get_quantumnumber_index(constraints[i].first,basis_[s].basis());
      boost::get<0>(qn) += basis_[s].basis()[k].global_min();
      boost::get<1>(qn) += basis_[s].basis()[k].global_max();
      if (basis_[s].basis()[k].global_increment().is_odd())
        boost::get<2>(qn) = 0.5;
    }
    ranges_.push_back(std::make_pair(constraints[i],qn));
    std::cerr << "Quantumnumber " << constraints[i].first << " going from " 
              << boost::get<0>(qn) << " to " << boost::get<1>(qn) << " with increment "
              << boost::get<2>(qn) << "\n";
  }
}


template <class T, class M>
void DiagMatrix<T,M>::print() const
{
  std::cout << "------------------------------------------------------------------------------------------\n";
  std::cout << "Eigenvectors for the sector with parameters";
  std::cout << this->get_parameters();
  std::cout << "\n\nBasis:\n";
  this->print_basis(std::cout);
  std::cout << "\n\nVectors:\n";
  this->print_eigenvectors(std::cout);
}


template <class T, class M>
void DiagMatrix<T,M>::perform_measurements()
{
  typedef std::pair<std::string,std::string> string_pair;
  
  if (this->print_vectors())
    print();
  
  alps::EigenvectorMeasurements<value_type> meas(*this);
  if (this->calc_averages()) {

    // perform measurements

    BOOST_FOREACH (string_pair const& ex, this->average_expressions) {
      //std::cerr << "Evaluating " << ex.first << "\n";
      alps::SiteOperator op(ex.second+"(i)/"+ boost::lexical_cast<std::string>(this->num_sites()),"i");
      this->substitute_operators(op,this->alps::scheduler::Task::parms);
      meas.average_values[ex.first] = calculate(op);
    }

    // calculate local measurements
    if (!this->uses_translation_invariance()) {
      BOOST_FOREACH (string_pair const& ex, this->local_expressions) {
        //std::cerr << "Evaluating " << ex.first << "\n";
        if (this->has_bond_operator(ex.second)) {
          for (bond_iterator bit=this->bonds().first; bit!=this->bonds().second;++bit) {
            std::vector<value_type> av = calculate(ex.second,*bit);
            meas.local_values[ex.first].resize(av.size());
            for (unsigned i=0;i<av.size();++i)
              meas.local_values[ex.first][i].push_back(av[i]);
          }
        }
        else {
          for (site_iterator sit=this->sites().first; sit!=this->sites().second;++sit) {
            std::vector<value_type> av = calculate(ex.second,*sit);
            meas.local_values[ex.first].resize(av.size());
            for (unsigned i=0;i<av.size();++i)
              meas.local_values[ex.first][i].push_back(av[i]);
          }
        }
      }
    }
    
    // calculate correlations
    typedef std::pair<std::string,std::pair<std::string,std::string> > string_string_pair_pair;
    BOOST_FOREACH (string_string_pair_pair const& ex, this->correlation_expressions) {
      //std::cerr << "Evaluating " << ex.first << "\n";
      std::vector<bool> done(this->num_distances(),false);
      for (site_iterator sit1=this->sites().first; sit1!=this->sites().second ; ++sit1) {
        for (site_iterator sit2=this->sites().first; sit2!=this->sites().second ; ++sit2) {
          std::size_t d = alps::scheduler::DiagTask<T>::distance(*sit1,*sit2);
          if (!done[d] || this->uses_translation_invariance()) {
            std::vector<value_type> av;
            if (*sit1 == *sit2) {
              alps::SiteOperator op(ex.second.first+"(i)*"+ex.second.second+"(i)","i");
              this->substitute_operators(op,this->alps::scheduler::Task::parms);
              av = calculate(op,*sit1);
            }
            else {
              alps::BondOperator op(ex.second.first+"(i)*"+ex.second.second+"(j)","i","j");
              this->substitute_operators(op,this->alps::scheduler::Task::parms);
              av = calculate(op,std::make_pair(*sit1,*sit2));
            }
            meas.correlation_values[ex.first].resize(av.size());
            for (unsigned i=0;i<av.size();++i) {
              if (meas.correlation_values[ex.first][i].size()<=d)
                meas.correlation_values[ex.first][i].resize(d+1);
              meas.correlation_values[ex.first][i][d] += (av[i] / 
                   static_cast<double>(this->uses_translation_invariance() ? this->multiplicities_[d] : 1.));
            }
            done[d]=true;
          }
        }
      }
    }



    // calculate structure factor
    BOOST_FOREACH (string_string_pair_pair const& ex, this->structurefactor_expressions) {
      //std::cerr << "Evaluating " << ex.first << "\n";
      boost::multi_array<std::vector<value_type>,2> corrs(boost::extents[this->num_sites()][this->num_sites()]);
      for (site_iterator sit1=this->sites().first; sit1!=this->sites().second ; ++sit1)
        for (site_iterator sit2=this->sites().first; sit2!=this->sites().second ; ++sit2)
          if (*sit1 == *sit2) {
            alps::SiteOperator op(ex.second.first+"(i)*"+ex.second.second+"(i)","i");
            this->substitute_operators(op,this->alps::scheduler::Task::parms);
            corrs[*sit1][*sit2] = calculate(op,*sit1);
          }
          else {
            alps::BondOperator op(ex.second.first+"(i)*"+ex.second.second+"(j)","i","j");
            this->substitute_operators(op,this->alps::scheduler::Task::parms);
              corrs[*sit1][*sit2] = calculate(op,std::make_pair(*sit1,*sit2));
          }
      
      // do Fourier-transformed emasurements
      for (typename alps::graph_helper<>::momentum_iterator mit=this->momenta().first; mit != this->momenta().second; ++mit) {
        std::vector<value_type> av;
        av.clear();
        for (site_iterator sit1=this->sites().first; sit1!=this->sites().second ; ++sit1)
          for (site_iterator sit2=this->sites().first; sit2!=this->sites().second ; ++sit2) {
          if (av.size() < corrs[*sit1][*sit2].size())
            av.resize(corrs[*sit1][*sit2].size());
          for (unsigned i=0;i<corrs[*sit1][*sit2].size();++i)
          {
            double phase1 = alps::numeric::scalar_product(this->momentum(*mit), alps::scheduler::DiagTask<T>::coordinate(*sit1));
            double phase2 = alps::numeric::scalar_product(this->momentum(*mit), alps::scheduler::DiagTask<T>::coordinate(*sit2));
            av[i] += std::real(corrs[*sit1][*sit2][i]
                      * std::conj( std::complex<double>(std::cos(phase1), std::sin(phase1)) )
                      * std::complex<double>(std::cos(phase2), std::sin(phase2))
                     ) / static_cast<double>(this->num_sites());
          }
        }
        if (meas.structurefactor_values[ex.first].size() < av.size())
          meas.structurefactor_values[ex.first].resize(av.size());
        for (unsigned i=0;i<av.size();++i)
          meas.structurefactor_values[ex.first][i].push_back(av[i]); 
      }
    }


  }
  this->measurements_.push_back(meas);
}


template <class T, class M>
template <class Op>
std::vector<T> DiagMatrix<T,M>::calculate(Op const& op) const
{
  
  return calculate(this->template operator_matrix<operator_matrix_type>(op));
}

template <class T, class M>
template <class Op, class D>
std::vector<T> DiagMatrix<T,M>::calculate(Op const& op, D const& s) const
{
  return calculate(this->template operator_matrix<operator_matrix_type>(op,s));
}

template <class T, class M>
template <class Op, class D>
std::vector<T> DiagMatrix<T,M>::calculate(Op const& op, std::pair<D,D> const& s) const
{
  return calculate(this->template operator_matrix<operator_matrix_type>(op,s.first,s.second));
}


