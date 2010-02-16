/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>
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

#include "matrix.h"
#include <alps/scheduler/measurement_operators.h>
#include <boost/regex.hpp> 

#ifdef ALPS_HAVE_HDF5
#include <alps/hdf5.hpp>
#endif


template <class T, class M>
class DiagMatrix : public HamiltonianMatrix<T,M>, public alps::MeasurementOperators
{
public:
  typedef HamiltonianMatrix<T,M> super_type;
  typedef T value_type;
  typedef typename super_type::magnitude_type magnitude_type;
  typedef typename super_type::matrix_type matrix_type;
  typedef typename super_type::site_iterator site_iterator;
  typedef typename super_type::site_descriptor site_descriptor;
  typedef typename super_type::bond_iterator bond_iterator;
  typedef boost::numeric::ublas::vector<value_type> vector_type;
  typedef boost::numeric::ublas::vector<magnitude_type> mag_vector_type;
  typedef typename super_type::half_integer_type half_integer_type;
  typedef boost::numeric::ublas::mapped_vector_of_mapped_vector<T, boost::numeric::ublas::row_major>  operator_matrix_type;
  
  DiagMatrix (const alps::ProcessList& where , const boost::filesystem::path& p);

#ifdef ALPS_HAVE_HDF5
  void serialize(alps::hdf5::oarchive &) const;
  void serialize(alps::hdf5::iarchive &);
#endif

protected:
  void write_xml_body(alps::oxstream&, const boost::filesystem::path&,bool) const;
  void handle_tag(std::istream& infile, const alps::XMLTag& tag); 

  std::vector<mag_vector_type> eigenvalues_;
  std::vector<unsigned int> multiplicities_;
  std::vector<alps::EigenvectorMeasurements<value_type> > measurements_;
  void perform_measurements();
private:
  virtual std::vector<value_type> calculate(operator_matrix_type const&) const=0;
  
  template <class Op> 
  std::vector<value_type> calculate(Op const& op) const;
  
  template <class Op, class D> 
  std::vector<value_type> calculate(Op const& op, D const&) const;
  
  template <class Op, class D> 
  std::vector<value_type> calculate(Op const& op, std::pair<D,D>  const&) const;
  
  bool read_hdf5_;
};

template <class T, class M>
DiagMatrix<T,M>::DiagMatrix(const alps::ProcessList& where , const boost::filesystem::path& p) 
    : super_type(where,p)
    , alps::MeasurementOperators(this->parms)
    , read_hdf5_(false)
{ 
  if (this->calc_averages())
    multiplicities_ = this->distance_multiplicities();
}


#ifdef ALPS_HAVE_HDF5
template <class T, class M>
void DiagMatrix<T,M>::serialize(alps::hdf5::iarchive & ar) {
  alps::scheduler::Task::serialize(ar);
  if (ar.is_group("/spectrum/sectors")) {
      std::vector<std::string> list = ar.list_children("/spectrum/sectors");
      measurements_.resize(list.size(),alps::EigenvectorMeasurements<value_type>(*this));
      for (unsigned i=0; i<list.size();++i) {
        std::string sectorpath = "/spectrum/sectors/"+list[i];
        
        // read quantum numbers
        std::vector<std::pair<std::string,std::string> > qnvals;
        if (ar.is_group(sectorpath+"/quantumnumbers")) {
          std::vector<std::string> qnlist = ar.list_children(sectorpath+"/quantumnumbers");
          for (std::vector<std::string>::const_iterator it = qnlist.begin(); it != qnlist.end(); ++it) {
              std::string v;
              ar >> alps::make_pvp(sectorpath+"/quantumnumbers/"+*it, v);
              qnvals.push_back(std::make_pair(*it,v));
          }
        }
        this->quantumnumbervalues_.push_back(qnvals);

        // read energies
        if (ar.is_data(sectorpath+"/energies")) {
          mag_vector_type evals_vector;
          ar >> alps::make_pvp(sectorpath+"/energies", evals_vector);
          eigenvalues_.push_back(evals_vector);
        }
        // read measurements
          ar >> alps::make_pvp(sectorpath,measurements_[i]);
      }
  }
  this->read_hdf5_ = true; // skip XML, once all is being read
}

template <class T, class M>
void DiagMatrix<T,M>::serialize(alps::hdf5::oarchive & ar) const {
  alps::scheduler::Task::serialize(ar);
  for (int i=0;i<eigenvalues_.size();++i) {
    std::string sectorpath = "/spectrum/sectors/" + boost::lexical_cast<std::string>(i);
    for (int j=0;j<this->quantumnumbervalues_[i].size();++j)
      ar << alps::make_pvp(sectorpath + "/quantumnumbers/" + this->quantumnumbervalues_[i][j].first,
                     this->quantumnumbervalues_[i][j].second);
      ar << alps::make_pvp(sectorpath + "/energies",eigenvalues_[i]);
    if (calc_averages() || this->parms.value_or_default("MEASURE_ENERGY",true))
      ar << alps::make_pvp(sectorpath,measurements_[i]);
  }
}
#endif


template <class T, class M>
void DiagMatrix<T,M>::write_xml_body(alps::oxstream& out, const boost::filesystem::path& name,bool writeallxml) const
{
  if (writeallxml) {
    for (int i=0;i<eigenvalues_.size();++i) {
      int num_eigenvalues = std::min(int(this->parms.value_or_default("NUMBER_EIGENVALUES",
                  eigenvalues_[i].size())),int(eigenvalues_[i].size()));
      out << alps::start_tag("EIGENVALUES") << alps::attribute("number",num_eigenvalues);
      for (int j=0;j<this->quantumnumbervalues_[i].size();++j)
        out << alps::start_tag("QUANTUMNUMBER") << alps::attribute("name",this->quantumnumbervalues_[i][j].first)
            << alps::attribute("value",this->quantumnumbervalues_[i][j].second) << alps::end_tag("QUANTUMNUMBER");
      for (int j=0;j<num_eigenvalues;++j)
        out << eigenvalues_[i][j] << "\n";
      out << alps::end_tag("EIGENVALUES");
    }

    if (calc_averages() || this->parms.value_or_default("MEASURE_ENERGY",true)) {
      for (int i=0;i<eigenvalues_.size();++i) {
         int num_eigenvalues = std::min(int(this->parms.value_or_default("NUMBER_EIGENVALUES",
                  eigenvalues_[i].size())),int(eigenvalues_[i].size()));
        out << alps::start_tag("EIGENSTATES") << alps::attribute("number",num_eigenvalues);
        for (int j=0;j<this->quantumnumbervalues_[i].size();++j)
          out << alps::start_tag("QUANTUMNUMBER") << alps::attribute("name",this->quantumnumbervalues_[i][j].first)
              << alps::attribute("value",this->quantumnumbervalues_[i][j].second) << alps::end_tag("QUANTUMNUMBER");
        for (int j=0;j<num_eigenvalues;++j) {
          out << alps::start_tag("EIGENSTATE") << alps::attribute("number",j);
          measurements_[i].write_xml_one_vector(out,name,j);
          out << alps::end_tag("EIGENSTATE");   
        }
        out << alps::end_tag("EIGENSTATES");
      }
    }
  }
}
   
template <class T, class M>
void DiagMatrix<T,M>::handle_tag(std::istream& infile, const alps::XMLTag& intag) 
{
  alps::XMLTag tag(intag);
  
  // we don't need to read the XML file if the HDF-5 file has already been read
  if (this->read_hdf5_) {
     skip_element(infile,tag);
    return;
  }
  
  if (intag.type==alps::XMLTag::SINGLE)
    return;
  if (intag.name=="EIGENVALUES") {
    // we don't need to read the XML file if the HDF-5 file has already been read

    std::vector<std::pair<std::string,std::string> > qnvals;
    std::vector<magnitude_type> evals;
    char c;
    infile >> c;
    while (c=='<' && infile) {
      infile.putback(c);
      tag=alps::parse_tag(infile);
      if (tag.name=="QUANTUMNUMBER") {
        qnvals.push_back(std::make_pair(tag.attributes["name"],tag.attributes["value"]));
      }
      else if (tag.name=="/EIGENVALUES")
        return;
      alps::skip_element(infile,tag);
      infile >> c;
    }
    do {
      infile.putback(c);
      magnitude_type ev;
      infile >> ev >> c;
      evals.push_back(ev);
    } while (c!='<' && infile);
    infile.putback(c);
    tag=alps::parse_tag(infile);
    if (tag.name!="/EIGENVALUES")
      boost::throw_exception(std::runtime_error("Encountered unexpected tag " + tag.name 
                                  + " while pasrsing <EIGENVALUES>\n"));
    mag_vector_type evals_vector(evals.size());
    std::copy(evals.begin(),evals.end(),evals_vector.begin());
    eigenvalues_.push_back(evals_vector);
    this->quantumnumbervalues_.push_back(qnvals);
  }
  else if (intag.name=="EIGENSTATES") {
    measurements_.push_back(alps::EigenvectorMeasurements<value_type>(*this));
    std::vector<std::pair<std::string,half_integer_type> > qnvals;
    alps::XMLTag tag=alps::parse_tag(infile);
    while (tag.name !="/EIGENSTATES") {
      if (tag.name=="QUANTUMNUMBER")
        alps::skip_element(infile,tag);
      else if (tag.name=="EIGENSTATE") {
        if (tag.type != alps::XMLTag::SINGLE) {
          tag = alps::parse_tag(infile);
          tag = measurements_.rbegin()->handle_tag(infile,tag);
          if (tag.name!="/EIGENSTATE") 
            boost::throw_exception(std::runtime_error("unexpected element " + tag.name + " inside <EIGENSTATE>"));
        }
      }
      tag=alps::parse_tag(infile);
    }
  }    
  else
    alps::skip_element(infile,intag);
    
}


template <class T, class M>
void DiagMatrix<T,M>::perform_measurements()
{
  typedef std::pair<std::string,std::string> string_pair;
  
  alps::EigenvectorMeasurements<value_type> meas(*this);
  if (this->calc_averages()) {

    // perform measurements

    BOOST_FOREACH (string_pair const& ex, this->average_expressions) {
      //std::cerr << "Evaluating " << ex.first << "\n";
      meas.average_values[ex.first] = calculate(ex.second);
    }

    // calculate local measurements
    if (!this->uses_translation_invariance()) {
      BOOST_FOREACH (string_pair const& ex, this->local_expressions) {
        //std::cerr << "Evaluating " << ex.first << "\n";
        if (this->has_bond_operator(ex.second)) {
          for (bond_iterator bit=this->bonds().first; bit!=this->bonds().second;++bit) {
            std::vector<value_type> av = calculate(ex.second,*bit);
            meas.local_values[ex.first].resize(av.size());
            for (int i=0;i<av.size();++i)
              meas.local_values[ex.first][i].push_back(av[i]);
          }
        }
        else {
          for (site_iterator sit=this->sites().first; sit!=this->sites().second;++sit) {
            std::vector<value_type> av = calculate(ex.second,*sit);
            meas.local_values[ex.first].resize(av.size());
            for (int i=0;i<av.size();++i)
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
          std::size_t d = distance(*sit1,*sit2);
          if (!done[d] || this->uses_translation_invariance()) {
            std::vector<value_type> av;
            if (*sit1 == *sit2) {
              alps::SiteOperator op(ex.second.first+"(i)*"+ex.second.second+"(i)","i");
              this->substitute_operators(op,this->parms);
              av = calculate(op,*sit1);
            }
            else {
              alps::BondOperator op(ex.second.first+"(i)*"+ex.second.second+"(j)","i","j");
              this->substitute_operators(op,this->parms);
              av = calculate(op,std::make_pair(*sit1,*sit2));
            }
            meas.correlation_values[ex.first].resize(av.size());
            for (int i=0;i<av.size();++i) {
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
            this->substitute_operators(op,this->parms);
            corrs[*sit1][*sit2] = calculate(op,*sit1);
          }
          else {
            alps::BondOperator op(ex.second.first+"(i)*"+ex.second.second+"(j)","i","j");
            this->substitute_operators(op,this->parms);
              corrs[*sit1][*sit2] = calculate(op,std::make_pair(*sit1,*sit2));
          }
      
      // do Fourier-transformed emasurements
      for (typename super_type::momentum_iterator mit=this->momenta().first; mit != this->momenta().second; ++mit) {
        std::vector<value_type> av;
        av.clear();
        for (site_iterator sit1=this->sites().first; sit1!=this->sites().second ; ++sit1)
          for (site_iterator sit2=this->sites().first; sit2!=this->sites().second ; ++sit2) {
          if (av.size() < corrs[*sit1][*sit2].size())
            av.resize(corrs[*sit1][*sit2].size());
          for (int i=0;i<corrs[*sit1][*sit2].size();++i)
            av[i] += std::real(corrs[*sit1][*sit2][i]
                      *std::conj(mit.phase(coordinate(*sit1)))
                      *mit.phase(coordinate(*sit2)))/static_cast<double>(this->num_sites());
        }
        if (meas.structurefactor_values[ex.first].size() < av.size())
          meas.structurefactor_values[ex.first].resize(av.size());
        for (int i=0;i<av.size();++i)
          meas.structurefactor_values[ex.first][i].push_back(av[i]); 
      }
    }


  }
  measurements_.push_back(meas);
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
