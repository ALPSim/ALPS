/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_SCHEDULER_DIAG_HPP
#define ALPS_SCHEDULER_DIAG_HPP

#include <alps/type_traits/norm_type.hpp>
#include <alps/scheduler/measurement_operators.h>
#include <alps/model/model_helper.h>
#include <alps/lattice/graph_helper.h>
#include <alps/scheduler/task.h>
#include <alps/params.hpp>
#include <alps/ngs/make_parameters_from_deprecated.hpp>

#include <alps/numeric/matrix/vector.hpp>
#include <alps/hdf5/numeric_vector.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/type_traits/is_complex.hpp>

#include <sstream>
#include <stdexcept>

namespace alps { namespace scheduler {

template <class T, class G = typename graph_helper<>::graph_type>
class DiagTask : public scheduler::Task, public graph_helper<G>,  public model_helper<>,  public MeasurementOperators
{
public:
  typedef T value_type;
  typedef typename norm_type<value_type>::type magnitude_type;
  typedef alps::numeric::vector<value_type> vector_type;
  typedef alps::numeric::vector<magnitude_type> mag_vector_type;
  typedef model_helper<>::half_integer_type half_integer_type;
  typedef boost::numeric::ublas::mapped_vector_of_mapped_vector<T, boost::numeric::ublas::row_major>  operator_matrix_type;
  
  // define_parameters defines the parameters the task reads: a derived
  // class passes its own, which calls this class's
  DiagTask (const ProcessList& where , const boost::filesystem::path& p,bool delay_construct=false,
            void (*define_parameters)(alps::params&) = &DiagTask::define_parameters);

  /// Defines the parameters DiagTask reads, with their defaults
  static void define_parameters(alps::params& parameters);

  /// A task's parameters, with the parameters `define_parameters` defines parsed
  /** @throws std::invalid_argument naming each of them that does not parse */
  static alps::params read_parameters(const Parameters& parms, void (*define_parameters)(alps::params&));

  void dostep() { boost::throw_exception(std::logic_error("Cannot call dostep on the base class DiagTask")); }

  void save(hdf5::archive &) const;
  void load(hdf5::archive &);
  

protected:
  void write_xml_body(oxstream&, const boost::filesystem::path&,bool) const;
  void handle_tag(std::istream& infile, const XMLTag& tag); 

  std::vector<mag_vector_type> eigenvalues_;
  std::vector<EigenvectorMeasurements<value_type> > measurements_;

  std::vector<std::vector<std::pair<std::string,std::string> > > quantumnumbervalues_;

  // per sector, filled when SAVE_EIGENVECTORS is set: the basis states as
  // row-major [state][site] local-state indices and the eigenvectors as
  // row-major [eigenvector][state] coefficients
  std::vector<std::vector<int> > basis_states_;
  std::vector<std::vector<value_type> > eigenvectors_;
  
  bool calc_vectors() const { return calc_averages() || print_vectors() || save_vectors();}
  bool print_vectors() const { return settings_["PRINT_EIGENVECTORS"].template as<bool>();}
  bool save_vectors() const { return settings_["SAVE_EIGENVECTORS"].template as<bool>();}
  bool measure_energy() const { return settings_["MEASURE_ENERGY"].template as<bool>();}
  // the parameters define_parameters defined, read when the task is constructed
  const alps::params& settings() const { return settings_;}
private:
  // NUMBER_EIGENVALUES, at most all n eigenvalues of a sector
  unsigned number_eigenvalues(std::size_t n) const;

  const alps::params settings_;
  bool skipped_vectors_;
  bool read_hdf5_;
};

template <class T, class G>
DiagTask<T,G>::DiagTask(const ProcessList& where , const boost::filesystem::path& p, bool delay_construct,
                        void (*define_parameters)(alps::params&))
    : scheduler::Task(where,p)
    , graph_helper<G>(this->get_parameters())
    , model_helper<>(this->get_parameters())
    , MeasurementOperators(this->get_parameters())
    , settings_(read_parameters(this->get_parameters(), define_parameters))
    , skipped_vectors_(false)
    , read_hdf5_(false)
{
  if (!delay_construct)
    this->construct();
}

template <class T, class G>
void DiagTask<T,G>::define_parameters(alps::params& parameters) {
  parameters.define<bool>("PRINT_EIGENVECTORS", false, "print the eigenvectors")
            .define<bool>("SAVE_EIGENVECTORS", false, "save the eigenvectors and their basis to the HDF5 output")
            .define<bool>("MEASURE_ENERGY", true, "write the energy of each eigenstate");
  if (parameters.supplied("NUMBER_EIGENVALUES"))
    parameters.define<int>("NUMBER_EIGENVALUES", "number of eigenvalues per sector");
}

template <class T, class G>
alps::params DiagTask<T,G>::read_parameters(const Parameters& parms, void (*define_parameters)(alps::params&)) {
  alps::params parameters = make_parameters_from_deprecated(parms);
  define_parameters(parameters);
  std::ostringstream errors;
  if (parameters.has_missing(errors))
    boost::throw_exception(std::invalid_argument(errors.str()));
  return parameters;
}

template <class T, class G>
unsigned DiagTask<T,G>::number_eigenvalues(std::size_t n) const {
  return settings_.exists("NUMBER_EIGENVALUES")
    ? (std::min)(unsigned(settings_["NUMBER_EIGENVALUES"].template as<int>()), unsigned(n)) : unsigned(n);
}


template <class T, class G>
void DiagTask<T,G>::load(hdf5::archive & ar) {
  scheduler::Task::load(ar);
  if (ar.is_group("/spectrum/sectors")) {
      std::vector<std::string> list = ar.list_children("/spectrum/sectors");
      measurements_.resize(list.size(),EigenvectorMeasurements<value_type>(*this));
      for (unsigned i=0; i<list.size();++i) {
        std::string sectorpath = "/spectrum/sectors/"+list[i];
        
        // read quantum numbers
        std::vector<std::pair<std::string,std::string> > qnvals;
        if (ar.is_group(sectorpath+"/quantumnumbers")) {
          std::vector<std::string> qnlist = ar.list_children(sectorpath+"/quantumnumbers");
          for (std::vector<std::string>::const_iterator it = qnlist.begin(); it != qnlist.end(); ++it) {
              std::string v;
              ar >> make_pvp(sectorpath+"/quantumnumbers/"+*it, v);
              qnvals.push_back(std::make_pair(*it,v));
          }
        }
        this->quantumnumbervalues_.push_back(qnvals);

        // read energies
        if (ar.is_data(sectorpath+"/energies")) {
          mag_vector_type evals_vector;
          ar >> make_pvp(sectorpath+"/energies", evals_vector);
          eigenvalues_.push_back(evals_vector);
        }
        // read measurements
          ar >> make_pvp(sectorpath,measurements_[i]);

        // read eigenvectors; readers that always use double (fulldiag_evaluate,
        // convert2xml) skip the vectors of a complex run instead of failing
        if (ar.is_group(sectorpath+"/eigenvectors")) {
          if (ar.is_complex(sectorpath+"/eigenvectors/values") == is_complex<value_type>::value) {
            basis_states_.resize(i+1);
            eigenvectors_.resize(i+1);
            ar >> make_pvp(sectorpath+"/basis/values", basis_states_[i]);
            ar >> make_pvp(sectorpath+"/eigenvectors/values", eigenvectors_[i]);
          }
          else
            skipped_vectors_ = true;
        }
      }
  }
  std::cerr << eigenvalues_.size() << " sectors\n";
  this->read_hdf5_ = true; // skip XML, once all is being read
}

template <class T, class G>
void DiagTask<T,G>::save(hdf5::archive & ar) const {
  if (skipped_vectors_)
    boost::throw_exception(std::runtime_error("the eigenvectors in this result file are not of this task's "
      "type (real or complex); rewriting the file would drop them"));
  scheduler::Task::save(ar);
  for (unsigned  i=0;i<eigenvalues_.size();++i) {
    std::string sectorpath = "/spectrum/sectors/" + boost::lexical_cast<std::string>(i);
    for (unsigned j=0;j<this->quantumnumbervalues_[i].size();++j)
      ar << make_pvp(sectorpath + "/quantumnumbers/" + this->quantumnumbervalues_[i][j].first,
                     this->quantumnumbervalues_[i][j].second);
    ar << make_pvp(sectorpath + "/energies",eigenvalues_[i]);
    if (calc_averages() || measure_energy())
      ar << make_pvp(sectorpath,measurements_[i]);
    if (i < eigenvectors_.size()) {
      std::vector<std::size_t> shape(2);
      shape[1] = this->num_sites();
      shape[0] = shape[1] ? basis_states_[i].size()/shape[1] : 0;
      ar << make_pvp(sectorpath + "/basis/shape", shape);
      ar << make_pvp(sectorpath + "/basis/values", basis_states_[i]);
      shape[1] = shape[0];
      shape[0] = shape[1] ? eigenvectors_[i].size()/shape[1] : 0;
      ar << make_pvp(sectorpath + "/eigenvectors/shape", shape);
      ar << make_pvp(sectorpath + "/eigenvectors/values", eigenvectors_[i]);
    }
  }
  if (!eigenvectors_.empty()) {
    // the local states that the basis indices refer to, as in the dmrg MPS output
    std::vector<int> types;
    for (unsigned s = 0; s < this->num_sites(); ++s)
      types.push_back(this->site_type(s));
    ar << make_pvp("/spectrum/site_type", types);
    for (int type = 0; type <= alps::maximum_vertex_type(this->graph()); ++type) {
      alps::site_basis<short> b(this->site_basis(type));
      for (std::size_t q = 0; q < this->site_basis(type).size(); ++q) {
        std::vector<double> values;
        for (std::size_t s = 0; s < b.size(); ++s)
          values.push_back(b[s][q].to_double());
        ar << make_pvp("/spectrum/site_basis/" + boost::lexical_cast<std::string>(type) + "/"
                       + hdf5_name_encode(this->site_basis(type)[q].name()), values);
      }
    }
  }
}


template <class T, class G>
void DiagTask<T,G>::write_xml_body(oxstream& out, const boost::filesystem::path& name,bool writeallxml) const
{
  if (writeallxml) {
    for (unsigned i=0;i<eigenvalues_.size();++i) {
      unsigned num_eigenvalues = number_eigenvalues(eigenvalues_[i].size());
      out << start_tag("EIGENVALUES") << attribute("number",num_eigenvalues);
      for (unsigned j=0;j<this->quantumnumbervalues_[i].size();++j)
        out << start_tag("QUANTUMNUMBER") << attribute("name",this->quantumnumbervalues_[i][j].first)
            << attribute("value",this->quantumnumbervalues_[i][j].second) << end_tag("QUANTUMNUMBER");
      for (unsigned j=0;j<num_eigenvalues;++j)
        out << eigenvalues_[i][j] << "\n";
      out << end_tag("EIGENVALUES");
    }

    if (calc_averages() || measure_energy()) {
      for (unsigned i=0;i<eigenvalues_.size();++i) {
         unsigned num_eigenvalues = number_eigenvalues(eigenvalues_[i].size());
        out << start_tag("EIGENSTATES") << attribute("number",num_eigenvalues);
        for (unsigned j=0;j<this->quantumnumbervalues_[i].size();++j)
          out << start_tag("QUANTUMNUMBER") << attribute("name",this->quantumnumbervalues_[i][j].first)
              << attribute("value",this->quantumnumbervalues_[i][j].second) << end_tag("QUANTUMNUMBER");
        for (unsigned j=0;j<num_eigenvalues;++j) {
          out << start_tag("EIGENSTATE") << attribute("number",j);
          measurements_[i].write_xml_one_vector(out,name,j);
          out << end_tag("EIGENSTATE");   
        }
        out << end_tag("EIGENSTATES");
      }
    }
  }
}
   
template <class T, class G>
void DiagTask<T,G>::handle_tag(std::istream& infile, const XMLTag& intag) 
{
  XMLTag tag(intag);
  
  // we don't need to read the XML file if the HDF-5 file has already been read
  if (this->read_hdf5_) {
     skip_element(infile,tag);
    return;
  }
  
  if (intag.type==XMLTag::SINGLE)
    return;
  if (intag.name=="EIGENVALUES") {
    // we don't need to read the XML file if the HDF-5 file has already been read

    std::vector<std::pair<std::string,std::string> > qnvals;
    std::vector<magnitude_type> evals;
    char c;
    infile >> c;
    while (c=='<' && infile) {
      infile.putback(c);
      tag=parse_tag(infile);
      if (tag.name=="QUANTUMNUMBER") {
        qnvals.push_back(std::make_pair(tag.attributes["name"],tag.attributes["value"]));
      }
      else if (tag.name=="/EIGENVALUES")
        return;
      skip_element(infile,tag);
      infile >> c;
    }
    do {
      infile.putback(c);
      magnitude_type ev;
      infile >> ev >> c;
      evals.push_back(ev);
    } while (c!='<' && infile);
    infile.putback(c);
    tag=parse_tag(infile);
    if (tag.name!="/EIGENVALUES")
      boost::throw_exception(std::runtime_error("Encountered unexpected tag " + tag.name 
                                  + " while pasrsing <EIGENVALUES>\n"));
    mag_vector_type evals_vector(evals.size());
    std::copy(evals.begin(),evals.end(),evals_vector.begin());
    eigenvalues_.push_back(evals_vector);
    this->quantumnumbervalues_.push_back(qnvals);
  }
  else if (intag.name=="EIGENSTATES") {
    measurements_.push_back(EigenvectorMeasurements<value_type>(*this));
    std::vector<std::pair<std::string,half_integer_type> > qnvals;
    XMLTag tag=parse_tag(infile);
    while (tag.name !="/EIGENSTATES") {
      if (tag.name=="QUANTUMNUMBER")
        skip_element(infile,tag);
      else if (tag.name=="EIGENSTATE") {
        if (tag.type != XMLTag::SINGLE) {
          tag = parse_tag(infile);
          tag = measurements_.rbegin()->handle_tag(infile,tag);
          if (tag.name!="/EIGENSTATE") 
            boost::throw_exception(std::runtime_error("unexpected element " + tag.name + " inside <EIGENSTATE>"));
        }
      }
      tag=parse_tag(infile);
    }
  }    
  else
    skip_element(infile,intag);
    
}

} } // name space scheduler

#endif // ALPS_SCHEDULER_DIAG_HPP

