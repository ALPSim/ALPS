/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2009 by Matthias Troyer <troyer@comp-phys.org>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <cassert>

#include <alps/config.h> // needed to set up correct bindings
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/lapack/driver/heev.hpp>

#include "../diag.h"
#include "measurementplots.h"

#include <alps/plot.h>

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/ublas_sparse_functions.hpp>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/upper.hpp>
#include <boost/numeric/bindings/lower.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/trim.hpp>

template <class T> 
struct pick_real_complex 
{
  template <class Matrix, class Vector>
  static int heev (const char jobz, Matrix const& m, Vector& v)
   { return boost::numeric::bindings::lapack::syev(jobz,m,v);}
};

template <class T>
struct pick_real_complex<std::complex<T> > 
{
  template <class Matrix, class Vector>
  static int heev (const char jobz, Matrix const& m, Vector& v)
   { return boost::numeric::bindings::lapack::heev(jobz,m,v);}
};

template <class T>
class FullDiagMatrix : public DiagMatrix<T,alps::numeric::matrix<T> >
{
public:
  typedef DiagMatrix<T,alps::numeric::matrix<T> > super_type;
  typedef T value_type;
  typedef typename super_type::magnitude_type magnitude_type;
  typedef typename super_type::matrix_type matrix_type;
  typedef typename super_type::vector_type vector_type;
  typedef typename super_type::size_type size_type;
  typedef typename super_type::mag_vector_type mag_vector_type;
  typedef typename super_type::half_integer_type half_integer_type;
  typedef typename super_type::site_iterator site_iterator;
  typedef typename super_type::site_descriptor site_descriptor;
  typedef typename super_type::operator_matrix_type operator_matrix_type;
  
  // settings, if given, are the parameters already read by read_parameters
  FullDiagMatrix (const alps::ProcessList& where , const boost::filesystem::path& p,
                  const alps::params* settings = 0);
  void evaluate(const alps::Parameters&, const std::string&) const;
  void print_eigenvectors(std::ostream& os) const;
  void append_eigenvectors(std::vector<value_type>& v) const;
private:
  // a temperature or field, which can be an expression in the model parameters
  /** @throws std::invalid_argument naming the parameter if it is empty or cannot be evaluated */
  static double parameter_value(const alps::Parameters& p, const std::string& name);
  magnitude_type groundstate_energy() const;
  void do_subspace();
  magnitude_type calculate_averages(MeasurementsPlot<magnitude_type>&) const;
  void write_xml_body(alps::oxstream&, const boost::filesystem::path&,bool) const;
  std::vector<value_type> calculate(operator_matrix_type const& m) const;

  mutable magnitude_type energy;
  mutable magnitude_type free_energy;
  mutable magnitude_type specific_heat;
  mutable magnitude_type entropy;
  mutable magnitude_type magnetization;
  mutable magnitude_type susceptibility;
  mutable magnitude_type beta;
  mutable magnitude_type field;
  mutable magnitude_type field0;
  mutable magnitude_type min_val;
  mutable bool have_Conserved;
  mutable bool measure_properties;

  // a few definitions to allow changing the coupling to the conserved quantity
  mutable std::string conserved_name;
};

template <class T>
FullDiagMatrix<T>::FullDiagMatrix (const alps::ProcessList& where , const boost::filesystem::path& p,
                                   const alps::params* settings)
 : super_type(where,p,true,&FullDiagMatrix::define_parameters,settings)
 , field(0.)
 , field0(0.) 
 , have_Conserved(false) 
{ 
  this->construct(); 
}

    
template <class T>
double FullDiagMatrix<T>::parameter_value(const alps::Parameters& p, const std::string& name)
{
  std::string const text = boost::algorithm::trim_copy(std::string(p[name]));
  if (text.empty())
    boost::throw_exception(std::invalid_argument("parameter '" + name + "' is empty"));
  try {
    return alps::evaluate<double>(text, p);
  }
  catch (std::exception& e) {
    boost::throw_exception(std::invalid_argument("Cannot evaluate parameter '" + name + "': " + e.what()));
  }
  return 0.;
}

template <class T>
void FullDiagMatrix<T>::evaluate(const alps::Parameters& parmsin,const std::string& basename) const
{
  alps::Parameters p(this->get_parameters());
  alps::graph_helper<> lattice(this->get_parameters());
  p << parmsin;
  alps::params settings = alps::make_parameters_from_deprecated(p);
  
  bool have_t_range=false;
  bool muN = false;        // coupling -mu N ?
  magnitude_type tmin;
  magnitude_type tmax;
  magnitude_type deltaT;
  
  conserved_name = "Sz";    // Default: Coupling -h Sz
  std::string field_name = "h";
  std::string field_label = "Magnetic Field";
  std::string conserved_label = "Magnetization";
  std::string second_moment_name = "Uniform Susceptibility";

  // If you wish to change the coupling to a conserved quantity, you have to change the following block:  
  settings.define<std::string>("couple", "", "mu or MU: couple -mu N instead of -h Sz")
          .define<std::string>("versus", "", "the field name, in either case: plot versus the field instead of T");
  std::string const couple = settings["couple"].as<std::string>();
  muN = (couple == "mu");
  if(couple == "MU")
     muN = true;
  if(muN)            // Switch to coupling -mu N
   {
    conserved_name = "N";
    field_name = "mu";
    field_label = "Chemical Potential";
    conserved_label = "Particle number";
    second_moment_name = "Compressibility";
   }
  // ... until here

  std::string field_name_uc = boost::algorithm::to_upper_copy(field_name);
  std::string conserved_label_lc = boost::algorithm::to_lower_copy(conserved_label);
  boost::algorithm::replace_all(conserved_label_lc, " ", "_");
  std::string second_moment_name_lc = boost::algorithm::to_lower_copy(second_moment_name);
  boost::algorithm::replace_all(second_moment_name_lc, " ", "_");

  // define only what this evaluation reads; the others need not be numbers.
  // The temperatures and fields can be expressions in the model parameters,
  // such as h = J/2, so they are evaluated below instead
  std::string const measure_name = muN?"MEASURE_CHARGE_PROPERTIES":"MEASURE_MAGNETIC_PROPERTIES";
  settings.define<bool>(measure_name, true, "measure the conserved quantity and its second moment")
          .define<bool>("DENSITIES", true, "divide by the number of sites");
  bool const have_t_range_parameters = p.defined("T_MAX") && p.defined("DELTA_T");
  bool const have_field_range = p.defined(field_name_uc+"_MAX") && p.defined("DELTA_"+field_name_uc);
  std::ostringstream errors;
  if (settings.has_missing(errors))
    boost::throw_exception(std::invalid_argument(errors.str()));

  measure_properties = settings[measure_name].as<bool>();
  bool densities = settings["DENSITIES"].as<bool>();

  std::string dstr = densities?" Density":"";
  std::string sstr = densities?" per Site":""; 

  magnitude_type perSite = densities?num_sites(lattice.graph()):1;

  std::string const versus = settings["versus"].as<std::string>();
  bool versus_field = (versus == field_name);
  if(versus == field_name_uc)
      versus_field = true;

  std::string vstr = versus_field?field_label:"Temperature";

  if (have_t_range_parameters) {
    have_t_range=true;
    tmin=p.defined("T_MIN") ? parameter_value(p, "T_MIN") : 0.;
    tmax=parameter_value(p, "T_MAX");
    deltaT=parameter_value(p, "DELTA_T");
  }
  
  if (p.defined("T"))
    beta = 1./parameter_value(p, "T");
  else if (p.defined("beta"))
    beta = parameter_value(p, "beta");
  else if (!have_t_range)
  {
    std::cerr << "DELTA_T and T_MAX or T need to be defined\n";
    return;
  }

  if (!have_t_range) {
   // define a T-range with 1 point
   tmin = tmax = 1/beta;
   deltaT = 1;
  }

  magnitude_type hmin;
  magnitude_type hmax;
  magnitude_type deltaH;
  
  if (have_field_range) {
    hmin=p.defined(field_name_uc+"_MIN") ? parameter_value(p, field_name_uc+"_MIN") : 0.;
    hmax=parameter_value(p, field_name_uc+"_MAX");
    deltaH=parameter_value(p, "DELTA_"+field_name_uc);
    field0 = p.defined(field_name) ? parameter_value(p, field_name) : 0.;
  }
  else {
   // define an h-range with 1 point
   field0 = hmin = hmax = p.defined(field_name) ? parameter_value(p, field_name) : 0.;
   deltaH = 1;
  }

  MeasurementsPlot<magnitude_type> plots (this->measurements_[0],vstr);

  alps::plot::Plot<magnitude_type> energy_plot("Energy"+dstr+" versus "+vstr);
  alps::plot::Plot<magnitude_type> free_energy_plot("Free Energy"+dstr+" versus "+vstr);
  alps::plot::Plot<magnitude_type> entropy_plot("Entropy"+dstr+" versus "+vstr);
  alps::plot::Plot<magnitude_type> specific_heat_plot("Specific Heat"+sstr+" versus "+vstr);
  alps::plot::Plot<magnitude_type> magnetization_plot(conserved_label+sstr+" versus "+vstr);
  alps::plot::Plot<magnitude_type> susceptibility_plot(second_moment_name+sstr+" versus "+vstr);

  // versus T
  if (! versus_field) {
    energy_plot.set_labels("T","Energy"+dstr);                                // add labels
    free_energy_plot.set_labels("T","Free Energy"+dstr);
    entropy_plot.set_labels("T","Entropy"+dstr);
    specific_heat_plot.set_labels("T","Specific Heat"+sstr);
    magnetization_plot.set_labels("T", conserved_label+sstr);
    susceptibility_plot.set_labels("T", second_moment_name+sstr);
    // we have to hide the legend in order to avoid trouble with the
    // plot2text tool - this isn't really nice and can hopefully be
    // removed in future releases
    energy_plot.show_legend(false);
    free_energy_plot.show_legend(false);
    entropy_plot.show_legend(false);
    specific_heat_plot.show_legend(false);
    magnetization_plot.show_legend(false);
    susceptibility_plot.show_legend(false);
    //
    for (field=hmin;field<hmax+deltaH/2.;field+=deltaH) {
        alps::plot::Set<magnitude_type> energies, free_energies, entropyP, specific_heatP,
                                        magnetizationP, susceptibilityP;
        energies << field_name+"="+boost::lexical_cast<std::string>(field);        // add legend
        free_energies << field_name+"="+boost::lexical_cast<std::string>(field);
        entropyP << field_name+"="+boost::lexical_cast<std::string>(field);
        specific_heatP << field_name+"="+boost::lexical_cast<std::string>(field);
        magnetizationP << field_name+"="+boost::lexical_cast<std::string>(field);
        susceptibilityP << field_name+"="+boost::lexical_cast<std::string>(field);

        // get ground state energy once
        min_val = groundstate_energy();
        for (magnitude_type t=tmin;t<tmax+deltaT/2.;t+=deltaT) {
          beta=1./t;
          plots.reset(t);
          double z = calculate_averages(plots);
          energies << boost::make_tuple(t, energy/perSite);
          free_energies << boost::make_tuple(t, free_energy/perSite);
          entropyP << boost::make_tuple(t, entropy/perSite);
          specific_heatP << boost::make_tuple(t, specific_heat/perSite);
          if (have_Conserved && measure_properties) {
            magnetizationP << boost::make_tuple(t, magnetization/perSite);
            susceptibilityP << boost::make_tuple(t, susceptibility/perSite);
          }
        }
        energy_plot << energies;
        free_energy_plot << free_energies;
        entropy_plot << entropyP;
        specific_heat_plot << specific_heatP;
        if (have_Conserved && measure_properties) {
          magnetization_plot << magnetizationP;
          susceptibility_plot << susceptibilityP;
        }
    }
  }
  // versus h
  else {
    energy_plot.set_labels(field_name, "Energy"+dstr);                                // add labels
    free_energy_plot.set_labels(field_name, "Free Energy"+dstr);
    entropy_plot.set_labels(field_name, "Entropy"+dstr);
    specific_heat_plot.set_labels(field_name, "Specific Heat"+sstr);
    magnetization_plot.set_labels(field_name, conserved_label+sstr);
    susceptibility_plot.set_labels(field_name, second_moment_name+sstr);
    // we have to hide the legend in order to avoid trouble with the
    // plot2text tool - this isn't really nice and can hopefully be 
    // removed in future releases
    energy_plot.show_legend(false);
    free_energy_plot.show_legend(false);
    entropy_plot.show_legend(false);
    specific_heat_plot.show_legend(false);
    magnetization_plot.show_legend(false);
    susceptibility_plot.show_legend(false);
    //
    for (magnitude_type t=tmin;t<tmax+deltaT/2.;t+=deltaT) {
        beta=1./t;
        alps::plot::Set<magnitude_type> energies, free_energies, entropyP, specific_heatP,
                                        magnetizationP, susceptibilityP;
        energies << "T="+boost::lexical_cast<std::string>(t);                // add legend
        free_energies << "T="+boost::lexical_cast<std::string>(t);
        entropyP << "T="+boost::lexical_cast<std::string>(t);
        specific_heatP << "T="+boost::lexical_cast<std::string>(t);
        magnetizationP << "T="+boost::lexical_cast<std::string>(t);
        susceptibilityP << "T="+boost::lexical_cast<std::string>(t);
        for (field=hmin;field<hmax+deltaH/2.;field+=deltaH) {
          min_val = groundstate_energy();
          plots.reset(field); // start new point
          calculate_averages(plots);
          energies << boost::make_tuple(field, energy/perSite);
          free_energies << boost::make_tuple(field, free_energy/perSite);
          entropyP << boost::make_tuple(field, entropy/perSite);
          specific_heatP << boost::make_tuple(field, specific_heat/perSite);
          if (have_Conserved && measure_properties) {
            magnetizationP << boost::make_tuple(field, magnetization/perSite);
            susceptibilityP << boost::make_tuple(field, susceptibility/perSite);
          }
        }
        energy_plot << energies;
        free_energy_plot << free_energies;
        entropy_plot << entropyP;
        specific_heat_plot << specific_heatP;
        if (have_Conserved && measure_properties) {
          magnetization_plot << magnetizationP;
          susceptibility_plot << susceptibilityP;
        }
    }
  }

  // write plots to files
  boost::filesystem::path enp(basename+".plot.energy.xml");
  alps::oxstream ens(enp);
  ens << energy_plot;
  boost::filesystem::path fenp(basename+".plot.free_energy.xml");
  alps::oxstream fens(fenp);
  fens << free_energy_plot;
  boost::filesystem::path entp(basename+".plot.entropy.xml");
  alps::oxstream ents(entp);
  ents << entropy_plot;
  boost::filesystem::path shp(basename+".plot.specific_heat.xml");
  alps::oxstream shs(shp);
  shs << specific_heat_plot;
  
  plots.write(basename);
  
  if (have_Conserved && measure_properties) {
     boost::filesystem::path magp(basename+".plot."+conserved_label_lc+".xml");
     alps::oxstream mags(magp);
     mags << magnetization_plot;
     boost::filesystem::path susp(basename+".plot."+second_moment_name_lc+".xml");
     alps::oxstream suss(susp);
     suss << susceptibility_plot;
    }
}


template <class T>
typename FullDiagMatrix<T>::magnitude_type FullDiagMatrix<T>::groundstate_energy() const
{
  magnitude_type min_val=std::numeric_limits<magnitude_type>::max();
  typename std::vector<mag_vector_type>::const_iterator vit;
  unsigned int i=0;
  have_Conserved=true;
  for (vit=this->eigenvalues_.begin(); vit!=this->eigenvalues_.end(); ++vit, ++i) {
    magnitude_type mval=0.;
    if (have_Conserved) {
      int j;
      for (j=0;j<this->quantumnumbervalues_[i].size() && this->quantumnumbervalues_[i][j].first!=conserved_name;++j)
        ; // search for QN
          
      if (j==this->quantumnumbervalues_[i].size())
        have_Conserved=false;
      else
        mval=alps::expression::Expression<double>(this->quantumnumbervalues_[i][j].second).value();
    }
    typename mag_vector_type::const_iterator it = std::min_element(vit->begin(), vit->end());
    if(it != vit->end() && *it - mval * (field-field0) < min_val)
      min_val = *it - mval * (field-field0);
  }
  return min_val;
}

template <class T>
typename FullDiagMatrix<T>::magnitude_type FullDiagMatrix<T>::calculate_averages(MeasurementsPlot<magnitude_type>& plots) const
{
  magnitude_type z=0.;
  magnitude_type en=0.;
  magnitude_type en2=0.;
  magnitude_type m=0.;
  magnitude_type m2=0.;
  
  

//  min_val = groundstate_energy();

  // loop over sectors
  unsigned int i=0;
  have_Conserved=true;

  
  for (typename std::vector<mag_vector_type>::const_iterator vit=this->eigenvalues_.begin();
       vit!=this->eigenvalues_.end();++vit,++i) {
    magnitude_type mval=0.;
    if (have_Conserved) {
      int j;
      for (j=0;j<this->quantumnumbervalues_[i].size() && this->quantumnumbervalues_[i][j].first != conserved_name; ++j)
        ; // search for QN
          
      if (j==this->quantumnumbervalues_[i].size())
        have_Conserved=false;
      else
        mval=alps::expression::Expression<double>(this->quantumnumbervalues_[i][j].second).value();
    }
    
    // now compute observables
    typename mag_vector_type::const_iterator it;
    int j;
    for (it=vit->begin(), j=0;it!=vit->end();++it,++j) {
      magnitude_type val = *it - mval * (field-field0);
      magnitude_type w = std::exp(-beta*(val-min_val));
      z += w;
      en += val * w;
      en2 += val * val * w;
      if (have_Conserved) {
        m += mval * w;
        m2 += mval * mval * w;
      }
      // now do custom measurements
      plots.add(w,this->measurements_[i],j);
    }
  }

  energy = en/z;
  free_energy = -std::log(z)/beta+min_val;
  specific_heat = (en2/z-(en/z)*(en/z))*beta*beta;
  entropy = beta*(energy-free_energy);
  if (have_Conserved) {
    magnetization = m/z;
    susceptibility = (m2/z-(m/z)*(m/z))*beta;
   }
   
   // now do all the custom measurements
  plots.normalize(z);
    
  return z;
}


template <class T>
void FullDiagMatrix<T>::write_xml_body(alps::oxstream& out, const boost::filesystem::path& p, bool writeallxml) const
{
    const alps::Parameters& pa = this->get_parameters();
    if (writeallxml && (pa.defined("T") || pa.defined("beta"))) {
        // T, beta, H and h can be expressions in the model parameters, such as h = J/2
        beta = (pa.defined("T") ? 1./parameter_value(pa, "T")
                           : parameter_value(pa, "beta"));
        field0 = field = (pa.defined("H") ? parameter_value(pa, "H")
                          : pa.defined("h") ? parameter_value(pa, "h") : 0.);
    conserved_name = "Sz";    // Here we support only a coupling -h Sz

    min_val = groundstate_energy();
    MeasurementsPlot<magnitude_type> noplots (alps::MeasurementLabels(*this,0),"");

    calculate_averages(noplots);
    // TODO: evaluate the quantities just for one value

    out << alps::start_tag("AVERAGES") 
        << alps::start_tag("SCALAR_AVERAGE") << alps::attribute("name","Energy") << alps::no_linebreak
        << alps::start_tag("MEAN") << alps::no_linebreak<< energy << alps::end_tag("MEAN")
        << alps::end_tag("SCALAR_AVERAGE");
    out << alps::start_tag("SCALAR_AVERAGE") << alps::attribute("name","Free energy") << alps::no_linebreak
        << alps::start_tag("MEAN") << alps::no_linebreak << free_energy << alps::end_tag("MEAN")
        << alps::end_tag("SCALAR_AVERAGE");
    out << alps::start_tag("SCALAR_AVERAGE") << alps::attribute("name","Specific Heat") << alps::no_linebreak
        << alps::start_tag("MEAN") << alps::no_linebreak << specific_heat << alps::end_tag("MEAN")
        << alps::end_tag("SCALAR_AVERAGE"); 
    out << alps::start_tag("SCALAR_AVERAGE") << alps::attribute("name","Entropy") << alps::no_linebreak
        << alps::start_tag("MEAN") << alps::no_linebreak << entropy << alps::end_tag("MEAN")
        << alps::end_tag("SCALAR_AVERAGE");

    if (have_Conserved) {
      out << alps::start_tag("SCALAR_AVERAGE") << alps::attribute("name","Magnetization") << alps::no_linebreak
          << alps::start_tag("MEAN") << alps::no_linebreak << magnetization << alps::end_tag("MEAN")
          << alps::end_tag("SCALAR_AVERAGE");
      out << alps::start_tag("SCALAR_AVERAGE") << alps::attribute("name","Susceptibility") << alps::no_linebreak
          << alps::start_tag("MEAN") << alps::no_linebreak << susceptibility << alps::end_tag("MEAN")
          << alps::end_tag("SCALAR_AVERAGE"); 
    }
    out << alps::end_tag("AVERAGES");
  }
  super_type::write_xml_body(out,p,writeallxml);
}
   
template <class T>
std::vector<T> FullDiagMatrix<T>::calculate(operator_matrix_type const& m) const
{
  std::vector<value_type> av;
  for (unsigned i=0; i < num_cols(this->matrix()); ++i) {
    alps::numeric::column_view<matrix_type const> const v(this->matrix(),i);
    //av.push_back(scalar_product(conj(v),m*v)); // scalar_product uses conj_mult [M. Pikulski]
    av.push_back(scalar_product(v,m*v));
  }
  return av;
}


template <class T>
void FullDiagMatrix<T>::append_eigenvectors(std::vector<value_type>& v) const
{
  for (unsigned i=0;i<num_cols(this->matrix());++i)
    for (unsigned j=0;j<num_rows(this->matrix());++j)
      v.push_back(this->matrix()(j,i));
}

template <class T>
void FullDiagMatrix<T>::print_eigenvectors(std::ostream& os) const
{
  for (unsigned i=0;i<num_cols(this->matrix());++i)
    os << "Eigenvector# " << i << ":\n" << alps::numeric::column_view<matrix_type const>(this->matrix(),i) << "\n";
}

template <class T>
void FullDiagMatrix<T>::do_subspace()
{
  using std::copy;
  this->build();
  if (this->dimension()) {
    mag_vector_type eigenvalues(this->dimension());
    pick_real_complex<T>::heev(this->calc_vectors() ? 'V' : 'N',boost::numeric::bindings::upper(this->matrix()),eigenvalues);
    this->perform_measurements();
    copy(eigenvalues.begin(),eigenvalues.end(),std::back_inserter(this->measurements_.rbegin()->average_values["Energy"]));
    this->eigenvalues_.push_back(eigenvalues);
  }
}
