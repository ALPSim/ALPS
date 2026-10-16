/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>
* Copyright (C) 2026 by the ALPS collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: dmrg.h 2552 2007-09-10 18:53:03Z afeiguin $ */

// The dmrg application's task: parameters, measurement definitions and
// output as before, solved by ALPS::classic_dmrg, which runs Adrian
// Feiguin's dmtk engine (formerly copied into this directory) on numeric
// operators. The XML catalog Hamiltonian and the measured operators are
// compiled with alps::models::CatalogModel (<alps/model/catalog.hpp>), the
// evaluator ALPS::ed is checked against, instead of being translated term by
// term into dmtk operators here. Parameters, measurement names and the output
// layout are those of the earlier task in this file; the routing itself was
// written in ALPS for this change and is not taken from alpsx.

#include <alps/dmrg/classic.hpp>
#include <alps/model/catalog.hpp>

#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/scheduler.h>
#include <alps/scheduler/measurement_operators.h>
#include <alps/numeric/real.hpp>
#include <alps/type_traits/is_complex.hpp>
#include <alps/utility/os.hpp>
#include <alps/hdf5.hpp>

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#define DMRG_VERSION "1.0.0"
#define DMRG_DATE "2006/10/02"

inline void print_dmrg_copyright(std::ostream& os)
{
  os << "ALPS/dmrg version " DMRG_VERSION " (" DMRG_DATE ")\n"
     << "  Density Matrix Renormalization Group algorithm\n"
     << "  for low-dimensional interacting systems.\n"
     << "  available from http://alps.comp-phys.org/\n"
     << "  copyright (c) 2006-2013 by Adrian E. Feiguin\n\n";
}

template<class value_type>
class DMRGTask
 : public alps::scheduler::Task
 , public alps::graph_helper<>
 , public alps::model_helper<>
 , protected alps::EigenvectorMeasurements<value_type >
{
public:
  DMRGTask<value_type>(const alps::ProcessList& , const boost::filesystem::path& );
  DMRGTask<value_type>(const alps::ProcessList& w, const alps::Parameters& p);
  void dostep();
  void write_xml_body(alps::oxstream&, const boost::filesystem::path&,bool) const;

  static void print_copyright(std::ostream& os = std::cout)
  {
    print_dmrg_copyright(os);
  }

  void save(alps::hdf5::archive &) const;

  // iteration measurements
  std::map<std::string,std::vector<double> > iteration_measurements;

private:
  typedef alps::models::CatalogModel<value_type> catalog_type;

  alps::SiteOperator make_site_term(std::string x)
  {
    if (x[x.size()-1]!=')')
      x += "(i)";
    alps::SiteOperator op(x,"i");
    substitute_operators(op,parms);
    return op;
  }

  void init();
  alps::dmrg::classic::Measurements requests(catalog_type const& model);
  void store(alps::dmrg::classic::Result const& result);

  int num_eigenvalues;
  // Conserved charges: those of CONSERVED_QUANTUMNUMBERS with a <name>_total.
  std::vector<std::string> conserved;
  alps::dmrg::classic::Options options;
};

// Measurement request names: the expression kind, its name and the location.
inline std::string dmrg_request(char kind, std::string const& name, std::size_t first = 0, std::size_t second = 0)
{
  return std::string(1, kind) + '\x1f' + name + '\x1f' + std::to_string(first) + '\x1f' + std::to_string(second);
}

template<class value_type>
DMRGTask<value_type>::DMRGTask(const alps::ProcessList& w,const boost::filesystem::path& fn)
  : alps::scheduler::Task(w,fn)
  , alps::graph_helper<>(parms)
  , alps::model_helper<>(*this,parms)
  , alps::EigenvectorMeasurements<value_type >(*this)
{
  init();
}

template<class value_type>
DMRGTask<value_type>::DMRGTask(const alps::ProcessList& w,const alps::Parameters& p)
  : alps::scheduler::Task(w,p)
  , alps::graph_helper<>(parms)
  , alps::model_helper<>(*this,parms)
  , alps::EigenvectorMeasurements<value_type >(*this)
{
  init();
}


template<class value_type>
void DMRGTask<value_type>::init()
{
  namespace classic = alps::dmrg::classic;
  num_eigenvalues = this->parms.value_or_default("NUMBER_EIGENVALUES",1);

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" ,");

  // read number of sweeps
  int const num_sweeps = parms.value_or_default("SWEEPS",4);
  // Resuming from the scratch files of an interrupted run is gone with the
  // application's private engine: the solver owns its scratch directory.
  if (static_cast<int>(parms.value_or_default("START_SWEEP",0)) != 0)
    boost::throw_exception(std::runtime_error("START_SWEEP (with START_DIR and START_ITER) is no longer supported: "
                                              "the dmrg application cannot resume an interrupted run"));
  int const verbose = parms.value_or_default("VERBOSE",0);

  // read number of states
  int nwarmup = 20;
  if (parms.defined("NUM_WARMUP_STATES")) {
    nwarmup = static_cast<int>(parms["NUM_WARMUP_STATES"]);
  }
  double error = -1.;
  if (parms.defined("TRUNCATION_ERROR")) {
    error = static_cast<double>(parms["TRUNCATION_ERROR"]);
  }
  double lanczos_tol = -1.;
  if (parms.defined("LANCZOS_TOLERANCE")) {
    lanczos_tol = static_cast<double>(parms["LANCZOS_TOLERANCE"]);
  }
  int maxstates = -1;
  std::vector<int> num_states;
  if (parms.defined("STATES")) {
    std::string states_string = parms["STATES"];
    tokenizer state_tokens(states_string, sep);
    for (tokenizer::const_iterator it = state_tokens.begin(); it !=state_tokens.end();++it)
      num_states.push_back(boost::lexical_cast<int>(*it));
    if (num_states.size() < 2*num_sweeps)
      boost::throw_exception(std::runtime_error("Need to specify either 2*SWEEPS different values in STATES, or one MAXSTATES value"));
  }
  else if (parms.defined("MAXSTATES")) {
    maxstates = static_cast<int>(parms["MAXSTATES"]);
    for (int i=0; i< 2*num_sweeps; ++i)
      num_states.push_back((i+1)*maxstates/(2*num_sweeps));
  }
  else if (parms.defined("NUMSTATES")) {
    maxstates = static_cast<int>(parms["NUMSTATES"]);
    for (int i=0; i< 2*num_sweeps; ++i)
      num_states.push_back(maxstates);
  }
  else
    boost::throw_exception(std::runtime_error("Need to specify either 2*SWEEPS different values in STATES, or one MAXSTATES value"));
  if (num_sweeps < 1)
    boost::throw_exception(std::runtime_error("SWEEPS must be positive"));
  for (int i = 0; i < 2*num_sweeps; ++i)
    if (num_states[i] < 1)
      boost::throw_exception(std::runtime_error("Every half sweep must keep at least one state; increase MAXSTATES or STATES"));

  // the quantum numbers of every site type
  std::set<std::string> qns;
  for (site_iterator it = sites().first ; it != sites().second ; ++it) {
    std::set<std::string> newqns = quantum_numbers(site_type(*it));
    qns.insert(newqns.begin(),newqns.end());
  }

  // Only the listed quantum numbers given a <name>_total are conserved; the
  // other names of CONSERVED_QUANTUMNUMBERS are ignored, as they always were
  // ("N,Sz" for the spin model).
  if (parms.defined("CONSERVED_QUANTUMNUMBERS")) {
    std::string qn_string = parms["CONSERVED_QUANTUMNUMBERS"];
    tokenizer qn_tokens(qn_string, sep);
    for (tokenizer::const_iterator it = qn_tokens.begin(); it != qn_tokens.end(); ++it) {
      std::string const name = *it;
      if (!parms.defined(name+"_total") || options.sector.count(name))
        continue;
      if (!qns.count(name))
        boost::throw_exception(std::runtime_error("Quantum number " + name + " is not defined in the model" ));
      options.sector[name] = alps::evaluate<double>(static_cast<std::string>(parms[name+"_total"]),parms);
      conserved.push_back(name);
    }
  }

  options.sweeps.clear();
  for (int i = 0; i < num_sweeps; ++i)
    options.sweeps.push_back(classic::SweepStates{static_cast<std::size_t>(num_states[2*i]),
                                                  static_cast<std::size_t>(num_states[2*i+1])});
  options.warmup_states = nwarmup > 0 ? static_cast<std::size_t>(nwarmup) : 20;
  options.num_eigenvalues = num_eigenvalues > 0 ? static_cast<std::size_t>(num_eigenvalues) : 1;
  if (error > 0.) {
    options.truncation_error = error;
    if (maxstates > 0)
      options.truncation_max_states = static_cast<std::size_t>(maxstates);
  }
  if (lanczos_tol > 0.)
    options.lanczos_tolerance = lanczos_tol;
  options.arithmetic = alps::is_complex<value_type>::value ? classic::Arithmetic::complex : classic::Arithmetic::real;
  options.verbose = verbose > 0 ? verbose : 0;
  options.log = &std::cout;
}

// One expectation value per site or bond (local), per expression (average)
// and per ordered site pair (correlations), compiled by the catalog model.
template<class value_type>
alps::dmrg::classic::Measurements
DMRGTask<value_type>::requests(catalog_type const& model)
{
  typedef std::pair<std::string,std::string> string_pair;
  typedef std::pair<std::string,std::pair<std::string,std::string> > string_string_pair_pair;
  typedef typename catalog_type::site_descriptor catalog_site;
  alps::dmrg::classic::Measurements result;
  auto& expectation = result.expectation;

  BOOST_FOREACH (string_pair const& ex, this->local_expressions) {
    if (has_bond_operator(ex.second)) {
      alps::BondOperator const op = this->get_bond_operator(ex.second,parms);
      std::size_t i = 0;
      for (auto it = model.graph().bonds().first; it != model.graph().bonds().second; ++it, ++i)
        expectation[dmrg_request('L', ex.first, i)] = model.compile_operator(op, *it);
    }
    else {
      alps::SiteOperator const op = make_site_term(ex.second);
      for (auto it = model.graph().sites().first; it != model.graph().sites().second; ++it)
        expectation[dmrg_request('L', ex.first, *it)] = model.compile_operator(op, *it);
    }
  }

  BOOST_FOREACH (string_pair const& ex, this->average_expressions) {
    if (has_bond_operator(ex.second))
      expectation[dmrg_request('A', ex.first)] = model.compile_operator(this->get_bond_operator(ex.second,parms));
    else
      expectation[dmrg_request('A', ex.first)] = model.compile_operator(make_site_term(ex.second));
  }

  BOOST_FOREACH (string_string_pair_pair const& ex, this->correlation_expressions) {
    alps::SiteOperator const onsite = make_site_term(ex.second.first+"(i)*"+ex.second.second+"(i)");
    alps::BondOperator pair(ex.second.first+"(i)*"+ex.second.second+"(j)","i","j");
    substitute_operators(pair,parms);
    // The pair tensor depends on the two site types only; evaluate it once per type pair.
    std::map<std::pair<int,int>, alps::multi_array<std::pair<value_type,bool>,4> > tensors;
    for (site_iterator s1 = sites().first; s1 != sites().second; ++s1)
      for (site_iterator s2 = sites().first; s2 != sites().second; ++s2) {
        catalog_site const first = static_cast<catalog_site>(*s1), second = static_cast<catalog_site>(*s2);
        std::string const name = dmrg_request('C', ex.first, *s1, *s2);
        if (*s1 == *s2) {
          expectation[name] = model.compile_operator(onsite, first);
          continue;
        }
        std::pair<int,int> const types(site_type(*s1), site_type(*s2));
        auto found = tensors.find(types);
        if (found == tensors.end())
          found = tensors.emplace(types, model.local_matrix(pair, first, second)).first;
        expectation[name] = model.compile_operator(found->second, first, second);
      }
  }
  return result;
}

template<class value_type>
void DMRGTask<value_type>::dostep()
{
  if (finished())
    return;

  // Select scratch storage when this task runs; other tasks may have been
  // constructed since init(). The solver works in a private directory below
  // it and removes that directory when it returns.
  if (parms.defined("TEMP_DIRECTORY"))
    options.scratch_parent = static_cast<std::string>(parms["TEMP_DIRECTORY"]);
  else
    options.scratch_parent = alps::temp_directory_path().string();

  // The scheduler seeded the disorder generator before this task built its
  // lattice; seed it again so the catalog model builds the same graph.
  alps::Disorder::seed(parms.value_or_default("DISORDERSEED",0));
  catalog_type const model(parms);
  if (model.graph().num_sites() != num_sites())
    boost::throw_exception(std::logic_error("dmrg: the catalog model and the task disagree on the lattice"));
  alps::models::CompiledOperator const hamiltonian = model.compile_hamiltonian();
  alps::models::PhysicalBasis const basis = model.physical_basis(conserved);

  alps::dmrg::classic::Result const result =
    alps::dmrg::classic::solve(hamiltonian, basis, options, requests(model));
  store(result);
  finish();
}


template<class value_type>
void
DMRGTask<value_type>::store(alps::dmrg::classic::Result const& result)
{
  typedef std::pair<std::string,std::string> string_pair;
  typedef std::pair<std::string,std::pair<std::string,std::string> > string_string_pair_pair;
  auto value = [&](std::string const& name, std::size_t target) {
    return result.measurements.at(name).values.at(target).at(0);
  };
  std::size_t const targets = result.energies.size();

  for (auto const& iteration : result.iterations) {
    // dmtk's LEFT2RIGHT is 0 and RIGHT2LEFT is 1.
    iteration_measurements["Direction"].push_back(
      iteration.direction == alps::dmrg::classic::Direction::right_to_left ? 1. : 0.);
    iteration_measurements["Iteration"].push_back(static_cast<double>(iteration.index));
    iteration_measurements["Energy"].push_back(iteration.energy);
    iteration_measurements["Truncation Error"].push_back(iteration.truncation_error);
    iteration_measurements["Entropy"].push_back(iteration.entropy);
  }

  std::vector<unsigned int> const multiplicity = distance_multiplicities();
  for (std::size_t target = 0; target < targets; ++target) {
    BOOST_FOREACH (string_pair const& ex, this->local_expressions) {
      std::vector<value_type> values;
      std::size_t const count = has_bond_operator(ex.second) ? num_bonds() : num_sites();
      for (std::size_t i = 0; i < count; ++i)
        values.push_back(value(dmrg_request('L', ex.first, i), target));
      this->local_values[ex.first].push_back(values);
    }
    BOOST_FOREACH (string_pair const& ex, this->average_expressions)
      this->average_values[ex.first].push_back(value(dmrg_request('A', ex.first), target));
    // The average over the site pairs at each distance.
    BOOST_FOREACH (string_string_pair_pair const& ex, this->correlation_expressions) {
      std::vector<value_type> values(num_distances(), value_type(0.));
      for (site_iterator s1 = sites().first; s1 != sites().second; ++s1)
        for (site_iterator s2 = sites().first; s2 != sites().second; ++s2) {
          std::size_t const d = distance(*s1,*s2);
          values[d] += value(dmrg_request('C', ex.first, *s1, *s2), target) / static_cast<double>(multiplicity[d]);
        }
      this->correlation_values[ex.first].push_back(values);
    }
  }
  for (std::size_t i = 0; i < targets; ++i)
    this->average_values["Energy"].push_back(result.energies[i]);
  this->average_values["Truncation error"].push_back(result.final_truncation_error);
}

#ifdef ALPS_HAVE_HDF5
template<class value_type>
void DMRGTask<value_type>::save(alps::hdf5::archive & ar) const
{
  using alps::numeric::real;
  alps::scheduler::Task::save(ar);
  typename std::map<std::string,std::vector<value_type> >::const_iterator it = this->average_values.find("Energy");
  if (it != this->average_values.end()) {
    std::vector<double> energies = real(it->second);
    ar["spectrum/energies"] << energies;
  }

  std::string context = ar.get_context();
  ar.set_context(ar.complete_path("spectrum"));
  alps::EigenvectorMeasurements<value_type>::save(ar);
  ar.set_context(context);

  typedef typename std::map<std::string,std::vector<double> >::const_iterator IT;
  for (IT it=iteration_measurements.begin(); it != iteration_measurements.end();++it)
      ar["simulation/results/Iteration "+alps::hdf5_name_encode(it->first)+"/mean/value"] << it->second;
}
#endif

template<class value_type>
void DMRGTask<value_type>::write_xml_body(alps::oxstream& out, const boost::filesystem::path& p, bool writeallxml) const
{
  if (writeallxml) {
    out << alps::start_tag("EIGENSTATES") << alps::attribute("number",num_eigenvalues);
    for (int j=0;j<num_eigenvalues;++j) {
      out << alps::start_tag("EIGENSTATE") << alps::attribute("number",j);
      this->write_xml_one_vector(out,p,j);
      out << alps::end_tag("EIGENSTATE");
    }
    out << alps::end_tag("EIGENSTATES");
  }
}
