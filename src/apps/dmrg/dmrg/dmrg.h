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

/* $Id: dmrg.h 2552 2007-09-10 18:53:03Z afeiguin $ */

#define WITH_LAPACK
//#define WITH_WARNINGS

#include "dmtk/dmtk.h"
#include "mps_export.h"
#include "checkpoint_cache.h"

#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/scheduler.h>
#include <alps/scheduler/measurement_operators.h>
#include <alps/numeric/real.hpp>
#include <alps/utilities/os.hpp>
#include <alps/scheduler.h>
#include <alps/params.hpp>
#include <alps/ngs/make_parameters_from_deprecated.hpp>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/archive/tmpdir.hpp>
#include <boost/foreach.hpp>

#include <optional>
#include <sstream>
#include <stdexcept>


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

#include <alps/hdf5.hpp>

template<class value_type>
class DMRGTask 
 : public alps::scheduler::Task
 , public alps::graph_helper<>
 , public alps::model_helper<>
 , protected alps::EigenvectorMeasurements<value_type >
{
public:  
  typedef alps::half_integer<short> half_integer_type;
  DMRGTask<value_type>(const alps::ProcessList& , const boost::filesystem::path& );
  DMRGTask<value_type>(const alps::ProcessList& w, const alps::Parameters& p);

  /// Defines the parameters dmrg reads, with their defaults
  static void define_parameters(alps::params& parameters);

  /// A task's parameters, with the parameters define_parameters defines parsed
  /** @throws std::invalid_argument naming each of them that does not parse */
  static alps::params read_parameters(const alps::Parameters& parms);

  void dostep();
  void write_xml_body(alps::oxstream&, const boost::filesystem::path&,bool) const;

  static void print_copyright(std::ostream& os = std::cout) 
  {
    print_dmrg_copyright(os);
  }

  void save(alps::hdf5::archive &) const;

  // iteration measurements
  std::map<std::string,std::vector<double> > iteration_measurements;

  void checkpoint_iteration(const dmtk::System<value_type >& S);

private:

  alps::SiteOperator make_site_term(std::string x)
  {
    if (x[x.size()-1]!=')')
      x += "(i)";
    alps::SiteOperator op(x,"i");
    substitute_operators(op,parms);
    return op;
  }
  
  
  void init();
  void use_temporary_directory();
  dmtk::BasicOp<value_type > create_site_operator(std::string const& name, alps::SiteOperator const& siteop, int type);
  void build_site_operator(alps::SiteOperator const& siteop, int site, dmtk::Hami<value_type > &);
  void build_bond_operator(alps::BondOperator const& bondop, bond_descriptor const& b, dmtk::Hami<value_type > &this_hami);  
  void build_2site_operator(std::pair<alps::SiteOperator,alps::SiteOperator> const& siteops, 
                            std::pair<int,int> sites, dmtk::Hami<value_type > &this_hami);
  void build_nsite_operator(std::vector<alps::nsite_product<value_type> > const& products,
                                alps::NSiteOperator const& op, alps::nsite_group const& p, dmtk::Hami<value_type > &this_hami);
  std::vector<alps::nsite_product<value_type> > split_nsite_operator(alps::NSiteOperator const& op, alps::nsite_group const& p);
 
  void save_results(); 
  void extract_mps();
  
  int num_sweeps;
  std::vector<int> num_states;
  std::vector<std::string> quantumnumber_names;
  std::vector<bool> conserved_quantumnumber;
  std::vector<half_integer_type> conserved_quantumnumber_value;
  int qnmask;

  dmtk::System<value_type > system;
  dmtk::Hami<value_type > hami;
  dmtk::Lattice lattice;
  std::vector<dmtk::Block<value_type > > site_block;
  // site operators (by site type and name) whose matrix is zero
  std::map<std::pair<int,std::string>,bool> zero_site_operator;
  bool warned_cutoff;
  
  int num_eigenvalues;
  int verbose;
  int nwarmup;
  int maxstates;
  double error;
  double lanczos_tol;

  // For resuming a previous run
  int start_sweep;
  int start_dir;
  int start_iter;

  // Failed runs keep their files and a rerun resumes from them
  dmrg_checkpoint::Cache checkpoint;
  bool in_sweeps = false;

  // SAVE_MPS: final state as tensors [left bond][local state][right bond]
  bool save_mps;
  std::vector<std::vector<value_type> > mps_tensors;
  std::vector<std::vector<std::size_t> > mps_shapes;

  std::optional<std::string> temp_directory;
  long seed;
};



template<class value_type>
bool
handler(dmtk::System<value_type>& S, size_t signal_id, void *data)
{
  DMRGTask<value_type> &task = * (DMRGTask<value_type> *)data;
  if(signal_id == dmtk::SYSTEM_SIGNAL_END_ITER){
    task.iteration_measurements["Direction"].push_back(S.get_dir());
    task.iteration_measurements["Iteration"].push_back(S.get_iter());
    task.iteration_measurements["Energy"].push_back(S.energy[0]);
    task.iteration_measurements["Truncation Error"].push_back(S.truncation_error());
    task.iteration_measurements["Entropy"].push_back(S.entropy());
    task.checkpoint_iteration(S);
  }
  return false;
}

template<class value_type>
void DMRGTask<value_type>::checkpoint_iteration(const dmtk::System<value_type >& S)
{
  // warmup and the final sweep are redone after a failure
  if (!in_sweeps || S.in_warmup() || !checkpoint.active())
    return;
  dmrg_checkpoint::Position at = {S.get_sweep(), int(S.get_dir()), int(S.get_iter()), dmtk::random_state(), seed};
  std::ostringstream note;
  note << "sweep " << at.sweep << "/" << num_sweeps << ", "
       << (at.dir == dmtk::LEFT2RIGHT ? "left-to-right" : "right-to-left") << " iteration " << at.iter;
  checkpoint.commit(dmrg_checkpoint::write_progress(at, note.str(), iteration_measurements));
}

template<class value_type>
DMRGTask<value_type>::DMRGTask(const alps::ProcessList& w,const boost::filesystem::path& fn)
  : alps::scheduler::Task(w,fn)
  , alps::graph_helper<>(parms) 
  , alps::model_helper<>(*this,parms,false,true)
  , alps::EigenvectorMeasurements<value_type >(*this)
{
  init();
}

template<class value_type>
DMRGTask<value_type>::DMRGTask(const alps::ProcessList& w,const alps::Parameters& p)
  : alps::scheduler::Task(w,p) 
  , alps::graph_helper<>(parms) 
  , alps::model_helper<>(*this,parms,false,true)
  , alps::EigenvectorMeasurements<value_type >(*this)
{
  init();
}


template<class value_type>
void DMRGTask<value_type>::define_parameters(alps::params& parameters)
{
  parameters.define<bool>("COMPLEX", false, "use complex arithmetic")
            .define<int>("NUMBER_EIGENVALUES", 1, "number of eigenstates to target")
            .define<int>("SWEEPS", 4, "number of sweeps")
            .define<int>("START_SWEEP", 0, "sweep to resume from")
            .define<int>("START_DIR", 0, "direction to resume in: 0 left to right, 1 right to left")
            .define<int>("START_ITER", 1, "iteration to resume from")
            .define<int>("VERBOSE", 0, "print progress")
            .define<bool>("SAVE_MPS", false, "save the final state as an MPS")
            .define<int>("NUM_WARMUP_STATES", 20, "states kept in the warmup")
            .define<long>("SEED", 0, "seed of the random Lanczos start vectors")
            .define<std::string>("CONSERVED_QUANTUMNUMBERS", "", "quantum numbers fixed by <name>_total");
  // without them, the number of states alone sets the truncation and the
  // Lanczos tolerance is DMTK's
  if (parameters.supplied("TRUNCATION_ERROR"))
    parameters.define<double>("TRUNCATION_ERROR", "truncation error to aim for, with at most MAXSTATES states");
  if (parameters.supplied("LANCZOS_TOLERANCE"))
    parameters.define<double>("LANCZOS_TOLERANCE", "Lanczos tolerance");
  // define only the one init() reads: STATES before MAXSTATES before NUMSTATES
  if (parameters.supplied("STATES"))
    parameters.define<std::string>("STATES", "states kept in each half sweep");
  else if (parameters.supplied("MAXSTATES"))
    parameters.define<int>("MAXSTATES", "states kept in the last half sweep, growing linearly to it");
  else if (parameters.supplied("NUMSTATES"))
    parameters.define<int>("NUMSTATES", "states kept in each half sweep");
  if (parameters.supplied("TEMP_DIRECTORY"))
    parameters.define<std::string>("TEMP_DIRECTORY", "directory for scratch files");
  // a run given TEMP_DIRECTORY keeps its files there, as before checkpoints
  parameters.define<bool>("CHECKPOINT", !parameters.supplied("TEMP_DIRECTORY"),
                          "keep the files of a failed run next to the task files to resume from");
}

template<class value_type>
alps::params DMRGTask<value_type>::read_parameters(const alps::Parameters& parms)
{
  alps::params parameters = alps::make_parameters_from_deprecated(parms);
  define_parameters(parameters);
  std::ostringstream errors;
  if (parameters.has_missing(errors))
    boost::throw_exception(std::invalid_argument(errors.str()));
  return parameters;
}

template<class value_type>
void DMRGTask<value_type>::init()
{
  warned_cutoff = false;
  alps::params const settings = read_parameters(parms);
  num_eigenvalues = settings["NUMBER_EIGENVALUES"].as<int>();
   
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" ,");

  dmtk::QN::init();
  // read number of sweeps
  num_sweeps = settings["SWEEPS"].as<int>();
  start_sweep = settings["START_SWEEP"].as<int>();
  start_dir = settings["START_DIR"].as<int>();
  if(start_dir != 0 && start_dir != 1){
    boost::throw_exception(std::runtime_error("START_DIR can assume the values 0 (left-to-right) or 1 (right-to-left)"));
  }
  start_iter = settings["START_ITER"].as<int>();
  verbose = settings["VERBOSE"].as<int>();
  save_mps = settings["SAVE_MPS"].as<bool>();
  checkpoint.configure(settings);

  // read number of states
  nwarmup = settings["NUM_WARMUP_STATES"].as<int>();
  error = settings.exists("TRUNCATION_ERROR") ? settings["TRUNCATION_ERROR"].as<double>() : -1.;
  lanczos_tol = settings.exists("LANCZOS_TOLERANCE") ? settings["LANCZOS_TOLERANCE"].as<double>() : -1.;
  maxstates = -1;
  if (settings.exists("STATES")) {
    std::string states_string = settings["STATES"].as<std::string>();
    tokenizer state_tokens(states_string, sep);
    for (tokenizer::const_iterator it = state_tokens.begin(); it !=state_tokens.end();++it)
      num_states.push_back(boost::lexical_cast<int>(*it));
    if (num_states.size() < 2*num_sweeps)
      boost::throw_exception(std::runtime_error("Need to specify either 2*SWEEPS different values in STATES, or one MAXSTATES value"));
  }
  else if (settings.exists("MAXSTATES")) {
    maxstates = settings["MAXSTATES"].as<int>();
    for (int i=0; i< 2*num_sweeps; ++i)
      num_states.push_back((i+1)*maxstates/(2*num_sweeps));
  }
  else if (settings.exists("NUMSTATES")) {
    maxstates = settings["NUMSTATES"].as<int>();
    for (int i=0; i< 2*num_sweeps; ++i)
      num_states.push_back(maxstates);
  }
  else 
    boost::throw_exception(std::runtime_error("Need to specify either 2*SWEEPS different values in STATES, or one MAXSTATES value"));
  if (settings.exists("TEMP_DIRECTORY"))
    temp_directory = settings["TEMP_DIRECTORY"].as<std::string>();
  seed = settings["SEED"].as<long>();
  
  // get all quantum numbers
  std::set<std::string> qns;
  for (site_iterator it = sites().first ; it != sites().second ; ++it) {
    std::set<std::string> newqns = quantum_numbers(site_type(*it));
    qns.insert(newqns.begin(),newqns.end());
  }

  
  // read quantum numbers and their total values
  std::copy(qns.begin(),qns.end(),std::back_inserter(quantumnumber_names));
  conserved_quantumnumber.resize(quantumnumber_names.size(),false);
  conserved_quantumnumber_value.resize(quantumnumber_names.size());
  std::string qn_string = settings["CONSERVED_QUANTUMNUMBERS"].as<std::string>();
  if (!qn_string.empty()) {
    tokenizer qn_tokens(qn_string, sep);
    std::vector<std::string> conserved_quantumnumber_names;
    std::copy(qn_tokens.begin(),qn_tokens.end(),std::back_inserter(conserved_quantumnumber_names));
    for (int i=0; i<conserved_quantumnumber_names.size(); i++) {
      if (parms.defined(conserved_quantumnumber_names[i]+"_total")) {
        int j = std::find(quantumnumber_names.begin(),quantumnumber_names.end(),conserved_quantumnumber_names[i])-quantumnumber_names.begin();
        if (j >= quantumnumber_names.size())
          boost::throw_exception(std::runtime_error("Quantum number " + conserved_quantumnumber_names[i] + " is not defined in the model" ));
        conserved_quantumnumber[j] = true;
        // <name>_total can be an expression in the other parameters
        conserved_quantumnumber_value[j] = alps::evaluate<double>(static_cast<std::string>(parms[conserved_quantumnumber_names[i]+"_total"]),parms);
      }
    }
  }
}

template <class SiteOp>
std::string simplify_name(const SiteOp &op)
{
  std::string term = op.term();
  std::string arg = "("+op.site()+")";
  boost::algorithm::replace_all(term,arg,"");
  return term;
}

template<class value_type>
void DMRGTask<value_type>::use_temporary_directory()
{
  if (temp_directory) {
    dmtk::tmp_files.set_temp_dir(temp_directory->c_str());
  } else {
    dmtk::tmp_files.set_temp_dir(alps::temp_directory_path().string().c_str());
  }
}

template<class value_type>
void DMRGTask<value_type>::dostep()
{
  if (finished())
    return;

  // Select scratch storage when this task runs; other tasks may have been
  // constructed since init().
  if (!checkpoint.enabled())
    use_temporary_directory();

  dmtk::Lattice l(num_sites(),dmtk::OBC);
  hami = dmtk::Hami<value_type >(l);
  site_block.resize(alps::maximum_vertex_type(graph())+1);

// define quantum numbers first: we need to know which are conserved before we build operators

  qnmask = 0;
  for (int type  = 0 ; type < alps::maximum_vertex_type(graph())+1 ; ++type) {
    // create quantum numbers for this site block
    /*if(conserved_quantumnumber.size() > 0)*/ {
      for (int qn = 0 ; qn < site_basis(type).size(); ++qn) {
        int idx = dmtk::QN::add_qn_index(site_basis(type)[qn].name(),site_basis(type)[qn].fermionic()); 
        qnmask |= (1 << qn);
      }
    }
  }
  dmtk::QN::set_qn_mask(qnmask);

  dmtk::QN qn;
  qnmask = 0;
  for(int i = 0; i < quantumnumber_names.size(); i++) {
    if(conserved_quantumnumber[i]){
      qn[quantumnumber_names[i]] = conserved_quantumnumber_value[i];
      qnmask |= (1 << (dmtk::QN::get_qn_index(quantumnumber_names[i])));
    }
  }

//  Iterating over sites: create site blocks
  for (site_iterator it = sites().first ; it != sites().second ; ++it) {
    hami.sites(*it) = &site_block[site_type(*it)];
    site_block[site_type(*it)].clear();
  }

  // iterate over all ste types
  for (int type  = 0 ; type < alps::maximum_vertex_type(graph())+1 ; ++type) {
    // create quantum numbers for this site block
    /*if(conserved_quantumnumber.size() > 0)*/ {
      for (int qn = 0 ; qn < site_basis(type).size(); ++qn) {
        int idx = dmtk::QN::add_qn_index(site_basis(type)[qn].name(),site_basis(type)[qn].fermionic()); 
      }
    }

    // create site basis for this block
    alps::site_basis<short> b(site_basis(type));
    dmtk::Basis basis(b.size());
    // iterate over basis states s
    for (int s=0 ; s<b.size();++s) {
      dmtk::QN real_qn;
      // extract values of the quantum numbers
      for (int qn = 0 ; qn < site_basis(type).size() ; ++qn){
        if (verbose)
          std::cout << site_basis(type)[qn].name() << "=" << b[s][qn] << "  ";
        int idx = dmtk::QN::get_qn_index(site_basis(type)[qn].name()); 
        if (idx < QN_MAX_SIZE)
          real_qn[idx] = b[s][qn];
      }
      // create the basis state for this block
      basis[s] = dmtk::State(s,real_qn);
    }
    basis.reorder();
    site_block[type].resize(basis);
    site_block[type].set_lattice(dmtk::Lattice(1,dmtk::OBC));
  }

  // create site terms
  for (site_iterator it = sites().first ; it != sites().second ; ++it)
    build_site_operator(site_term(site_type(*it)),*it, hami);

  // create bond terms
  for (bond_iterator it = bonds().first ; it != bonds().second ; ++it)
    build_bond_operator(bond_term(bond_type(*it)),*it, hami);

  // create n-site terms on the n-site groups. The split depends only on the
  // group type and the site types, so do it once per combination.
  typedef std::pair<unsigned int,std::vector<int> > nsite_group_key;
  std::map<nsite_group_key,std::vector<alps::nsite_product<value_type> > > split_cache;
  for (std::size_t q = 0; q < nsite_groups().size(); ++q) {
    alps::nsite_group const& p = nsite_groups()[q];
    alps::NSiteOperator term = model().nsite_term(p.type);
    if (term.term().empty())
      continue;
    nsite_group_key key(p.type,std::vector<int>());
    for (std::size_t r = 0; r < p.sites.size(); ++r)
      key.second.push_back(site_type(p.sites[r]));
    if (split_cache.find(key) == split_cache.end())
      split_cache[key] = split_nsite_operator(term, p);
    build_nsite_operator(split_cache[key], term, p, hami);
  }

//-------------------------------------------------------------

  if (verbose)
    std::cout << hami.description() << endl;

  // set up measurements
  
  typedef std::pair<std::string,std::string> string_pair;
 
  dmtk::Hami<value_type > meas_terms(hami);
  meas_terms.clear();

  // calculate local measurements
  BOOST_FOREACH (string_pair const& ex, this->local_expressions) {
    if (has_bond_operator(ex.second)) {
      int i=0;
      for (bond_iterator bit=bonds().first; bit!=bonds().second;++bit,++i) {
        dmtk::Hami<value_type > meas;
        build_bond_operator(get_bond_operator(ex.second,parms),*bit,meas);
        dmtk::Term<value_type > tmeas = meas[0];
//        tmeas.set_name((ex.first+"["+boost::lexical_cast<std::string>(i)+"]").c_str());
//        meas_terms += dmtk::BasicOp<value_type >(tmeas);
        meas_terms += tmeas;
      }
    }
    else {
      for (site_iterator sit=sites().first; sit!=sites().second;++sit) {
        dmtk::Hami<value_type > meas;
        build_site_operator(make_site_term(ex.second),*sit,meas);
        dmtk::Term<value_type > tmeas = meas[0];
//        tmeas.set_name((ex.first+"["+boost::lexical_cast<std::string>(*sit)+"]").c_str());
//        meas_terms += dmtk::BasicOp<value_type >(tmeas);
        meas_terms += tmeas[0];
      }
    }
  }
  
  // average measurements will be identical loops, but all terms added instead of stored separately

   BOOST_FOREACH (string_pair const& ex, this->average_expressions) {
    dmtk::Hami<value_type > meas;
    if (has_bond_operator(ex.second)) {
      for (bond_iterator bit=bonds().first; bit!=bonds().second;++bit)
        build_bond_operator(get_bond_operator(ex.second,parms),*bit,meas);
      meas.set_name(ex.first.c_str());
//meas.set_name(meas.description().c_str());
      meas_terms += dmtk::BasicOp<value_type >(meas);
    }
    else {
      for (site_iterator sit=sites().first; sit!=sites().second;++sit)
        build_site_operator(make_site_term(ex.second),*sit,meas);
      meas.set_name(ex.first.c_str());
//meas.set_name(meas.description().c_str());
//      meas_terms += dmtk::BasicOp<value_type >(meas);
      meas_terms += dmtk::BasicOp<value_type >(meas);
    }
    // store into average_values instead of local_values
    // local_values[ex.first].push_back(av);
  }
  
  // correlations

  std::vector<unsigned int> distance_mult = distance_multiplicities();
  
  // calculate correlations
  typedef std::pair<std::string,std::pair<std::string,std::string> > string_string_pair_pair;
  BOOST_FOREACH (string_string_pair_pair const& ex, this->correlation_expressions) {
    alps::SiteOperator ops1 = make_site_term(ex.second.first+"(i)");
    alps::SiteOperator ops2 = make_site_term(ex.second.second+"(i)");
    alps::SiteOperator ops = make_site_term(ex.second.first+"(i)*"+ex.second.second+"(i)");

    // use num_distances() for retrieval
    for (site_iterator sit1=this->sites().first; sit1!=this->sites().second ; ++sit1)
      for (site_iterator sit2=this->sites().first; sit2!=this->sites().second ; ++sit2) {
        std::size_t d = distance(*sit1,*sit2);
        // loop over all terms in ops1 and ops2 as above where we loop over all terms in ops
        // build terms like above
        if (*sit1 == *sit2) {
          // create matrices for combined term
          // use ops
          dmtk::Hami<value_type > meas;
          build_site_operator(ops,*sit1,meas);
          dmtk::Term<value_type > tmeas = meas[0];
          tmeas *= 1./double(distance_mult[d]);
//          tmeas.set_name((ex.first+"["+boost::lexical_cast<std::string>(d)+"]").c_str());
          meas_terms += tmeas[0];
        }
        else {
          // use ops1 and ops2
          dmtk::Hami<value_type > meas;
          build_2site_operator(std::make_pair(ops1,ops2),std::pair<int,int>(*sit1,*sit2),meas);
          dmtk::Term<value_type > tmeas = meas[0];
          tmeas *= 1./double(distance_mult[d]);
//          tmeas.set_name((ex.first+"["+boost::lexical_cast<std::string>(d)+"]").c_str());
          meas_terms += dmtk::BasicOp<value_type >(tmeas);
        }
    }
  }
  typename dmtk::Hami<value_type>::iterator iter;
  if (verbose) {
    for(iter = meas_terms.begin(); iter != meas_terms.end(); iter++)
      cout << iter->name() << " " << iter->description() << endl;

    std::cout << meas_terms.description() << "\n";
  }
  
///////////////////////////////////////////////////////////////
// Simulation
///////////////////////////////////////////////////////////////
  // the Lanczos start vectors are random: seed them from the task's SEED
  dmtk::seed_random(seed);
  hami = hami.reorder_terms();
  if(verbose)
    cout << hami.description() << endl;
  if (checkpoint.enabled()) {
    if (!checkpoint.open(dmrg_checkpoint::fingerprint(parms, hami.description(), DMRG_VERSION), infilename))
      use_temporary_directory();
    else if (temp_directory)
      std::cout << "ALPS DMRG: TEMP_DIRECTORY is not used while checkpoints are on; set CHECKPOINT=0 to use it" << std::endl;
  }

  try {
  this->system = dmtk::System<value_type >(hami,l,"ALPS");
  dmtk::System<value_type > &S = this->system;
  S.set_store_products(false);
  S.signal_connect(handler, dmtk::SYSTEM_SIGNAL_END_ITER, this);
  if(error > 0.0) S.set_error(error, maxstates);
  if(lanczos_tol > 0.0) S.set_lanczos_tolerance(lanczos_tol);
  
  S.set_calc_gap(num_eigenvalues-1); 
  dmtk::Matrix<size_t> nstates(2,num_sweeps);

  S.qnt = qn;
  S.set_qn_mask(qnmask);
  S.set_use_hc(false);
//  S.set_store_products(false);
  S.set_grow_symmetric(false);
  S.set_verbose(verbose);
  for(int i = 0; i < num_sweeps; i++)
    for(int j = 0; j < 2; j++) {
      nstates(j,i) = num_states[i*2+j];
    }
  S.start(num_sweeps, nstates);
  dmrg_checkpoint::Position at;
  std::string note;
  in_sweeps = true;
  if(checkpoint.active() && !checkpoint.record().empty()
     && dmrg_checkpoint::read_progress(checkpoint.record(), at, note, iteration_measurements)) {
    // continue after the last committed iteration
    int sweep_max = S.full_sweep() ? num_sites()-2 : num_sites()-3;
    int sweep = at.sweep, dir = at.dir, iter = at.iter + 1;
    if (iter >= sweep_max) {
      iter = 1;
      if (dir == dmtk::RIGHT2LEFT) {
        dir = dmtk::LEFT2RIGHT;
      } else {
        dir = dmtk::RIGHT2LEFT;
        ++sweep;
      }
    }
    std::cout << "ALPS DMRG: resuming from checkpoint after " << note << std::endl;
    // the results follow the SEED the checkpoint was started with, so the
    // output lists that one
    if (at.seed != seed) {
      std::cout << "ALPS DMRG: the checkpoint was started with SEED=" << at.seed
                << ", which this run continues with instead of SEED=" << seed << std::endl;
      seed = at.seed;
      parms["SEED"] = seed;
    }
    dmtk::random_state() = at.random_state;
    S.resume(sweep, dir, iter);
  } else if(start_sweep != 0) {
    S.resume(start_sweep, start_dir, start_iter);
  } else {
    S.run(nwarmup);
  }
  in_sweeps = false;
  if(S.store_products()) S.set_full_sweep(true);
  S.corr += meas_terms;
  S.final_sweep(num_states[num_states.size()-1], dmtk::RIGHT2LEFT, 1, true); 
  if(!S.store_products()) S.measure();
  if(save_mps) extract_mps();
  save_results();
  for(int i = 1; i < S._target.size(); i++){
    S.gs = S._target[i];
    S.measure();
    save_results();
  }
  for(int i = 0; i < S.energy.size(); i++) {
    this->average_values["Energy"].push_back(S.energy[i]);
  }
  this->average_values["Truncation error"].push_back(S.truncation_error());
  finish();
  } catch (...) {
    in_sweeps = false;
    checkpoint.failed();
    dmtk::tmp_files.cleanup();
    throw;
  }
  dmtk::tmp_files.cleanup();
  checkpoint.succeeded();
}


template<class value_type>
void
DMRGTask<value_type>::save_results()
{
  dmtk::System<value_type > &S = this->system;
  typename dmtk::Hami<value_type>::iterator iter = S.corr.begin();
  typedef std::pair<std::string,std::string> string_pair;
  typedef std::pair<std::string,std::pair<std::string,std::string> > string_string_pair_pair;
  using alps::numeric::real;
  
  // store local measurements
  BOOST_FOREACH (string_pair const& ex, this->local_expressions) {
    std::vector<value_type> av;
    for (int i=0; i< (has_bond_operator(ex.second) ? num_bonds() : num_sites());++i){
      cout << iter->name() << " " << iter->description() << " " << iter->value() << endl;
      av.push_back(real(iter++->value()));   
    }
    this->local_values[ex.first].push_back(av);
  }
  
  // average measurements will be identical loops, but all terms added instead of stored separately

   BOOST_FOREACH (string_pair const& ex, this->average_expressions) {
    cout << iter->name() << " " << iter->description() << " " << iter->value() << endl;
    this->average_values[ex.first].push_back(real(iter++->value()));
  }
  
  // correlations
  BOOST_FOREACH (string_string_pair_pair const& ex, this->correlation_expressions) {
    std::vector<value_type> av;
    for (int i=0; i<num_distances();++i){
      cout << iter->name() << " " << iter->description() << " " << iter->value() << endl;
      av.push_back(real(iter++->value()));   
    }
    this->correlation_values[ex.first].push_back(av);
  }

  if (iter != S.corr.end())
    std::cerr << "Did not get right number of measurements\n";
}

// SAVE_MPS: assemble the ground state of the final sweep as an MPS.
//
// The last step of final_sweep (with rotation) leaves the wave function on
// left block n = L/2-1, sites n and n+1, and right block L-n-2 in
// gs_<name>_<n>_l, and the transformation that built left (right) block k
// from block k-1 and one site in rho_<name>_<k>_l (_r). The basis stored
// with each transformation is that of its rows, labelled (old block, site)
// on the left and (site, old block) on the right. Left tensors are left
// orthonormal, right tensors right orthonormal, and the two-site center is
// split by SVD, so the MPS is in mixed canonical form with its center on
// site n+1. Local states are in the order of the ALPS site basis.
template<class value_type>
void
DMRGTask<value_type>::extract_mps()
{
  dmtk::System<value_type > &S = this->system;
  const int L = num_sites();
  if (L < 4)
    boost::throw_exception(std::runtime_error("SAVE_MPS needs at least 4 sites"));
  const int n = L/2 - 1;

  // dmtk sorts each site basis by quantum numbers; State::i1 keeps the ALPS index
  std::vector<std::vector<std::size_t> > alps_index(L);
  for (int i = 0; i < L; ++i) {
    const dmtk::Basis &b = site_block[site_type(i)].basis();
    for (std::size_t p = 0; p < b.size(); ++p)
      alps_index[i].push_back(b[p].i1);
  }
  std::vector<std::vector<value_type> > tensors(L);
  std::vector<std::vector<std::size_t> > shapes(L, std::vector<std::size_t>(3));

  // the single-site blocks at either end: bond index = dmtk site state
  for (int i = 0; i < L; i += L-1) {
    std::size_t d = alps_index[i].size();
    shapes[i][0] = (i == 0 ? 1 : d); shapes[i][1] = d; shapes[i][2] = (i == 0 ? d : 1);
    tensors[i].assign(d*d, value_type(0));
    for (std::size_t p = 0; p < d; ++p) {
      if (i == 0)
        tensors[i][alps_index[i][p]*d + p] = 1.;   // [1][s][bond]
      else
        tensors[i][p*d + alps_index[i][p]] = 1.;   // [bond][s][1]
    }
  }

  // block k on the left adds site k-1, block k on the right adds site L-k
  for (int k = 2; k <= std::max(n, L-n-2); ++k)
    for (int side = 0; side < 2; ++side) {
      int position = (side == 0 ? dmtk::LEFT : dmtk::RIGHT);
      int site = (side == 0 ? k-1 : L-k);
      if ((side == 0 && k > n) || (side == 1 && k > L-n-2))
        continue;
      dmtk::BMatrix<value_type > u;
      dmtk::Basis rows;
      S.read_rho(u, rows, k, position);
      std::size_t old_dim = (side == 0 ? shapes[site-1][2] : shapes[site+1][0]);
      std::size_t d = alps_index[site].size();
      std::size_t new_dim = 0;
      for (typename dmtk::BMatrix<value_type>::const_iterator it = u.begin(); it != u.end(); ++it)
        new_dim = std::max<std::size_t>(new_dim, it->col_range().end()+1);
      if (side == 0) {
        shapes[site][0] = old_dim; shapes[site][1] = d; shapes[site][2] = new_dim;
      } else {
        shapes[site][0] = new_dim; shapes[site][1] = d; shapes[site][2] = old_dim;
      }
      std::vector<value_type> &t = tensors[site];
      t.assign(old_dim*d*new_dim, value_type(0));
      for (typename dmtk::BMatrix<value_type>::const_iterator it = u.begin(); it != u.end(); ++it) {
        std::size_t r0 = it->row_range().begin(), c0 = it->col_range().begin();
        for (std::size_t i = 0; i < it->row_range().size(); ++i) {
          const dmtk::State &st = rows[r0+i];
          std::size_t old_state = (side == 0 ? st.i1 : st.i2);
          std::size_t s = alps_index[site][side == 0 ? st.i2 : st.i1];
          if (old_state >= old_dim)
            boost::throw_exception(std::runtime_error("SAVE_MPS: block transformations of the final sweep do not chain"));
          for (std::size_t j = 0; j < it->col_range().size(); ++j) {
            std::size_t c = c0 + j;
            if (side == 0)
              t[(old_state*d + s)*new_dim + c] = (*it)(j, i);
            else
              t[(c*d + s)*old_dim + old_state] = (*it)(j, i);
          }
        }
      }
    }

  // two-site wave function of the last step
  char file[255];
  snprintf(file, sizeof(file), "gs_%s_%i_l.dat", "ALPS", n);
  std::ifstream in(dmtk::tmp_files.get_read_filename(file), std::ios::in|std::ios::binary);
  if (!in)
    boost::throw_exception(std::runtime_error(std::string("SAVE_MPS: could not read ") + file));
  dmtk::VectorState<value_type > gs;
  gs.read(in);
  // read() rebuilds the subspaces only when a quantum number is conserved,
  // but the values are stored by subspace in either case (cf. System::read_gs)
  gs.resize(gs.qn_mask());
  std::size_t da = shapes[n-1][2], ds = alps_index[n].size(), dt = alps_index[n+1].size(), db = shapes[n+2][0];
  if (gs.b1().dim() != da || gs.b2().dim() != ds || gs.b3().dim() != dt || gs.b4().dim() != db)
    boost::throw_exception(std::runtime_error(std::string("SAVE_MPS: ") + file + " does not match the blocks of the final sweep"));
  std::vector<value_type> psi(da*ds*dt*db, value_type(0));
  for (typename dmtk::VectorState<value_type>::const_iterator sp = gs.subspace_begin(); sp != gs.subspace_end(); ++sp) {
    const dmtk::StateSpace &ss = *sp;
    std::size_t idx = ss.start();
    dmtk::SubSpace r1 = ss[1], r2 = ss[2], r3 = ss[3], r4 = ss[4];
    for (std::size_t a = r1.begin(); a < r1.begin() + r1.dim(); ++a)
      for (std::size_t s = r2.begin(); s < r2.begin() + r2.dim(); ++s)
        for (std::size_t t = r3.begin(); t < r3.begin() + r3.dim(); ++t)
          for (std::size_t b = r4.begin(); b < r4.begin() + r4.dim(); ++b)
            psi[((a*ds + alps_index[n][s])*dt + alps_index[n+1][t])*db + b] = gs[idx++];
  }
  dmrg_mps::split_two_site(psi, da, ds, dt, db, tensors[n], shapes[n], tensors[n+1], shapes[n+1]);

  mps_tensors.swap(tensors);
  mps_shapes.swap(shapes);
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
  
//  ar["spectrum"] << static_cast<const alps::EigenvectorMeasurements<value_type >&>(*this);

  typedef typename std::map<std::string,std::vector<double> >::const_iterator IT;
  for (IT it=iteration_measurements.begin(); it != iteration_measurements.end();++it)
      ar["simulation/results/Iteration "+alps::hdf5_name_encode(it->first)+"/mean/value"] << it->second;

  if (!mps_tensors.empty()) {
    // tensors are row-major [left bond][local state][right bond]
    ar["mps/length"] << int(mps_tensors.size());
    ar["mps/center"] << int(mps_tensors.size()/2);
    std::vector<int> types;
    for (std::size_t i = 0; i < mps_tensors.size(); ++i) {
      std::string site = "mps/tensors/" + boost::lexical_cast<std::string>(i);
      ar[site + "/shape"] << mps_shapes[i];
      ar[site + "/values"] << mps_tensors[i];
      types.push_back(site_type(i));
    }
    ar["mps/site_type"] << types;
    // quantum numbers of the local states, in the order of the ALPS site basis
    for (int type = 0; type <= alps::maximum_vertex_type(graph()); ++type) {
      alps::site_basis<short> b(site_basis(type));
      for (std::size_t q = 0; q < site_basis(type).size(); ++q) {
        std::vector<double> values;
        for (std::size_t s = 0; s < b.size(); ++s)
          values.push_back(b[s][q].to_double());
        ar["mps/site_basis/" + boost::lexical_cast<std::string>(type) + "/"
           + alps::hdf5_name_encode(site_basis(type)[q].name())] << values;
      }
    }
  }
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


template<class value_type>
dmtk::BasicOp<value_type > 
DMRGTask<value_type>::create_site_operator(std::string const& name, alps::SiteOperator const& siteop, int type)
{ 
  dmtk::Block<value_type > &block = site_block[type];
  dmtk::BasicOp<value_type > *op = block(name.c_str(),0);
  if(!op) {
    // we add new operator to single-site block 
    boost::multi_array<value_type,2> orig = alps::get_matrix(value_type(),siteop,site_basis(type),parms,true);
    dmtk::Matrix<value_type > dest(orig.shape()[0],orig.shape()[1]);
    dmtk::QN dqn;
    bool first = true;
    for(size_t i = 0; i < orig.shape()[0]; i++){
      for(size_t j = 0; j < orig.shape()[1]; j++){
        dest[j][i] = orig[i][j];
        if(abs(dest[j][i]) > 1.e-10){
          dmtk::QN new_dqn;
          new_dqn = block.basis()[j].qn() - block.basis()[i].qn();      
          if(!first && new_dqn != dqn) 
            boost::throw_exception(std::runtime_error("Matrix elements are inconsistent with change in quantum numbers"));
          dqn = new_dqn;
          first = false;
        } 
      }
    }
    dmtk::BasicOp<value_type > new_op(name.c_str(),0);
    new_op.dqn = dqn;
    new_op.resize(block.basis());
    new_op = dest;
    new_op.set_fermion(dqn.fermion_sign() == -1);
    block.push_back(new_op);
    return new_op;
  }

  return *op;
}
    
    
template<class value_type>
void
DMRGTask<value_type>::build_site_operator(alps::SiteOperator const& siteop, int site, dmtk::Hami<value_type > &this_hami)
{
  typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator> > V;
  V  ops = siteop.template templated_split<value_type>();
  int type = site_type(site);
  alps::expression::ParameterEvaluator<value_type> coords(coordinate_as_parameter(site));
  for (typename V::iterator it=ops.begin(); it!=ops.end();++it) {
    std::string name = simplify_name(it->template get<1>());
    dmtk::BasicOp<value_type> op1(name.c_str(),site);
    op1 = create_site_operator(name,it->template get<1>(),type);
    op1.set_site(site);
    dmtk::Term<value_type> real_t = op1;
    it->template get<0>().partial_evaluate(coords);
    real_t.coef() = it->template get<0>().value();
    this_hami += real_t;
  }
}

template<class value_type>
void
DMRGTask<value_type>::build_2site_operator(std::pair<alps::SiteOperator,alps::SiteOperator> const& siteops, 
                            std::pair<int,int> sites, dmtk::Hami<value_type > &this_hami)
{
  typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator> > V;
  V  ops1 = siteops.first.template templated_split<value_type>();
  V  ops2 = siteops.second.template templated_split<value_type>();
  for (typename V::const_iterator tit1=ops1.begin(); tit1!=ops1.end();++tit1)
    for (typename V::const_iterator tit2=ops2.begin(); tit2!=ops2.end();++tit2) {
      std::string s_name1 = simplify_name(tit1->template get<1>());
      std::string s_name2 = simplify_name(tit2->template get<1>());
      dmtk::BasicOp<value_type > op1(s_name1.c_str(),sites.first); 
      dmtk::BasicOp<value_type > op2(s_name2.c_str(),sites.second); 
      op1 = create_site_operator(s_name1,tit1->template get<1>(),site_type(sites.first));
      op2 = create_site_operator(s_name2,tit2->template get<1>(),site_type(sites.second));
      op1.set_site(sites.first);
      op2.set_site(sites.second);
      dmtk::Term<value_type > real_t = op1*op2;
      real_t.coef() = tit1->template get<0>().value()*tit2->template get<0>().value();
      this_hami += real_t;
   }
}


template<class value_type>
void
DMRGTask<value_type>::build_bond_operator(alps::BondOperator const& bondop, bond_descriptor const& b, dmtk::Hami<value_type > &this_hami)
{
  typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator,alps::SiteOperator > > V;
  alps::expression::ParameterEvaluator<value_type> coords(coordinate_as_parameter(b));
  alps::SiteBasisDescriptor<short> b1 = basis().site_basis(site_type(source(b)));
  alps::SiteBasisDescriptor<short> b2 = basis().site_basis(site_type(target(b)));
  
  V  ops = bondop.template templated_split<value_type>(b1,b2);
  for (typename V::iterator tit=ops.begin(); tit!=ops.end();++tit) {
    std::string s_name1 = simplify_name(tit->template get<1>());
    std::string s_name2 = simplify_name(tit->template get<2>());
    dmtk::BasicOp<value_type > op1(s_name1.c_str(),source(b)); 
    dmtk::BasicOp<value_type > op2(s_name2.c_str(),target(b)); 
    op1 = create_site_operator(s_name1,tit->template get<1>(),site_type(source(b)));
    op2 = create_site_operator(s_name2,tit->template get<2>(),site_type(target(b)));
    op1.set_site(source(b));
    op2.set_site(target(b));
    tit->template get<0>().partial_evaluate(coords);
    dmtk::Term<value_type > real_t = op1*op2;
    if (s_name1=="0")    
      real_t = op2;
    else if (s_name2=="0")
      real_t = op1;
    real_t.coef() = tit->template get<0>().value();
    this_hami += real_t;
  }
}

// Split an n-site term into products of site operators (nsiteoperator.h).
// As in exact diagonalization, a term that changes a conserved quantum number
// or has an odd number of fermion operators is refused (building the local
// table checks both), and one that changes a quantity capped by a cutoff
// (Nmax) warns.
template<class value_type>
std::vector<alps::nsite_product<value_type> >
DMRGTask<value_type>::split_nsite_operator(alps::NSiteOperator const& op, alps::nsite_group const& p)
{
  std::vector<alps::SiteBasisDescriptor<short> > bases;
  for (std::size_t r = 0; r < p.sites.size(); ++r)
    bases.push_back(basis().site_basis(site_type(p.sites[r])));
  alps::nsite_matrix<value_type> mat = op.template matrix<value_type>(bases,parms);
  for (std::size_t i = 0; i < quantumnumber_names.size(); ++i)
    if (conserved_quantumnumber[i] && mat.nonconserved.count(quantumnumber_names[i]))
      boost::throw_exception(std::runtime_error("<NSITETERM> on n-site groups of type "
        + boost::lexical_cast<std::string>(p.type) + " changes the total " + quantumnumber_names[i]
        + ", which is in CONSERVED_QUANTUMNUMBERS. Remove it from CONSERVED_QUANTUMNUMBERS or fix the term."));
  if (!warned_cutoff && !mat.nonconserved.empty()) {
    std::string warning = alps::nsite_cutoff_warning(mat.nonconserved,bases,parms);
    if (!warning.empty()) {
      std::cerr << warning;
      warned_cutoff = true;
    }
  }
  return op.template split<value_type>(bases,parms);
}

// One dmtk term per product of the n-site term: coef * O_0 * O_1 * ... in
// role order, one site operator per site of the group. dmtk reorders the
// operators by site with the fermion sign (Hami::reorder_terms).
template<class value_type>
void
DMRGTask<value_type>::build_nsite_operator(std::vector<alps::nsite_product<value_type> > const& products,
                                               alps::NSiteOperator const& op, alps::nsite_group const& p, dmtk::Hami<value_type > &this_hami)
{
  for (std::size_t t = 0; t < products.size(); ++t) {
    dmtk::Term<value_type > real_t;
    bool zero = false;
    for (std::size_t r = 0; r < p.sites.size() && !zero; ++r) {
      if (products[t].ops[r].empty())
        continue;
      int type = site_type(p.sites[r]);
      alps::SiteOperator siteop(products[t].ops[r],op.sites()[r]);
      std::string name = simplify_name(siteop);
      std::pair<int,std::string> key(type,name);
      if (zero_site_operator.find(key) == zero_site_operator.end()) {
        boost::multi_array<value_type,2> m = alps::get_matrix(value_type(),siteop,site_basis(type),parms,true);
        bool z = true;
        for (std::size_t i = 0; i < m.shape()[0] && z; ++i)
          for (std::size_t j = 0; j < m.shape()[1] && z; ++j)
            if (std::abs(m[i][j]) > 1.e-10)
              z = false;
        zero_site_operator[key] = z;
      }
      if (zero_site_operator[key]) {
        zero = true;
        break;
      }
      dmtk::BasicOp<value_type > op1 = create_site_operator(name,siteop,type);
      op1.set_site(p.sites[r]);
      real_t = real_t*op1;
    }
    if (zero || real_t.size() == 0)
      continue;
    real_t.coef() = products[t].coef;
    this_hami += real_t;
  }
}
