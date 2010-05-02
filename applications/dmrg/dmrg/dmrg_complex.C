/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*                             Matthias Troyer <troyer@itp.phys.ethz.ch>
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

/* $Id: dmrg.C 2896 2008-07-26 15:56:05Z troyer $ */


#define WITH_LAPACK
//#define WITH_WARNINGS
#include <dmtk/dmtk.h>
#include "dmrg_complex.h"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/archive/tmpdir.hpp>
#include <alps/numeric/real.hpp>
#include <alps/scheduler.h>
#include <boost/foreach.hpp>

DMRGTask::DMRGTask(const alps::ProcessList& w,const boost::filesystem::path& fn)
  : alps::scheduler::Task(w,fn)
  , alps::graph_helper<>(parms) 
  , alps::model_helper<>(*this,parms)
  , alps::EigenvectorMeasurements<value_type >(*this)
{
  init();
}

DMRGTask::DMRGTask(const alps::ProcessList& w,const alps::Parameters& p)
  : alps::scheduler::Task(w,p) 
  , alps::graph_helper<>(parms) 
  , alps::model_helper<>(*this,parms)
  , alps::EigenvectorMeasurements<value_type >(*this)
{
  init();
}


void DMRGTask::init()
{
  if (parms.defined("TEMP_DIRECTORY")) {
    std::string temp_dir = parms["TEMP_DIRECTORY"];
    dmtk::tmp_files.set_temp_dir(temp_dir.c_str());
  } else {
    dmtk::tmp_files.set_temp_dir(boost::archive::tmpdir());
  }

  num_eigenvalues = this->parms.value_or_default("NUMBER_EIGENVALUES",1);
   
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" ,");

  dmtk::QN::init();
  // read number of sweeps
  num_sweeps = parms.value_or_default("SWEEPS",4);
  start_sweep = parms.value_or_default("START_SWEEP",0);
  start_dir = parms.value_or_default("START_DIR",0);
  if(start_dir != 0 && start_dir != 1){
    boost::throw_exception(std::runtime_error("START_DIR can assume the values 0 (left-to-right) or 1 (right-to-left)"));
  }
  start_iter = parms.value_or_default("START_ITER",1);
  verbose = parms.value_or_default("VERBOSE",0);

  // read number of states
  nwarmup = 20;
  if (parms.defined("NUM_WARMUP_STATES")) {
    nwarmup = static_cast<int>(parms["NUM_WARMUP_STATES"]);
  }
  if (parms.defined("STATES")) {
    std::string states_string = parms["STATES"];
    tokenizer state_tokens(states_string, sep);
    for (tokenizer::const_iterator it = state_tokens.begin(); it !=state_tokens.end();++it)
      num_states.push_back(boost::lexical_cast<int>(*it));
    if (num_states.size() < 2*num_sweeps)
      boost::throw_exception(std::runtime_error("Need to specify either 2*SWEEPS different values in STATES, or one MAXSTATES value"));
  }
  else if (parms.defined("MAXSTATES")) {
    int maxstates = static_cast<int>(parms["MAXSTATES"]);
    for (int i=0; i< 2*num_sweeps; ++i)
      num_states.push_back((i+1)*maxstates/(2*num_sweeps));
  }
  else if (parms.defined("NUMSTATES")) {
    int maxstates = static_cast<int>(parms["NUMSTATES"]);
    for (int i=0; i< 2*num_sweeps; ++i)
      num_states.push_back(maxstates);
  }
  else 
    boost::throw_exception(std::runtime_error("Need to specify either 2*SWEEPS different values in STATES, or one MAXSTATES value"));
  
  // read quantum numbers and their total values
  quantumnumber_names.resize(0);
  conserved_quantumnumber.resize(0);
  conserved_quantumnumber_value.resize(0);
  if (parms.defined("CONSERVED_QUANTUMNUMBERS")) {
    std::string qn_string = parms["CONSERVED_QUANTUMNUMBERS"];
    tokenizer qn_tokens(qn_string, sep);
    std::copy(qn_tokens.begin(),qn_tokens.end(),std::back_inserter(quantumnumber_names));
    conserved_quantumnumber.resize(quantumnumber_names.size(),false);
    conserved_quantumnumber_value.resize(quantumnumber_names.size());
    for (int i=0; i<quantumnumber_names.size(); i++) {
      if (parms.defined(quantumnumber_names[i]+"_total")) {
        conserved_quantumnumber[i] = true;
        conserved_quantumnumber_value[i] = alps::evaluate<double>(static_cast<std::string>(parms[quantumnumber_names[i]+"_total"]),parms);
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

void DMRGTask::dostep() 
{
  if (finished()) 
    return;
  
  dmtk::Lattice l(num_sites(),dmtk::OBC);
  hami = dmtk::Hami<value_type >(l);
  site_block.resize(alps::maximum_vertex_type(graph())+1);
  
//  Iterating over sites: create site blocks
  for (site_iterator it = sites().first ; it != sites().second ; ++it) {
    hami.sites(*it) = &site_block[site_type(*it)];
    site_block[site_type(*it)].clear();
  }

  // iterate over all ste types
  for (int type  = 0 ; type < alps::maximum_vertex_type(graph())+1 ; ++type) {
    // create quantum numbers for this site block
    if(conserved_quantumnumber.size() > 0) {
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

//-------------------------------------------------------------

  if (verbose)
    std::cout << hami.description() << endl;

  // set up measurements
  
  typedef std::pair<std::string,std::string> string_pair;
 
  dmtk::Hami<value_type > meas_terms(hami);
  meas_terms.clear();

  // calculate local measurements
  BOOST_FOREACH (string_pair const& ex, local_expressions) {
    if (has_bond_operator(ex.second)) {
      int i=0;
      for (bond_iterator bit=bonds().first; bit!=bonds().second;++bit,++i) {
        dmtk::Hami<value_type > meas;
        build_bond_operator(get_bond_operator(ex.second,parms),*bit,meas);
        meas.set_name((ex.first+"["+boost::lexical_cast<std::string>(i)+"]").c_str());
//meas.set_name(meas.description().c_str());
        meas_terms += dmtk::BasicOp<value_type >(meas);
      }
    }
    else {
      for (site_iterator sit=sites().first; sit!=sites().second;++sit) {
        dmtk::Hami<value_type > meas;
        build_site_operator(make_site_term(ex.second),*sit,meas);
        meas.set_name((ex.first+"["+boost::lexical_cast<std::string>(*sit)+"]").c_str());
//meas.set_name(meas.description().c_str());
        meas_terms += dmtk::BasicOp<value_type >(meas);
      }
    }
  }
  
  // average measurements will be identical loops, but all terms added instead of stored separately

   BOOST_FOREACH (string_pair const& ex, average_expressions) {
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
    std::vector<dmtk::Hami<value_type > > corr_meas(num_distances());

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
          meas.set_name((ex.first+"["+boost::lexical_cast<std::string>(d)+"]").c_str());
          meas *= 1./double(distance_mult[d]);
          corr_meas[d] += meas;
        }
        else {
          // use ops1 and ops2
          dmtk::Hami<value_type > meas;
          build_2site_operator(std::make_pair(ops1,ops2),std::pair<int,int>(*sit1,*sit2),meas);
          meas *= 1./double(distance_mult[d]);
          corr_meas[d] += meas;
        }
    }
    for (unsigned d=0; d<num_distances();++d) {
      corr_meas[d].set_name((ex.first+"["+boost::lexical_cast<std::string>(d)+"]").c_str());
      meas_terms += dmtk::BasicOp<value_type >(corr_meas[d]);
    }
  }

  dmtk::Hami<value_type >::iterator iter;
  int i = 0;
  
  if (verbose) {
    for(iter = meas_terms.begin(); iter != meas_terms.end(); i++, iter++)
      cout << i << " " << iter->name() << " " << iter->description() << endl;

    std::cout << meas_terms.description() << "\n";
  }
  
///////////////////////////////////////////////////////////////
// Simulation
///////////////////////////////////////////////////////////////
  hami = hami.reorder_terms();
  if(verbose)
    cout << hami.description() << endl;
  this->system = dmtk::System<value_type >(hami,l,"ALPS");
  dmtk::System<value_type > &S = this->system;
  
  S.set_calc_gap(num_eigenvalues-1); 
  dmtk::Matrix<size_t> nstates(2,num_sweeps);
  dmtk::QN qn;
  qnmask = 0;
  for(int i = 0; i < quantumnumber_names.size(); i++) {
    if(conserved_quantumnumber[i]){
      qn[quantumnumber_names[i]] = conserved_quantumnumber_value[i];
      qnmask |= (1 << (dmtk::QN::get_qn_index(quantumnumber_names[i])));
    }
  }
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
  if(start_sweep != 0) {
    S.resume(start_sweep, start_dir, start_iter);
  } else {
    S.run(nwarmup);
  }
  S.corr = meas_terms;
  S.final_sweep(num_states[num_states.size()-1], dmtk::RIGHT2LEFT, 1, true); 
  S.measure();
  save_results();
  for(int i = 1; i < S._target.size(); i++){
    S.gs = S._target[i];
    S.measure();
    save_results();
  }
  for(int i = 0; i < S.energy.size(); i++) {
    average_values["Energy"].push_back(S.energy[i]);
  }
  average_values["Truncation error"].push_back(S.truncation_error());
  finish();
}


void
DMRGTask::save_results()
{
  dmtk::System<value_type > &S = this->system;
  dmtk::Hami<value_type >::iterator iter = S.corr.begin();
  typedef std::pair<std::string,std::string> string_pair;
  typedef std::pair<std::string,std::pair<std::string,std::string> > string_string_pair_pair;

  // store local measurements
  BOOST_FOREACH (string_pair const& ex, local_expressions) {
    std::vector<value_type> av;
    for (int i=0; i< (has_bond_operator(ex.second) ? num_bonds() : num_sites());++i)
      av.push_back(iter++->value().real());   
    local_values[ex.first].push_back(av);
  }
  
  // average measurements will be identical loops, but all terms added instead of stored separately

   BOOST_FOREACH (string_pair const& ex, average_expressions) {
    average_values[ex.first].push_back(iter++->value().real());
  }
  
  // correlations
  BOOST_FOREACH (string_string_pair_pair const& ex, this->correlation_expressions) {
    std::vector<value_type> av;
    for (int i=0; i<num_distances();++i)
      av.push_back(iter++->value().real());   
    correlation_values[ex.first].push_back(av);
  }

  if (iter != S.corr.end())
    std::cerr << "Did not get right number of measurements\n";
}
    
#ifdef ALPS_HAVE_HDF5
void DMRGTask::serialize(alps::hdf5::oarchive & ar) const
{
  alps::scheduler::Task::serialize(ar);
  std::map<std::string,std::vector<value_type> >::const_iterator it = average_values.find("Energy");
  if (it != average_values.end()) {
    std::vector<double> energies = alps::numeric::real(it->second);
    ar << alps::make_pvp("spectrum/energies",energies);
  }
  ar << alps::make_pvp("spectrum",static_cast<const alps::EigenvectorMeasurements<std::complex<double> >&>(*this));

  for (std::size_t i=0;i<iteration_measurements.size();++i) {
    ar << alps::make_pvp("spectrum/iteration/"+boost::lexical_cast<std::string>(i),iteration_measurements[i]);
    std::map<std::string,std::vector<value_type> >::const_iterator it = iteration_measurements[i].average_values.find("Energy");
    if (it != average_values.end())
      ar << alps::make_pvp("spectrum/iteration/"+boost::lexical_cast<std::string>(i)+"/energies",it->second);
  }
}
#endif

void DMRGTask::write_xml_body(alps::oxstream& out, const boost::filesystem::path& p, bool writeallxml) const
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


dmtk::BasicOp<DMRGTask::value_type > DMRGTask::create_site_operator(std::string const& name, alps::SiteOperator const& siteop, int type)
{ 
  dmtk::Block<value_type > &block = site_block[type];
  dmtk::BasicOp<value_type > *op = block(name.c_str(),0);
  if(!op) {
    // we add new operator to single-site block 
    boost::multi_array<double,2> orig = alps::get_matrix(0.,siteop,site_basis(type),parms,true);
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
    block.push_back(new_op);
    return new_op;
  }

  return *op;
}
    
    
void
DMRGTask::build_site_operator(alps::SiteOperator const& siteop, int site, dmtk::Hami<value_type > &this_hami)
{
  typedef std::vector<boost::tuple<alps::Term,alps::SiteOperator> > V;
  V  ops = siteop.split();
  int type = site_type(site);
  alps::expression::ParameterEvaluator<value_type > coords(coordinate_as_parameter(site));
  for (V::iterator it=ops.begin(); it!=ops.end();++it) {
    std::string name = simplify_name(it->get<1>());
    dmtk::BasicOp<value_type> op1(name.c_str(),site);
    op1 = create_site_operator(name,it->get<1>(),type);
    op1.set_site(site);
    dmtk::Term<value_type> real_t = op1;
    it->get<0>().partial_evaluate(coords);
    real_t.coef() = it->get<0>().value().real();
    this_hami += real_t;
  }
}

void
DMRGTask::build_2site_operator(std::pair<alps::SiteOperator,alps::SiteOperator> const& siteops, 
                            std::pair<int,int> sites, dmtk::Hami<value_type > &this_hami)
{
  typedef std::vector<boost::tuple<alps::Term,alps::SiteOperator> > V;
  V  ops1 = siteops.first.split();
  V  ops2 = siteops.second.split();
  for (V::const_iterator tit1=ops1.begin(); tit1!=ops1.end();++tit1)
    for (V::const_iterator tit2=ops2.begin(); tit2!=ops2.end();++tit2) {
      std::string s_name1 = simplify_name(tit1->get<1>());
      std::string s_name2 = simplify_name(tit2->get<1>());
      dmtk::BasicOp<value_type > op1(s_name1.c_str(),sites.first); 
      dmtk::BasicOp<value_type > op2(s_name2.c_str(),sites.second); 
      op1 = create_site_operator(s_name1,tit1->get<1>(),site_type(sites.first));
      op2 = create_site_operator(s_name2,tit2->get<1>(),site_type(sites.second));
      op1.set_site(sites.first);
      op2.set_site(sites.second);
      dmtk::Term<value_type > real_t = op1*op2;
      real_t.coef() = tit1->get<0>().value().real()*tit2->get<0>().value().real();
      this_hami += real_t;
   }
}


void
DMRGTask::build_bond_operator(alps::BondOperator const& bondop, bond_descriptor const& b, dmtk::Hami<value_type > &this_hami)
{
  typedef std::vector<boost::tuple<alps::Term,alps::SiteOperator,alps::SiteOperator > > V;
  alps::expression::ParameterEvaluator<value_type > coords(coordinate_as_parameter(b));
  alps::SiteBasisDescriptor<short> b1 = basis().site_basis(source(b));
  alps::SiteBasisDescriptor<short> b2 = basis().site_basis(target(b));
  
  V  ops = bondop.split(b1,b2);
  for (V::iterator tit=ops.begin(); tit!=ops.end();++tit) {
    std::string s_name1 = simplify_name(tit->get<1>());
    std::string s_name2 = simplify_name(tit->get<2>());
    dmtk::BasicOp<value_type > op1(s_name1.c_str(),source(b)); 
    dmtk::BasicOp<value_type > op2(s_name2.c_str(),target(b)); 
    op1 = create_site_operator(s_name1,tit->get<1>(),site_type(source(b)));
    op2 = create_site_operator(s_name2,tit->get<2>(),site_type(target(b)));
    op1.set_site(source(b));
    op2.set_site(target(b));
    dmtk::Term<value_type > real_t = op1*op2;
    tit->get<0>().partial_evaluate(coords);
    real_t.coef() = tit->get<0>().value().real();
    this_hami += real_t;

  }
}


int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
   
   return alps::scheduler::start(argc,argv,alps::scheduler::SimpleFactory<DMRGTask>());

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
}

