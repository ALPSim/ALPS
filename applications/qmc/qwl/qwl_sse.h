/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2004-2005 by Stefan Wessel <wessel@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef ALPS_QWL_SSE_H
#define ALPS_QWL_SSE_H

#include "qwl_histogram.h"
#include <boost/lexical_cast.hpp>
#include <alps/scheduler/montecarlo.h>
#include <alps/alea.h>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>

using namespace std;
using namespace alps;

template<typename T> T sqr(T x) {return x*x;}


class QWL_SSE_Simulation : public scheduler::LatticeModelMCRun<> {
 public:

  QWL_SSE_Simulation(const ProcessList&,const Parameters&,int);
  void save(ODump&) const;
  void load(IDump&);
  void dostep();
  bool is_thermalized() const;
  double work_done() const;
  bool change_parameter(const std::string& name, const StringValue& value);
  static void print_copyright(std::ostream&);
  
 private:

  typedef boost::uint32_t statetype;
  typedef boost::uint32_t bond_state_type;
    
  unsigned int L;      // Size of the Operator String
  unsigned int non0;   // Number of non-unity operators
  double beta;     // Inverse Temperature
  double offset;
  
  unsigned int sweeps;
  
  histogram<double> g;
  histogram<double> histo;
  histogram<double> histoup;
  histogram<double> nmeasurements;
  histogram<double> uniform_structure_factor;
  histogram<double> staggered_structure_factor;
  histogram<double> transition_prob;

  unsigned int norder;
  unsigned int norder_min;
  unsigned int norder_max;
  boost::uint32_t includeLfactor;
  double logf;
  unsigned int logf_steps_total;
  unsigned int logf_step;
  double minimum_histogram;
  double flatness_treshold;
  boost::uint32_t use_zhou_bhatt;
  unsigned int block_sweeps_total;
  unsigned int block_sweeps;
  boost::uint32_t doing_multicanonical;
  boost::uint32_t upwalker;
  unsigned int thistime;
  boost::uint32_t all_done;
  boost::uint32_t measure_magnetics;

  struct vertextype {
    vertextype()
      : op(0) {
    };
    boost::uint32_t op;
    boost::uint32_t bond;
    bond_state_type vs;
    int         vv[4];
    boost::uint32_t  vn[4];
    boost::uint32_t vvisited[4];
  };
  
  vector<vertextype> operator_string;
  
  vector <statetype>  state;   
  vector <int>   statev;  //pointing to previous vertex
  vector <boost::uint32_t> staten;  //pointing to previous node number

  vector <boost::uint32_t> is_ferromagnetic;
  vector <double> matrix_factor;
  
  void init_measurements();
  void init_tables();
  void diagonal_update();
  void diagonal_update_multicanonical();
  void vertexbuild();
  void loopupdate();
  void measure();
  void store_histograms(string);
  boost::uint32_t exit_leg(vertextype& , boost::uint32_t);  
  double matrix_element(const int, const int, const int);

  statetype get_state(bond_state_type vs, int n) const {
    return ((vs>>n) & 1);
  }

  void set_statevector_in(bond_state_type& vs, statetype s0, statetype s1) {
    vs=s0 | (s1<<1);
  }

  void set_statevector_out(bond_state_type& vs, statetype s2, statetype s3) {
    vs=(vs | (s2<<2) ) | (s3<<3);
  }

  double decode_spin(statetype s) const {
    return ( (s) ? 0.5 : -0.5);
  }
  
  statetype random_state() const {
    return (random_01()>0.5 ? 1 : 0); 
  }
  
  void transform(bond_state_type& vs, int n)  {
    vs^=(1<<n);
  }

};

QWL_SSE_Simulation::QWL_SSE_Simulation(const ProcessList& where,const Parameters& p,int node)
  : scheduler::LatticeModelMCRun<>(where,p,node), 
    non0(0),
    sweeps(0),
    norder(0),
    logf_step(1),
    block_sweeps(0),
    doing_multicanonical(0),
    upwalker(1),
    thistime(0),
    all_done(0)
{
  state.resize(num_sites());
  statev.resize(num_sites());
  staten.resize(num_sites());
  for (site_iterator it=sites().first;it!=sites().second;++it) 
    state[*it]=random_state();
  beta=1;
  logf_steps_total=parms.value_or_default("NUMBER_OF_WANG_LANDAU_STEPS",16);
  norder_min=parms.value_or_default("EXPANSION_ORDER_MINIMUM",0);
  norder_max=parms.value_or_default("EXPANSION_ORDER_MAXIMUM",parms.value_or_default("CUTOFF",500));
  int g_left=norder_min;
  int g_size=norder_max-g_left+1;
  if (g_left<0)
    boost::throw_exception(std::runtime_error("EXPANSION_ORDER_MINIMUM must not be negative"));
  g.resize(g_size,g_left);
  histo.resize(g_size,g_left);
  histoup.resize(g_size,g_left);
  nmeasurements.resize(g_size,g_left);
  transition_prob.resize(g_size,g_left);
  L=parms.value_or_default("CUTOFF",norder_max);
  if (L<norder_max)
    boost::throw_exception(std::runtime_error("EXPANSION_ORDER_MAXIMUM must not exceed CUTOFF"));
  operator_string.resize(L);
  measure_magnetics=parms.value_or_default("MEASURE_MAGNETIC_PROPERTIES",1);
  if (measure_magnetics) {
    uniform_structure_factor.resize(g_size,g_left);
    if (is_bipartite()) staggered_structure_factor.resize(g_size,g_left);
  }    
  includeLfactor=parms.value_or_default("INCLUDE_COMBINATORICS_FACTORS",1);
  use_zhou_bhatt=parms.value_or_default("USE_ZHOU_BHATT_METHOD",1.);
  if (use_zhou_bhatt) {
    block_sweeps_total=g.size();
    logf=log(double(parms.value_or_default("INITIAL_MODIFICATION_FACTOR",exp(1.))));
    minimum_histogram=1./logf;
    flatness_treshold=parms.value_or_default("FLATNESS_TRESHOLD",1E10);
  }
  else {
    block_sweeps_total=parms.value_or_default("BLOCK_SWEEPS",10000);
    logf=log(double(parms.value_or_default("INITIAL_INCREASE_FACTOR",exp(L*log((double)num_sites())/block_sweeps_total))));
    minimum_histogram=0.;
    flatness_treshold=parms.value_or_default("FLATNESS_TRESHOLD",0.2);
  }
  init_tables();
  init_measurements();
}

void QWL_SSE_Simulation::save(ODump& dump) const {  
  g.save(dump);
  histo.save(dump);
  histoup.save(dump);
  nmeasurements.save(dump);
  transition_prob.save(dump);
  if (measure_magnetics) {
    uniform_structure_factor.save(dump);    
    staggered_structure_factor.save(dump);    
  }
  dump << sweeps;
  dump << L << non0;
  dump << logf << logf_step;
  dump << block_sweeps << block_sweeps_total;
  dump << norder;
  dump << minimum_histogram;
  dump << flatness_treshold;
  dump << doing_multicanonical << upwalker << thistime;
  dump << all_done;
  dump << state;
  for (unsigned int i=0;i<L;++i)
    dump << operator_string[i].op << operator_string[i].bond;    
}

void QWL_SSE_Simulation::load(IDump& dump) {
  g.load(dump);
  histo.load(dump);
  histoup.load(dump);
  nmeasurements.load(dump);
  transition_prob.load(dump);
  if (measure_magnetics) {
    uniform_structure_factor.load(dump);    
    staggered_structure_factor.load(dump);    
  }
  dump >> sweeps;
  dump >> L >> non0;
  dump >> logf >> logf_step;
  dump >> block_sweeps >> block_sweeps_total;
  dump >> norder;
  dump >> minimum_histogram;
  dump >> flatness_treshold;
  dump >> doing_multicanonical >> upwalker >> thistime;
  dump >> all_done;
  if(!where.empty()) { // skip this if we are just evaluating
    dump >> state;  
    operator_string.resize(L);
    for (unsigned int i=0;i<L;++i) 
      dump >> operator_string[i].op >> operator_string[i].bond;
    vertexbuild();
  }      
}

void QWL_SSE_Simulation::dostep() {
  ++sweeps; 
  ++block_sweeps;
  if (doing_multicanonical)
    diagonal_update_multicanonical();
  else
    diagonal_update();
  if (non0) {
    vertexbuild();
    loopupdate();
  }
  else
    for (site_iterator it=sites().first;it!=sites().second;++it) 
      state[*it]=random_state();
  if (doing_multicanonical) measure();
  if (block_sweeps>=block_sweeps_total) {
    if (doing_multicanonical) {
      if (!all_done) {
        block_sweeps=0;
        if (histo.min()>=minimum_histogram || parms.defined("SWEEPS")) {
          store_histograms("multicanonical");
          cerr << "[" << norder_min << "-" << norder_max << "]  all done." << endl;
          all_done=1;
        }
      }
    }
    else { 
      g.subtract();
      block_sweeps=0;
      double histomin=histo.min();
      double histoflatness=histo.flatness();
      cerr << "[" << norder_min << "-" << norder_max << "]  step "  << logf_step << ", ln[f]="<< logf 
           << " : flatness="<< histoflatness << " ratio=" << histomin/minimum_histogram <<endl;
      if (histoflatness<flatness_treshold && histomin>=minimum_histogram) {
        store_histograms(boost::lexical_cast<string>(logf_step));
        histo.fill(0.0);  
        if (logf_step==logf_steps_total) { 
          block_sweeps_total=parms.value_or_default("SWEEPS",block_sweeps_total);
          for (unsigned int i=g.left();i<g.right();++i)
            transition_prob[i]=exp(g[i]-g[i+1]);
          transition_prob[g.right()]=0.;
          cerr << "[" << norder_min << "-" << norder_max << "]  continuing using final weights..." << endl;
          doing_multicanonical=1;
        }
        else {
          ++logf_step; 
          logf/=2.0;
          if (use_zhou_bhatt) {
            block_sweeps_total=static_cast<unsigned int>(block_sweeps_total*1.41);
            minimum_histogram*=2;
          }
        }
      }
    }
  }
} 


bool QWL_SSE_Simulation::is_thermalized() const {
  return sweeps;
}


double QWL_SSE_Simulation::work_done() const {
  return all_done? 1. : ( doing_multicanonical ? (histo.min()>=minimum_histogram || parms.defined("SWEEPS") ? double(block_sweeps)/double(block_sweeps_total) : 0. ) : 0.);
}


bool QWL_SSE_Simulation::change_parameter(const std::string& name, const StringValue& value) {
  if(name=="SWEEPS")
    block_sweeps_total=static_cast<unsigned int>(value);
  else
    return false; // cannot do it
  return true; // could do it
}


void QWL_SSE_Simulation::diagonal_update() { // fixed length (vector) representation in WL
  unsigned int Lmnon0=L-non0;
  for (vector<vertextype>::iterator operator_str_iterator=operator_string.begin();operator_str_iterator!=operator_string.end();++operator_str_iterator) {
    switch(operator_str_iterator->op) {
    // insert staebchenspiel here
    case 0 : // identity
      if (norder<g.right()) {
        operator_str_iterator->bond=random_int(0,num_bonds()-1);
        bond_descriptor this_bond=bond(operator_str_iterator->bond);
        double probability=exp(g[norder]-g[norder+1]);
        probability*=matrix_element(inhomogeneous_bond_type(this_bond),state[source(this_bond)],state[target(this_bond)]);
        probability*=num_bonds();
        probability*=beta;
        if (includeLfactor) 
          probability/=Lmnon0;
        if (probability!=0. && ( (probability>=1) || (random_01()<probability) )) {
          operator_str_iterator->op=1;
          --Lmnon0;
          ++norder;
        }
      }
      break;
    case 1 : // diagonal
      if (norder>g.left()) {
        bond_descriptor this_bond=bond(operator_str_iterator->bond);
        double probability=exp(g[norder]-g[norder-1]);
        probability/=matrix_element(inhomogeneous_bond_type(this_bond),state[source(this_bond)],state[target(this_bond)]);
        probability/=num_bonds();
        probability/=beta;
        if (includeLfactor) 
          probability*=Lmnon0+1;
        if (probability!=0. && ( (probability>=1) || (random_01()<probability) )) {
            operator_str_iterator->op=0;
          ++Lmnon0;
          --norder;
        }
      }
      break;
    default : // off-diagonal
      state[source(bond(operator_str_iterator->bond))]=get_state(operator_str_iterator->vs,2);
      state[target(bond(operator_str_iterator->bond))]=get_state(operator_str_iterator->vs,3);
    };
    if (norder==g.right()) upwalker=0;
    else if (norder==g.left())  upwalker=1;
    g[norder]+=logf;
    ++histo[norder];
  };
  non0=L-Lmnon0;
}


void QWL_SSE_Simulation::diagonal_update_multicanonical() {
  unsigned int Lmnon0=L-non0;
  for (vector<vertextype>::iterator operator_str_iterator=operator_string.begin();operator_str_iterator!=operator_string.end();++operator_str_iterator) {
    switch(operator_str_iterator->op) {
    case 0 : 
      if (norder<g.right()) {
        operator_str_iterator->bond=random_int(0,num_bonds()-1);
        bond_descriptor this_bond=bond(operator_str_iterator->bond);
        double probability=transition_prob[norder];
        probability*=matrix_element(inhomogeneous_bond_type(this_bond),state[source(this_bond)],state[target(this_bond)]);
        probability*=num_bonds();
        probability*=beta;
        if (includeLfactor) 
          probability/=Lmnon0;
        if (probability!=0. && ( (probability>=1) || (random_01()<probability) )) {
          operator_str_iterator->op=1;
          --Lmnon0;
          ++norder;
        }
      }
      break;
    case 1 : 
      if (norder>g.left()) {
        bond_descriptor this_bond=bond(operator_str_iterator->bond);
        double probability=1./transition_prob[norder-1];
        probability/=matrix_element(inhomogeneous_bond_type(this_bond),state[source(this_bond)],state[target(this_bond)]);
        probability/=num_bonds();
        probability/=beta;
        if (includeLfactor) 
          probability*=Lmnon0+1;
        if (probability!=0. && ( (probability>=1) || (random_01()<probability) )) {
            operator_str_iterator->op=0;
          ++Lmnon0;
          --norder;
        }
      }
      break;
    default :
      state[source(bond(operator_str_iterator->bond))]=get_state(operator_str_iterator->vs,2);
      state[target(bond(operator_str_iterator->bond))]=get_state(operator_str_iterator->vs,3);
    };
    ++thistime;
    if (upwalker && norder==g.right()) {
      measurements["Time Up"]    << (double)thistime;
      measurements["Time Total"] << (double)thistime;
      thistime=0;
     upwalker=0;
    }
    else if (!upwalker && norder==g.left()) {
      measurements["Time Down"]  << (double)thistime;
      measurements["Time Total"] << (double)thistime;
      thistime=0;
      upwalker=1;    
    }
    ++histo[norder];
    if (upwalker) ++histoup[norder];
  };
  non0=L-Lmnon0;
}


void QWL_SSE_Simulation::vertexbuild() {
  for (unsigned int j=0;j<num_sites();++j) 
    statev[j]=-(j+1);
  unsigned int i=0;
  for (vector <vertextype>::iterator operator_str_iterator=operator_string.begin();operator_str_iterator!=operator_string.end();++operator_str_iterator) {
    if (operator_str_iterator->op) { 
      unsigned int bn0=source(bond(operator_str_iterator->bond));
      unsigned int bn1=target(bond(operator_str_iterator->bond));
      set_statevector_in(operator_str_iterator->vs,state[bn0],state[bn1]);
      operator_str_iterator->vv[0]=statev[bn0];
      operator_str_iterator->vv[1]=statev[bn1];
      operator_str_iterator->vn[0]=staten[bn0];
      operator_str_iterator->vn[1]=staten[bn1];
      if (statev[bn0] >= 0) { 
        operator_string[statev[bn0]].vv[staten[bn0]]=i;
        operator_string[statev[bn0]].vn[staten[bn0]]=0;
      }
      if (statev[bn1] >= 0) {
        operator_string[statev[bn1]].vv[staten[bn1]]=i;
        operator_string[statev[bn1]].vn[staten[bn1]]=1;
      }
      if (operator_str_iterator->op!=1) {
        if (operator_str_iterator->op == 2) {
          ++state[bn0];
          --state[bn1];
        }
        else {
          --state[bn0];
          ++state[bn1];
        }
      }
      set_statevector_out(operator_str_iterator->vs,state[bn0],state[bn1]);
      statev[bn0]=i;
      statev[bn1]=i;
      staten[bn0]=2;
      staten[bn1]=3;
    }
    ++i;
  };
  i=0;
  for (vector <vertextype>::iterator operator_str_iterator=operator_string.begin();operator_str_iterator!=operator_string.end();++operator_str_iterator) {
    if (operator_str_iterator->op) {
      int j=-(operator_str_iterator->vv[0]+1);
      if (j >= 0) {
        operator_str_iterator->vv[0]=statev[j];
        operator_str_iterator->vn[0]=staten[j];
        operator_string[statev[j]].vv[staten[j]]=i;
        operator_string[statev[j]].vn[staten[j]]=0;
      }
      j=-(operator_string[i].vv[1]+1);
      if (j >= 0) {
        operator_str_iterator->vv[1]=statev[j]; 
        operator_str_iterator->vn[1]=staten[j];
        operator_string[statev[j]].vv[staten[j]]=i;
        operator_string[statev[j]].vn[staten[j]]=1;
      }
      for (int l=0;l<4;++l)
        operator_str_iterator->vvisited[l]=0;
    }
    else
      for (int l=0;l<4;++l)
        operator_str_iterator->vvisited[l]=1;
    ++i;
  };
}
  
void QWL_SSE_Simulation::loopupdate() {
  static unsigned int i0;
  static unsigned int ir;
  static boost::uint32_t n0;
  static boost::uint32_t nr;
  static boost::uint32_t nex;
  i0=0;
  n0=0;
  int flipflag;
  do {
    while ( ( (i0!=L) ?  operator_string[i0].vvisited[n0] : 0) ) {
      if (n0==3) {
        ++i0;
        n0=0;
      }
      else 
       ++n0;
    };
    if (i0==L) break;
    flipflag=(random_01()>0.5);
    ir=i0;
    nr=n0;
    do {
      operator_string[ir].vvisited[nr]=1;
      if (flipflag)
        transform(operator_string[ir].vs,nr);
      nex=exit_leg(operator_string[ir],nr);
      operator_string[ir].vvisited[nex]=1;
      if (flipflag)
        transform(operator_string[ir].vs,nex);
      if (ir==i0 && nex==n0)
        break;
      nr=operator_string[ir].vn[nex];
      ir=operator_string[ir].vv[nex];
    } while((ir!=i0)||(nr!=n0));
  } while (1);
  for (vector <vertextype>::iterator operator_str_iterator=operator_string.begin();operator_str_iterator!=operator_string.end();++operator_str_iterator)
    if (operator_str_iterator->op)
      switch (get_state(operator_str_iterator->vs,2)-get_state(operator_str_iterator->vs,0)) {
        case 0:
          operator_str_iterator->op=1;
          break;
        case 1:
          operator_str_iterator->op=2;
          break;
        default:
          operator_str_iterator->op=3;
      }
  for (unsigned int j=0;j<num_sites();++j)
    state[j]= (statev[j]>=0) ? get_state(operator_string[statev[j]].vs,staten[j]) : random_state();
}


inline boost::uint32_t QWL_SSE_Simulation::exit_leg(vertextype&  opn, boost::uint32_t nr) {
  if (is_ferromagnetic[inhomogeneous_bond_type(bond(opn.bond))])
    return 3-nr;
  if (!nr)
    return 1;
  if (nr==1)
    return 0;
  if (nr==2)
    return 3;
  return 2;   
} 


inline double QWL_SSE_Simulation::matrix_element(const int b_t, const int s0, const int s1) {
  if (is_ferromagnetic[b_t])
    return ( (s0==s1) ? matrix_factor[b_t] : 0. );
  return ( (s0==s1) ? 0. : matrix_factor[b_t] );
}


void QWL_SSE_Simulation::init_tables() {
  int max_bond_type=0;
  for (bond_iterator it=bonds().first; it!=bonds().second;++it) {
    int this_bond_type=inhomogeneous_bond_type(*it);
    if (this_bond_type>=max_bond_type) max_bond_type=this_bond_type;
  }
  is_ferromagnetic.resize(max_bond_type+1);
  matrix_factor.resize(max_bond_type+1);
  offset=0;
  for (bond_iterator it=bonds().first; it!=bonds().second;++it) {
    Parameters p(parms);
    if (inhomogeneous()) 
      throw_if_xyz_defined(parms,*it);
    if (inhomogeneous_sites()) {
      p << coordinate_as_parameter(source(*it));
      p << coordinate_as_parameter(target(*it));
    }
    if (inhomogeneous_bonds()) 
      p << coordinate_as_parameter(*it);
    int this_bond_type=inhomogeneous_bond_type(*it);
    int source_site_type=inhomogeneous_site_type(source(*it));
    int target_site_type=inhomogeneous_site_type(target(*it));
    int source_num_states = model().basis().site_basis(source_site_type).num_states();
    int target_num_states = model().basis().site_basis(target_site_type).num_states();
    if (source_num_states!=2 || target_num_states!=2) 
      boost::throw_exception(std::runtime_error("This model cannot be simulated with this code"));
    boost::multi_array<double,2> source_Szmatrix_symbolic =
      alps::get_matrix(double(),
      SiteOperator("Sz"),
      model().basis().site_basis(source_site_type),
      parms);
    boost::multi_array<double,2> target_Szmatrix_symbolic =
      alps::get_matrix(double(),
      SiteOperator("Sz"),
      model().basis().site_basis(target_site_type),
      parms);
    int source_spin_up_state=(source_Szmatrix_symbolic[0][0]>0.) ? 0 : 1;
    int target_spin_up_state=(target_Szmatrix_symbolic[0][0]>0.) ? 0 : 1;
    boost::multi_array<double,4> bondhamiltonian =
      alps::get_matrix(double(),
      model().bond_term(this_bond_type),
      model().basis().site_basis(source_site_type),
      model().basis().site_basis(target_site_type),
      parms);
    is_ferromagnetic[this_bond_type]=
      bondhamiltonian[source_spin_up_state][target_spin_up_state]
                     [source_spin_up_state][target_spin_up_state]<0? 1 : 0;
    matrix_factor[this_bond_type]=
      fabs(bondhamiltonian[source_spin_up_state][target_spin_up_state]
                          [source_spin_up_state][target_spin_up_state])*2.;
    offset+=
      fabs(bondhamiltonian[source_spin_up_state][target_spin_up_state]
                          [source_spin_up_state][target_spin_up_state]);
  }
}


void QWL_SSE_Simulation::init_measurements() {
  measurements << SimpleRealObservable("Offset");
  for (int pos=parms.value_or_default("START_STORING",logf_steps_total);pos<=logf_steps_total;++pos) {
    measurements << SimpleRealVectorObservable("Coefficients "+boost::lexical_cast<string>(pos));
    measurements << SimpleRealObservable("Total Sweeps "+boost::lexical_cast<string>(pos));
  }
  measurements << SimpleRealVectorObservable("Coefficients");
  measurements << SimpleRealVectorObservable("Histogram");
  measurements << SimpleRealVectorObservable("Fraction");
  measurements << SimpleRealObservable("Total Sweeps");
  measurements << RealObservable("Time Up");
  measurements << RealObservable("Time Down");
  measurements << RealObservable("Time Total");
  if (measure_magnetics) {
    measurements << SimpleRealVectorObservable("Uniform Structure Factor Coefficients");
    if (is_bipartite())
      measurements << SimpleRealVectorObservable("Staggered Structure Factor Coefficients");
  }   
}

 
void QWL_SSE_Simulation::measure() {
  if (measure_magnetics) {
    ++nmeasurements[norder];
    double mag=0;
    double smag=0;
    for (unsigned int i=0;i<num_sites();++i) {
      double local_state=decode_spin(state[i]);
      mag+=local_state;
      if (is_bipartite()) 
        smag+=parity(i)*local_state;
    }
    uniform_structure_factor[norder]+=sqr(mag);    
    if (is_bipartite()) 
      staggered_structure_factor[norder]+=sqr(smag);
  }   
}


void QWL_SSE_Simulation::store_histograms(string pos) {
  if (pos=="multicanonical") {
    valarray<double> gval=g.getvalarray(norder_min,norder_max);
    valarray<double> histoval=histo.getvalarray(norder_min,norder_max);
    valarray<double> histoupval=histoup.getvalarray(norder_min,norder_max);
    histoupval/=histoval;
    histoval/=block_sweeps_total;    
    double histovalsum=histoval.sum();
    histoval/=histovalsum;
    double red=gval[0]+log(histoval[0])-num_sites()*log(2.);
    for (int i=0;i<gval.size();++i)
      gval[i]+=log(histoval[i])-red;
    measurements["Coefficients"] << gval;
    measurements["Histogram"] << histoval;
    measurements["Fraction"] << histoupval;
    measurements["Total Sweeps"] << (double)sweeps; 
    measurements["Offset"] << offset;
    if (measure_magnetics) {
      valarray<double> nmeasurementsval=nmeasurements.getvalarray(norder_min,norder_max);
      valarray<double> uniform_structure_factor_val=uniform_structure_factor.getvalarray(norder_min,norder_max);
      uniform_structure_factor_val/=nmeasurementsval;
      measurements["Uniform Structure Factor Coefficients"] << uniform_structure_factor_val;
      if (is_bipartite()) {
        valarray<double> staggered_structure_factor_val=staggered_structure_factor.getvalarray(norder_min,norder_max);
        staggered_structure_factor_val/=nmeasurementsval;
        measurements["Staggered Structure Factor Coefficients"] << staggered_structure_factor_val;
      }
    }
  }
  else {
    int posmin=parms.value_or_default("START_STORING",logf_steps_total);
    if (posmin>logf_steps_total)
      posmin=logf_steps_total;
    if (atoi(pos.c_str())>=posmin) {
      valarray<double> gval=g.getvalarray(norder_min,norder_max);
      double red=gval[0]-num_sites()*log(2.);
      gval-=red;
      measurements["Coefficients "+pos] << gval;
      measurements["Total Sweeps "+pos] << (double)sweeps;
    }
  }
}

void QWL_SSE_Simulation::print_copyright(std::ostream& out)
{
  out << "Quantum Wang-Landau stochastic series expansion program\n"
      << "  copyright(c) 2004 Stefan Wessel <wessel@comp-phys.org>\n"
      << "  for details see the publications:\n"     
      << "  M. Troyer, S. Wessel, and F. Alet, Phys. Rev. Lett. 90, 120201 (2003).\n"
      << "  A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).\n\n";
}

#endif
