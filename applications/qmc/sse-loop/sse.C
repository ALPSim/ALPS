/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabien Alet <alet@comp-phys.org>
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

#include "sse.h"
#include "vertex.h"
#include <alps/osiris/os.h>
#include <alps/alea.h>
#include <algorithm>
#include <alps/osiris/comm.h>
#include <boost/ref.hpp>

#include <list>
#include <cmath>
//#include <stdlib.h>

using namespace std;
using namespace alps;
//#define TIMINGS

template <class T> T sqr(T x) { return x*x;}

//=======================================================================
// SSERun
//-----------------------------------------------------------------------

void SSERun::print_copyright(std::ostream& out)
{
  out << "Quantum Monte simulations using the loop algorithm (SSE version)\n"
      << "  available from http://alps.comp-phys.org/\n"
      << "  copyright(c) 1999-2008 by Fabien Alet <alet@comp-phys.org>,\n"
      << "                            Synge Todo <wistaria@comp-phys.org>,\n"
      << "                            Matthias Troyer <troyer@comp-phys.org>\n"
      << "                            and Jean=David Picon <picon@comp-phys.org>\n"
      << " for details see the publication:\n"     
      << " A.F. Albuquerque et al., J. of Magn. and Magn. Materials 310, 1187 (2007).\n\n";
}

SSERun::SSERun(const ProcessList& where,const alps::Parameters& p,int node)
  : alps::scheduler::LatticeMCRun<>(where,p,node),
    beta(parms.defined("BETA") ? double(parms["BETA"]) : 1./double(parms["T"])),
    sweeps(0),
    thermalization(int(parms["THERMALIZATION"])),
    initial_state(num_sites()),
    vertex_string(int(parms.defined("CUTOFF") ? int(parms["CUTOFF"]) : 0.5 * num_sites() * beta)), // TODO: variable length starts with 0
    windings(dimension()),
    issigned(!is_bipartite() && double(parms.value_or_default("Jxy",double(parms.value_or_default("J",1.))))>0. ) // TODO: later have more than one J
    // issigned(has_sign_problem()); // TODO: use this once we actually use the model library
{
  if (inhomogeneous())
    boost::throw_exception(std::runtime_error("Disordered lattices not supported by the simple SSE loop program.\n"));

  if (issigned)
    measurements << RealObservable("Sign");
    
  // since "Sign" is the default name for the sign observable, we do  not need to 
  // to pass the name to alps::make_observable
  
  measurements << alps::make_observable(RealObservable("Energy"),issigned);
  measurements << alps::make_observable(RealObservable("n"),issigned);
  measurements << alps::make_observable(RealObservable("n^2"),issigned);
  measurements << alps::make_observable(RealObservable("Magnetization^2"),issigned);
  measurements << alps::make_observable(RealObservable("Magnetization^4"),issigned);
  measurements << alps::make_observable(RealObservable("Susceptibility"),issigned);

  if (is_bipartite()) {
    measurements << alps::make_observable(RealObservable("Staggered Magnetization^2"),issigned);
    measurements << alps::make_observable(RealObservable("Staggered Magnetization^4"),issigned);
  }
  
  measurements << alps::make_observable(RealVectorObservable("Windings"),issigned);
  measurements << alps::make_observable(RealObservable("Spin stiffness"),issigned);
  measurements << alps::make_observable(RealObservable("Density"),issigned);
  
  // these are update statistics and not physical observables.
  // they thus do not have a sign problem
  measurements << SimpleIntObservable("loops");
  measurements << SimpleRealObservable("loop length");

  setup_model();
}

void SSERun::load(IDump& dump)
{
 alps::scheduler::MCRun::load(dump);
 dump >> sweeps;
  if(!where.empty()) {
    dump  >> initial_state;
    dump  >> vertex_string;
  }
}

void SSERun::start()
{
  vertex_string.init();
}

void SSERun::save(ODump& dump) const
{
  alps::scheduler::MCRun::save(dump);
  dump << sweeps << initial_state << vertex_string; 
}

bool SSERun::is_thermalized() const
{
  return (sweeps >= thermalization); // thermalized
}

void SSERun::dostep()
{
  //  if (!(sweeps%1000)) { cout << sweeps << endl;}

#ifdef TIMINGS
  double t = -dclock();
#endif
  sweeps++;
  do_diagonal();

#ifdef TIMINGS
  t += dclock();
  cerr << "Diagonal time: " << t << "\n";
  t = - dclock();
#endif  
    do_measure();
#ifdef TIMINGS
  t += dclock();
  cerr << "Measure time: " << t << "\n";
  t = - dclock();
#endif 
  do_web();
#ifdef TIMINGS
  t += dclock();
  cerr << "Web time: " << t << "\n";
  t = - dclock();
#endif 
  do_loops();
  assert(is_self_consistent());

#ifdef TIMINGS
  t += dclock();
  cerr << "Loop time: " << t << "\n";
  t = - dclock();
#endif
  
}

void SSERun::do_diagonal()
{
    for (int i=0;i<windings.size();++i)
      windings=0.;
    
    StateVector state(initial_state); // this gets propagated

    VertexString new_string; // TODO: for variable length
    sign=1.;
    for (vertex_iterator it=vertex_string.begin(); it !=vertex_string.end() ; ++it)
    {
      if(it->is_identity()) { // TODO: replace by code for insertions
        // id -> diag
        bond_descriptor bond=*(bonds().first+random_int(0,num_bonds()-1)); // find a random bond
        
        //p like probability                                                        //this is the bond-type (bt)
        double p=diagonal_matrix_element(state[source(bond)],state[target(bond)],bond_type(bond))
             *num_bonds()*beta / (vertex_string.size()-vertex_string.n());
        if (p>=1. || random_real()<p)  {
          it->make_diagonal(state[source(bond)],state[target(bond)]);
          it->set_sites(bond,graph());
        }
      }
      else if (it->is_diagonal()) { // TODO: replace by new rmoving/copying 
        double p = double(vertex_string.size()-long(vertex_string.n())+1) /
                         (diagonal_matrix_element(state[it->site(0)],state[it->site(1)],bond_type(bond(it->bond_number())))
                           *num_bonds()*beta);
        if (p >= 1. || random_real()<p)
          it->make_identity();
      }
      else if (it->is_nondiagonal()) { // TODO: replace by copying
        vector_type v=bond_vector_relative(bond(it->bond_number()));
        double spin = (state[it->site(0)] ? 1 : -1);
        for (int i=0;i<v.size() && i<dimension();++i)
          windings[i]+=spin*v[i];
        assert(state[it->site(0)]!=state[it->site(1)]);
        state.flip(it->site(0));
        state.flip(it->site(1));
        if (issigned) //depend on bond type in general
          sign =-sign;
      }
      else
        boost::throw_exception(std::logic_error("Illegal vertex type"));
      
      // make the breakup, if necessary
      if (!it->is_identity())
        it->make_breakup(*this,state[it->site(0)],state[it->site(1)],random_real(),bond_type(bond(it->bond_number())));
    }
  
    if (!is_thermalized() && vertex_string.n() > vertex_string.size()*0.8)  { // TODO: remove
      vertex_string.grow(100+int(6.*std::sqrt(double(vertex_string.size()))),random_01);
    }
    if (vertex_string.size()==vertex_string.n()) // TODO: remove
      boost::throw_exception(std::runtime_error("Operator string exceeded maximum length\n"));
    assert((state==initial_state));
    // TODO: vertex_string.swap(new_string);
}

  // build vertex web
void SSERun::do_web() // TODO: nothing, unchanged, build the links
{  
  vector<alps::uint32_t> current_vertex(num_sites(),vertex_string.size());
  vector<alps::uint32_t> first_vertex(num_sites(),vertex_string.size());
  vertex_iterator it;
  alps::uint32_t i=0;
  for (it=vertex_string.begin(); it !=vertex_string.end() ; ++it, ++i)  {
    if(!it->is_identity())  {
      for (int leg=0; leg<2;++leg)  {
        alps::uint32_t site=it->site(leg);
        if (current_vertex[site]==vertex_string.size())
          first_vertex[site]=i; // I'm the first vertex on this site
        else  {
          // we can make links
          vertex_iterator prev = vertex_string.begin()+current_vertex[site];
          
          if(prev<vertex_string.begin() || prev>=vertex_string.end()) 
            boost::throw_exception(std::runtime_error("Out of range in do_web\n"));
          
          int prev_leg = ( (prev->site(0)==site) ? 0 : 1);
          assert(prev->site(prev_leg)==site);
          prev->set_upper(prev_leg,i,leg);
          it->set_lower(leg,current_vertex[site],prev_leg);
        }
        current_vertex[site]=i;
      }
    }
  }
 
  // fix ends around time-boundary
  for (alps::uint32_t site = 0 ; site < num_sites() ; ++ site)
    if (current_vertex[site]!=vertex_string.size())  {
      // was there at least one vertex connected to this site?
      // connect it around boundary
      assert(first_vertex[site]!=vertex_string.size());
      int first_leg = ((vertex_string[first_vertex[site]].site(0)==site) ? 0 : 1);
      int last_leg = ((vertex_string[current_vertex[site]].site(0)==site) ? 0 : 1);
      assert(vertex_string[first_vertex[site]].site(first_leg)==site);
      assert(vertex_string[current_vertex[site]].site(last_leg)==site);
      vertex_string[first_vertex[site]].set_lower(first_leg,current_vertex[site],last_leg);
      vertex_string[current_vertex[site]].set_upper(last_leg,first_vertex[site],first_leg);
    }
    else // free spin flips for random spins
      initial_state[site]=random_int(2);
}


// loop updates
void SSERun::do_loops()
{ 
  #ifdef TEST_MODE
  std::cout << "Loop Update Step:\n";
  diagnostics();
  #endif //TEST_MODE

  for_each(vertex_string.begin(),vertex_string.end(),mem_fun_ref(&Vertex::reset));
  int loops=0;
  vertex_iterator start = vertex_string.begin();
  alps::uint32_t loop_length(0);
  
  do {
    if (!start->is_identity() ) { 
      int start_leg;
      Direction start_dir;
      bool found=false;
      vertex_iterator here=start;
      bool flip=random_real()<0.5; // TODO for WIESE: we flip by leaving the last spin configuration after all local cluster updates
      for (int leg=0; !found && leg<2 ;++leg)
        for (Direction dir=UP; !found && dir>= DOWN; dir=Direction(int(dir)-1)) {
          if(!here->visited(leg,dir)){
            start_leg=leg;
            start_dir=dir;
            found = true;
          }
        }
      if(!found){      
        ++start;
        continue;
      }

      if(here<vertex_string.begin() || here>=vertex_string.end()) 
        boost::throw_exception(std::runtime_error("Out of range in ::do_loops\n"));
      if(start<vertex_string.begin() || start>=vertex_string.end()) 
        boost::throw_exception(std::runtime_error("Out of range in ::do_loops\n"));
      
      vertex_iterator next=here;//only for loop condition
      int leg=start_leg;
      Direction dir=start_dir;
      ++loops;
      do {
        loop_length++;
        assert(!here->visited(leg,dir));
        here->visit(leg,dir);
        here->pass(leg,dir,flip); 
        assert(!here->visited(leg,Direction(1-dir)));
        here->visit(leg,Direction(1-dir));
        next_vertex(here,leg, dir, flip);
      } while (leg!=start_leg || here !=next || dir !=start_dir);
      
      // TODO WIESE: keep track of internal vertices, then sample this loop cluster!
    }
    else{//start->is_identity()==true // TODO: remove
      ++start;
      continue;
    }
  } while(start != vertex_string.end());
  measurements["loops"] << loops;
  measurements["loop length"] << (loops ? double(loop_length)/double(loops) : 0.);
}  


// perform measurements
void SSERun::do_measure()
{
  assert(count_if(vertex_string.begin(),vertex_string.end(),not1(std::mem_fun_ref(&Vertex::is_identity)))==vertex_string.n());
   
  double staggered_magnetization=0;
  double magnetization=0;
  for (site_iterator it=sites().first;it!=sites().second;++it) {
    double m = 0.5*(initial_state[*it] ? 1 : -1);
    magnetization += m;
    staggered_magnetization+=parity(*it)*m;
  }
  
  double winding=0.;
  for (int i=0;i<windings.size();++i) {
    windings[i]=windings[i]*windings[i]/beta;
    winding += windings[i];
  }
  
  // normalize measurements and add them to the observables
  if (issigned)
    measurements["Sign"] << sign;
  else if (sign<0.)
    boost::throw_exception(std::logic_error("Sign problem encountered unexpectedly"));
  measurements["Energy"] << sign*(-double(vertex_string.n())/beta
                     +offset()*num_bonds());
  measurements["n"] << sign*vertex_string.n();
  measurements["n^2"] << sign*vertex_string.n()*vertex_string.n();
  windings *=sign;
  measurements["Windings"] << windings;
  measurements["Spin stiffness"] << sign*winding/dimension();
    

  measurements["Magnetization^2"] << sign*sqr(magnetization/num_sites());
  measurements["Susceptibility"] << sign*sqr(magnetization)/num_sites()*beta;
  measurements["Magnetization^4"] << sign*sqr(sqr(magnetization/num_sites()));
  if (is_bipartite()) {
    measurements["Staggered Magnetization^2"] << sign*sqr(0.5*staggered_magnetization/num_sites());
    measurements["Staggered Magnetization^4"] << sign*sqr(sqr(0.5*staggered_magnetization/num_sites()));
  }
}


bool SSERun::is_self_consistent() const
{
  int i;
  StateVector state(initial_state);
  assert((state == initial_state));
  for (i=0; i<vertex_string.size(); i++)
    if (vertex_string[i].is_nondiagonal())  {
      state.flip(vertex_string[i].site(0));
      state.flip(vertex_string[i].site(1));
    }
  return ( initial_state == state );
}

double SSERun::work_done() const
{
  return (is_thermalized() ? (sweeps-thermalization)/double(parms["SWEEPS"]) :0.);
}


void SSERun::setup_model() // TODO later: we want Heisenberg models with different bond types
{ 
  double J=parms.value_or_default("J",1.);
  
  double epsilon=0.;

  for(int i=0;i<MAX_LOCAL_CONFIGS;i++) {
    probabilities_[i]=1.;
    path_[i][0]=MAX_NUMBER_OF_PATHES;
    path_[i][1]=MAX_NUMBER_OF_PATHES;
  }
    
  if(J< 0.) {
    path_[NONDIAG][0]=JUMP;
    path_[DIAG_FM][0]=JUMP;
  }
  // HBAFM
  else  {
    epsilon = parms.value_or_default("EPSILON", 0.);
    path_[NONDIAG][0]=TURN;
    path_[DIAG_FM][0]=STRAIGHT;
    path_[DIAG_AFM][0]=TURN;
    path_[DIAG_AFM][1]=STRAIGHT;
    probabilities_[DIAG_AFM]=1./(1.+epsilon);
  }

  offset_ = std::abs(J)/4.+epsilon*J/2;
  diag_me[0][0]=diag_me[1][1]= -J/4. + offset_;
  diag_me[1][0]=diag_me[0][1]= +J/4. + offset_;
}


