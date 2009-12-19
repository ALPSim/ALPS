/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2009 by Matthias Troyer <troyer@comp-phys.org>,
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

#include "vertex.h"
#include <alps/scheduler/montecarlo.h>
#include <alps/lattice/parity.h>
#include <boost/property_map/vector_property_map.hpp>
#include <vector>

//creates lots of output to trace configuation changes
//#define TEST_MODE

using namespace std;

//=======================================================================
// SSERun
//-----------------------------------------------------------------------

class SSERun : public alps::scheduler::LatticeMCRun<>
{
public:
  static void print_copyright(std::ostream&);
  
  SSERun(const alps::ProcessList&,const alps::Parameters&,int);
  SSERun(const alps::ProcessList&,alps::IDump&,int);

  void save(alps::ODump&) const;
  void load(alps::IDump&);
  
  void dostep();
  bool is_thermalized() const;
    
  void start();
  
private:
  typedef VertexString::iterator vertex_iterator;

  void setup_model();
  double diagonal_matrix_element(uint8_t l, uint8_t r, int) const { return diag_me[l][r]; }
  double offset() const { return offset_; }
  SSEPath breakup(int ,const LocalConfig c, double r) const {return r < probabilities_[c] ? path_[c][0] : path_[c][1];}
  friend class ::Vertex;
  
  double work_done() const;
  bool is_self_consistent() const; //to be used in assert(is_self_consistent())
                                  //verifies that all states are set correctly
                                 //in the vertex string true is good as the name suggests
  void do_diagonal();
  void do_web();
  void do_loops();
  void do_measure();
  inline bool next_vertex(vertex_iterator&, int& leg, Direction dir, bool flip, uint32_t start=0);
  
  double beta;
  uint32_t sweeps;
  uint32_t thermalization;
  StateVector initial_state;
  VertexString vertex_string;
  std::valarray<double> windings;
  bool issigned;
  double sign;
  double offset_;
  double diag_me[2][2];
  //probabilities for the different breakups
  double probabilities_[MAX_LOCAL_CONFIGS];
  SSEPath path_[MAX_LOCAL_CONFIGS][2];
};


bool SSERun::next_vertex(vertex_iterator& here, int& leg, Direction dir, bool flip, uint32_t start)
{
  assert(here>=vertex_string.begin() && here<vertex_string.end());

  uint32_t next, newleg;
  uint32_t here_i = here-vertex_string.begin();
  bool passed;
  if(dir==UP) { 
    next=here->upper_vertex(leg);
    if (next <= here_i) {
      if (flip)
        initial_state.flip(here->site(leg));
      passed = (here_i <= start || next > start);
    }
    else
      passed = (here_i <= start && next > start);
    newleg = here->upper_leg(leg);
    assert(vertex_string[next].lower_vertex(newleg)==here-vertex_string.begin());
    assert(vertex_string[next].lower_leg(newleg)==leg);
  }
  else {
    next=here->lower_vertex(leg);
    if (next >= here_i) {
      if (flip)
        initial_state.flip(here->site(leg));
      passed = (here_i > start || next <= start);
    }  
    else
      passed = (here_i > start && next <= start);
    newleg = here->lower_leg(leg);
    assert(vertex_string[next].upper_vertex(newleg)==here-vertex_string.begin());
    assert(vertex_string[next].upper_leg(newleg)==leg);
  }
  here = vertex_string.begin()+next;
  leg=newleg;
  return passed;
}

