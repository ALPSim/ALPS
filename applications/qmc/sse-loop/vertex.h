/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2001-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Fabien Alet <alet@comp-phys.org>,
*                            Andreas Lange <alange@phys.ethz.ch>
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

#ifndef VERTEX_H
#define VERTEX_H

#include <alps/lattice/graphproperties.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/vector.h>
#include <iostream>
#include <stdexcept>

enum SSEPath { JUMP, TURN, FROZEN, STRAIGHT, MAX_NUMBER_OF_PATHES};
enum Direction { DOWN, UP };
enum LocalConfig { DIAG_FM, DIAG_AFM, NONDIAG, MAX_LOCAL_CONFIGS};


class Vertex;

typedef unsigned int bondtype_type;
typedef unsigned int site_number_type;
typedef std::vector<Vertex>::size_type vertex_number_type;

class StateVector : public std::vector<alps::uint8_t>
{
public:
  typedef std::vector<alps::uint8_t> super_type;
  typedef super_type::size_type size_type;
  StateVector(size_type s=0) : std::vector<alps::uint8_t>(s) {}
  void flip(site_number_type i) { (*this)[i]=1-(*this)[i];}
  
  void save(alps::ODump& od) const { od << static_cast<const super_type&>(*this); }
  void load(alps::IDump& id) { id >> static_cast<super_type&>(*this); }
};


std::ostream& operator<< (std::ostream& ost, const StateVector& stv);

namespace alps {
inline alps::ODump& operator<<(alps::ODump&, const Vertex&);
inline alps::IDump& operator>>(alps::IDump&, Vertex&);
}

class Vertex {
public:
  Vertex() : type_(0) {}
  // state
  void make_identity() { type_=0; Vertex::n_--;} // TODO: remove
  void make_diagonal(int, int) { if(is_identity()) Vertex::n_++; type_=1;} // TODO: remove
  bool is_identity() const { return type_==0;} // TODO: remove
  bool is_diagonal() const { return type_==1;} 
  bool is_nondiagonal() const { return type_==2;}
  void flip() { assert(!is_identity()); type_ = 3-type_;}

  template<class E, class G> void set_sites(const E& e, const G& g) { 
    sites_[0]=boost::source(e,g); 
    sites_[1]=boost::target(e,g);
    bond_number_=boost::get(alps::edge_index_t(),g,e);
  }
  
  site_number_type site(int lr) const { return sites_[lr];}
  unsigned int bond_number() const { return bond_number_; }
  
  void reset() { std::fill(visited_,visited_+4,false); }
  void visit(int leg, int dir) { visited_[leg+2*dir]=true;}
  bool visited(int leg, int dir) const { return visited_[leg+2*dir];}

  vertex_number_type lower_vertex(int lr) const { return down_[lr];}
  vertex_number_type upper_vertex(int lr) const { return up_[lr];}
  int lower_leg(int lr) const { return down_leg_[lr];}
  int upper_leg(int lr) const { return up_leg_[lr];}
  void set_lower(int lr, vertex_number_type vertex, int leg)
  { down_[lr]=vertex; down_leg_[lr]=leg;}
  void set_upper(int lr, vertex_number_type vertex, int leg)
  { up_[lr]=vertex; up_leg_[lr]=leg;}

  template<class MODEL> void make_breakup(const MODEL& model,int s1, int s2, double r,int bt)
  {
    LocalConfig lc;
    if(is_diagonal()) 
      lc =(s1==s2 ? DIAG_FM : DIAG_AFM);
    else if(is_nondiagonal())
      lc=NONDIAG;

    breakup_ = model.breakup(bt,lc,r);
  }

  SSEPath path() { return breakup_;}

  void pass (int& leg, Direction& dir, bool toflip)
  {
    if (path() != STRAIGHT) {
      leg = 1 - leg;  
      if(path() == TURN )
        dir = Direction(1 - dir); 
    }
    if(toflip)
      flip();
  }

  // dumping
  friend alps::ODump& alps::operator<<(alps::ODump&, const Vertex&);
  friend alps::IDump& alps::operator>>(alps::IDump&, Vertex&);
  static vertex_number_type n_; // TODO: not needed
  
private:
  alps::uint8_t type_; // 0 = identity, 1= diagonal 2=offdiagonal TODO: change, we obly need 0=diagonbal 1=offdiagonal
  site_number_type sites_[2]; 
  SSEPath breakup_;
  unsigned int bond_number_;
  bool visited_[4];            // flags for loop building. TODO: can combine all 4 into 4 bits of a single alps::uint8_t
  vertex_number_type down_[2]; // links for the web
  vertex_number_type up_[2];  // TODO: can combine down_ and 1 bit for down_leg_ into one integer, saves memory!
  int down_leg_[2];
  int up_leg_[2];
};

class VertexString : public std::vector<Vertex>
{
public:
  typedef std::vector<Vertex> super_type;
  typedef super_type::size_type size_type;
  typedef super_type::iterator iterator;
  typedef super_type::const_iterator const_iterator;
  VertexString(size_type s=0) : super_type(s) {} // TODO: don't need sizes
  
  size_type n() const { return Vertex::n_;} //TODO: don't need n
  void init(); 

  template <class RNG> void grow(size_type s, RNG& rnd); // TODO: not needed

  void save(alps::ODump& od) const { od << static_cast<const super_type&>(*this); }
  void load(alps::IDump& id) { id >> static_cast<super_type&>(*this); }
};


namespace alps {

inline alps::ODump& operator<<(alps::ODump& dump, const Vertex& v) 
{
  return dump << v.sites_[0] << v.sites_[1] << v.type_ << v.bond_number_;
}

inline alps::IDump& operator>>(alps::IDump& dump, Vertex& v)
{
  return dump >> v.sites_[0] >> v.sites_[1] >> v.type_ >> v.bond_number_;
}

}


std::ostream& operator<<(std::ostream& out, const Vertex& v);

template <class RNG>
void VertexString::grow(size_type s, RNG& rng) // TODO: not needed
{
  super_type tmp(super_type::size()+s);
  double p = double(s)/double(tmp.size());
  // fill in id randomly
  iterator l = tmp.begin();
  iterator r = super_type::begin();
  size_type to_insert = s;
  while (to_insert && r!=super_type::end())
  {
    if (rng() < p) { // insert identity
      ++l;
      --to_insert;
    }
    else { // copy
      if( l<tmp.begin() || l>=tmp.end())
        boost::throw_exception(std::runtime_error("Out of range in VertexString::grow (l)\n"));
      if( r<super_type::begin() || r>=super_type::end() ) 
        boost::throw_exception(std::runtime_error("Operator string exceeded maximum length (r)\n"));
      *l++ = *r++;
    }
  }
  std::copy(r,super_type::end(),l);
  super_type::swap(tmp);
}


#endif
