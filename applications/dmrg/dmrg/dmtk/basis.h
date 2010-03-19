/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
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

#ifndef __DMTK_BASIS_H__
#define __DMTK_BASIS_H__

#include <iostream>
#include <vector>
#include <iosfwd>
#include "vector.h"
#include "subspace.h"
#include "qn.h"
#include "lanczos.cc"

namespace dmtk
{

class PackedBasis: public std::vector<SubSpace>
{
  public:
    typedef std::vector<SubSpace> _V;
    typedef std::vector<SubSpace>::iterator iterator;
    typedef std::vector<SubSpace>::const_iterator const_iterator;
  
    PackedBasis() {}
    PackedBasis(const PackedBasis &b): _V(b) {}

    PackedBasis operator=(const PackedBasis &b) 
      { _V::operator=(b); return *this; }
  
    SubSpace operator()(size_t i) const
      {
        const_iterator iter;
        for(iter = begin(); iter != end(); iter++){
          if(i >= (*iter).begin() && i <= (*iter).end()) return *iter;
        }  
        return SubSpace(QN(),DMTK_ERROR,DMTK_ERROR);
      }

    SubSpace operator() (const QN &qn) const
      {
        int i = get_index(qn);
        if(i >= 0) return _V::operator[](i);

//        cout << "*** WARNING 2: Subspace not found\n";
        return SubSpace(QN(),DMTK_ERROR,DMTK_ERROR);
      }

    int get_index(const QN &qn) const
      {
// WARNING: MAY NEED MORE TESTING
/*
        PackedBasis::const_iterator iter;
        int i = 0;
        for(iter = begin(); iter != end(); iter++, i++) 
          if((*iter).qn() == qn) return i;
*/

        int origin = 0;
        int end = size() - 1;
        int index = -1;

        while(origin <= end){
          int index_old = index;
          index = (origin + end) / 2;

          const SubSpace &s = _V::operator[](index);

          if(s.qn() == qn) return index;
          if(s.qn() > qn)
            end = index;
          else
            origin = index;

          if(index == index_old){
            if(index == end)
              end = end - 1;
            else
              origin = origin + 1;
          }
        }

        return -1;
      }

    size_t dim() const
      {
         const_iterator iter;
         size_t _dim = 0;
         for(iter = begin(); iter != end(); iter++) _dim += (*iter).dim();
         return _dim;
      }

//  Streams

    void write(std::ostream &s) const
    {
      size_t _size = size();
      s.write((const char *)&_size, sizeof(size_t));

      const_iterator iter;
      for(iter = begin(); iter != end(); iter++) (*iter).write(s);
    }

    void read(std::istream &s)
    {
      size_t _size = 0;
      s.read((char *)&_size, sizeof(size_t));
      resize(_size);

      iterator iter;
      for(iter = begin(); iter != end(); iter++) (*iter).read(s);
    }

};


class State
{
  private:
    QN _qn1;
    QN _qn2;
  public:
    size_t i1;
    size_t i2;

    State():i1(0),i2(0) {}
    State(size_t _i1, QN qn1, size_t _i2 = 0, QN qn2 = QN()): _qn1(qn1),_qn2(qn2),i1(_i1),i2(_i2) {}
    State(const State &s): _qn1(s._qn1),_qn2(s._qn2),i1(s.i1),i2(s.i2) {}

    State& operator=(const State &s) 
      { i1 = s.i1; i2 = s.i2; _qn1 = s._qn1; _qn2 = s._qn2; return *this; }

    QN qn() const { return _qn1+_qn2; } 
    QN qn1() const { return _qn1; }
    QN qn2() const { return _qn2; }
    QN& qn1() { return _qn1; }
    QN& qn2() { return _qn2; }

    bool operator==(const State& s) const { QN qaux = _qn1+_qn2; return (qaux == s.qn()); }
    bool operator!=(const State& s) const { QN qaux = _qn1+_qn2; return (qaux != s.qn()); }
    bool operator>(const State& s) const { QN qaux = _qn1+_qn2; return (qaux > s.qn()); }
    bool operator>=(const State& s) const { QN qaux = _qn1+_qn2; return (qaux >= s.qn()); }
    bool operator<(const State& s) const {  return (qn() < s.qn()); }
    bool operator<=(const State& s) const { QN qaux = _qn1+_qn2; return (qaux <= s.qn()); }

    // Streams

    void write(std::ostream &s) const
    {
      _qn1.write(s);
      _qn2.write(s);
      s.write((const char *)&i1, sizeof(size_t));
      s.write((const char *)&i2, sizeof(size_t));
    }

    void read(std::istream &s)
    {
      _qn1.read(s);
      _qn2.read(s);
      s.read((char *)&i1, sizeof(size_t));
      s.read((char *)&i2, sizeof(size_t));
    }
};


class Basis: public Vector<State>
{
  private:
    Basis& init(const Basis&, const Basis &);
    Basis& init(const PackedBasis&, const PackedBasis &);
    PackedBasis _subspace;
    bool ordered;
  public:
    typedef Vector<State> _V;
    typedef Vector<State>::iterator iterator;
    typedef Vector<State>::const_iterator const_iterator;

    Basis(): Vector<State>(), ordered(false) {}
    Basis(size_t n): Vector<State>(n), ordered(false) 
       { for(int i = 0; i < n; i++) operator[](i) = State(i,QN()); }
    Basis(const Vector<State> &v): Vector<State>(v), ordered(false) {}
    Basis(const Basis &b): Vector<State>(b),_subspace(b._subspace), ordered(b.ordered) {}
    Basis(const Basis &b1, const Basis &b2): ordered(false) { init(b1,b2); }    
    Basis(const PackedBasis &b1, const PackedBasis &b2): ordered(false) { init(b1,b2); }    

    Basis& resize(const Basis &b1, const Basis &b2) 
      { ordered = false; init(b1,b2); return *this; }    

    Basis& resize(size_t n) 
      { 
        Vector<State>::resize(n); 
        if(ordered) { ordered = false; reorder(); }
        return *this; 
      }

    Basis& operator=(const Vector<State>& v)
      { Vector<State>::operator=(v); return *this; }
    Basis& operator=(const Basis& b)
      { Vector<State>::operator=(b); _subspace = b._subspace; ordered = b.ordered; return *this; }  

    Basis& pack();
    Basis& reorder(bool ascending = true);


    PackedBasis::iterator subspace_begin() { return _subspace.begin(); }
    PackedBasis::const_iterator subspace_begin() const { return _subspace.begin(); }
    PackedBasis::iterator subspace_end() { return _subspace.end(); }
    PackedBasis::const_iterator subspace_end() const { return _subspace.end(); }

    const PackedBasis& subspaces() const { return _subspace; }

    State& operator() (size_t n) { return Vector<State>::operator()(n); }
    State operator() (size_t n) const { return Vector<State>::operator()(n); }
    State& operator[] (size_t n) { return Vector<State>::operator[](n); }
    State operator[] (size_t n) const { return Vector<State>::operator[](n); }

    SubSpace operator() (const QN &qn) const
      { return _subspace(qn); }

    size_t dim() const { return size(); }

    // Streams

    void write(std::ostream &s) const
    {
      s.write((const char *)&ordered, sizeof(bool));

      size_t l = size();
      s.write((const char *)&l, sizeof(size_t));

      Vector<State>::const_iterator iter;
      for(iter = begin(); iter != end(); iter++)
        (*iter).write(s);
    }

    void read(std::istream &s)
    {
      s.read((char *)&ordered, sizeof(bool));

      size_t l;
      s.read((char *)&l, sizeof(size_t));
      Vector<State>::resize(l);

      Vector<State>::iterator iter;
      for(iter = begin(); iter != end(); iter++){
        (*iter).read(s);
      }

      if(ordered){ ordered = false; reorder(); }
    }
};

Basis &
Basis::init(const Basis &b1, const Basis &b2)
{
  int dim = b1.size()*b2.size();
  resize(dim);
  int i = 0;
  for(int i1 = 0; i1 < b1.size(); i1++){
    for(int i2 = 0; i2 < b2.size(); i2++){
      Vector<State>::operator[](i) = State(i1,b1[i1].qn(),i2,b2[i2].qn()); 
      i++;
    }
  }

  return *this;
}

Basis &
Basis::init(const PackedBasis &b1, const PackedBasis &b2)
{
  int dim = b1.dim()*b2.dim();
  resize(dim);
  int i = 0;
  PackedBasis::const_iterator iter1;
  PackedBasis::const_iterator iter2;
  for(iter1 = b1.begin(); iter1 != b1.end(); iter1++)
    for(iter2 = b2.begin(); iter2 != b2.end(); iter2++)
      for(int i1 = (*iter1).begin(); i1 <= (*iter1).end(); i1++)
        for(int i2 = (*iter2).begin(); i2 <= (*iter2).end(); i2++){
          Vector<State>::operator[](i) = State(i1,(*iter1).qn(),i2,(*iter2).qn()); 
          i++;
        }

  return *this;
}

Basis& 
Basis::pack()
{
  _subspace.clear();
  SubSpace s(operator[](0).qn(),0,0); 
  for(int i = 1; i < size(); i++) {
    if(s.qn() != operator[](i).qn()){
      _subspace.push_back(SubSpace(s.qn(),s.begin(),i-1));
      s = SubSpace(operator[](i).qn(),i,i);
    }
  }
  _subspace.push_back(SubSpace(s.qn(),s.begin(),size()-1));

  return *this;
}

Basis&
Basis::reorder(bool ascending)
{
  if(ordered) return *this;
  ordered = true;

  sort<State>(size(), *this, ascending);

  // Calculate subspaces

  pack();

  return *this;
}



} //namespace dmtk

#endif // __DMTK_BASIS_H__
