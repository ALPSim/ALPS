/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2010 by Lode Pollet <pollet@itp.phys.ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id: lowa.hpp,v 1.1 2006/09/09 9:21:44 pollet Exp $ */

#ifndef lowa_Element_HPP
#define lowa_Element_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include <list>

using namespace std;
const int dim=3;
const int zcoord=2*dim;
typedef int32_t SiteType;                  // do not make unsigned!!!
typedef double TimeType;


class Element
{
 public:
  Element() {};
  Element(const int n0, const int n1, const SiteType s0, const SiteType s1, const TimeType t, const list<Element>::iterator v[zcoord] ) {
    mTime = t; 
    mBefore = n0; 
    mAfter = n1;
    mFrom = s0;
    mTo = s1;
    for (int i = 0; i < zcoord; i++) mAssoc[i] = v[i];
  }
  Element(const int n0, const int n1, const SiteType s0, const SiteType s1, const TimeType t ) {
    mTime = t; 
    mBefore = n0; 
    mAfter = n1;
    mFrom = s0;
    mTo = s1;
  }
  void init(const int n0, const int n1, const SiteType s0, const SiteType s1, const TimeType t) {
    time(t); 
    before(n0); 
    after(n1);
    from(s0);
    to(s1);
  }
  void init(const int n0, const int n1, const SiteType s0, const SiteType s1, const TimeType t, const list<Element>::iterator  v[zcoord]) {
    time(t); 
    before(n0); 
    after(n1);
    assoc(v);
    from(s0);
    to(s1);
  }
  ~Element() {}
  Element(const Element& src);
  Element& operator=(const Element& rhs);
  friend std::ostream &operator<<( std::ostream &os, const Element& rhs) {
    os << rhs.time() << "\t" 
       << rhs.from() << "\t"
       << rhs.to() << "\t"
       << rhs.before() << "\t" 
       << rhs.after();
    return (os);
  }
  friend std::istream& operator>>(std::istream& is, Element& e) {
    is >> e.mTime >> e.mFrom >> e.mTo >> e.mBefore >> e.mAfter;
    return (is);
  }
  friend bool operator== (const Element& lhs, const Element& rhs) {
    bool is_equal;
    ((lhs.mFrom == rhs.mFrom) && (lhs.to() == rhs.to()) && (lhs.time() == rhs.time()) && (lhs.before() == rhs.before()) && (lhs.after() == rhs.after()) ) ? is_equal=true : is_equal=false;
    return (is_equal);
  }    


  TimeType time() const {return (mTime);}
  int before() const {return (mBefore);}
  int after() const {return (mAfter);}
  int name() const { return (mName);}
  SiteType from() const {return mFrom;}
  SiteType to() const { return mTo;}

  void time(const TimeType t) {mTime=t;};
  void before(const int t) {mBefore = t;};
  void after(const int t) {mAfter = t;};
  void from(const SiteType s) {mFrom = s;}
  void to(const SiteType s) { mTo = s;}
  void name(const int n) {mName = n;}
  void assoc(const list<Element>::iterator  v[zcoord]) {
    for (int i = 0; i < zcoord; i++) mAssoc[i] = v[i];
  }

  void print() {
    cout << "\n" << this << "\t" << mTime << "\t" << mFrom << "\t" << mTo << "\t" << mBefore << "\t" << mAfter;
    for (SiteType j = 0; j < zcoord; j++) cout << "\t" << &(*mAssoc[j]) << "\t" << mAssoc[j]->time();
  }
  void get_assoc(list<Element>::iterator* v) {
    for (SiteType i = 0; i < zcoord; i++) { 
      v[i] = mAssoc[i]; 
      //cout << "\n Getasso : " << &(*mAssoc[i]);
    }
  }
  list<Element>::iterator& get_assoc(SiteType s) { return mAssoc[s];}
  void set_assoc(const SiteType j, const list<Element>::iterator v) {
    mAssoc[j] = v;
  }


 private:
  TimeType mTime;  // time of the interaction
  int mBefore;     // occupation on the site before the interaction
  int mAfter;      // occupation on the site after the interaction
  SiteType mFrom;  // site from which the particle hops
  SiteType mTo;    // site to which the particle hops
  list<Element>::iterator mAssoc[zcoord];   // associations; i.e. iterator to the kink on nb sites that are equal or just greater in time
  int mName;       // numerical name of an element; used to save/load associations of elements to/from file
 };

inline Element::Element(const Element& src) : mTime(src.time()), mBefore(src.before()), mAfter(src.after()), mFrom(src.from()), mTo(src.to()), mName(-1)
{
  for (int i = 0; i < zcoord; i++) mAssoc[i] = src.mAssoc[i];
}

inline Element& Element::operator=(const Element& rhs)
{
  //std::cout << "\n assigning Element\n";
  if (this == &rhs)
  {
    return (*this);
  }
  mTime = rhs.mTime;
  mBefore = rhs.mBefore;
  mAfter = rhs.mAfter;
  mFrom = rhs.mFrom;
  mTo = rhs.mTo;
  for (int i = 0; i < zcoord; i++) mAssoc[i] = rhs.mAssoc[i];
}



template <typename S, typename T>
class Matrix
{
private:
  T* mCell;
  S mRow, mCol, mZ, mSize;
public:
  Matrix() {}
  void init(S x, S y) {mRow = x; mCol = y;mZ = 1; mSize = mRow*mCol*mZ; mCell = new T [mRow*mCol];}
  void init(S x, S y, S z) {mRow = x; mCol = y;mZ = z; mSize = mRow*mCol*mZ; mCell = new T [mRow*mCol*mZ];}
  Matrix(S x, S y, S z) {mRow = x; mCol = y; mZ = z; mSize = mRow*mCol*mZ; mCell = new T [mRow*mCol*mZ];}
  Matrix(S x, S y) {mRow = x; mCol = y; mZ = 1; mSize = mRow*mCol*mZ; mCell = new T [mRow*mCol];}

  Matrix(const Matrix& rhs) {
    mRow = rhs.mRow;
    mCol = rhs.mCol;
    mZ   = rhs.mZ;
    mCell = new T [mRow*mCol*mZ];
    for (S i = 0; i < mSize; i++) {
      mCell[i] = rhs.mCell[i];
    }
  }

  Matrix& operator=(const Matrix& rhs) {
    // just copy elements, no delete-new
    if (rhs.mSize != mSize) {
      std::cout << "\nError in Matrix assignment "  << rhs.mSize << "\t" << mSize;
    }
    for (S i = 0; i < mSize; i++) {
      mCell[i] = rhs.mCell[i];
    }
  }


  S size() const {return mSize;}
  S dim1() const {return mRow;}
  S dim2() const {return mCol;}
  S dim3() const {return mZ;}

  T& operator() (const S x, const S y) const {return mCell[x*mCol +y];}    //CHANGED TO FORTRAN STYLE!!!
  //T& operator() (const S x, const S y, const S z) const {return mCell[(z*mCol +y)*mRow +x];}
};



#endif
