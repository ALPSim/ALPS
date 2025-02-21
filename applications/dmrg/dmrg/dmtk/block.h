/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
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

#ifndef __DMTK_BLOCK_H__
#define __DMTK_BLOCK_H__

#include <iosfwd>
#include <list>
#include "enums.h"
#include "vector.h"
#include "matrix.h"
#include "qn.h"
#include "state.h"
#include "subspace.h"
#include "basis.h"
#include "block_matrix.h"
#include "operators.h"
#include "lattice.h"

namespace dmtk
{

template<class T> class System;

template <class T>  
class Block: public std::list<BasicOp<T> >
{
  private:
    void resize_basis(const Basis &b);
    Basis _basis;
    Lattice _lattice;
    size_t _n_orbitals; // number of orbitals or internal sites;
    bool _single_site; // single site block ?
  public:
    typedef typename std::list<BasicOp<T> > _V;
    typedef typename std::list<BasicOp<T> >::iterator iterator;
    typedef typename std::list<BasicOp<T> >::const_iterator const_iterator;

//  Methods
    Block(size_t norbitals = 1): _n_orbitals(norbitals), _single_site(true) {};  
    Block(const Basis &basis, size_t norbitals = 1): _n_orbitals(norbitals), _single_site(true) { resize_basis(basis); }
    Block(const Basis &b1, const Basis &b2, size_t norbitals = 1): _n_orbitals(norbitals), _single_site(false) { resize_basis(Basis(b1,b2)); }
    Block(const Block<T> &b1, const Block<T> &b2): _n_orbitals(b1._n_orbitals+b2._n_orbitals), _single_site(false) 
      { resize_basis(Basis(b1.basis,b2.basis)); } 
    Block(const Block<T> &b): _V((const _V&)b), _basis(b._basis), _lattice(b._lattice), _n_orbitals(b._n_orbitals), _single_site(b._single_site) {}

    Block<T>& operator=(const Block<T>& b)
      { 
        _V::operator=((const _V&)b); 
        _basis = b._basis; 
        _lattice = b._lattice; 
        _n_orbitals = b._n_orbitals;
        _single_site = b._single_site;
        return *this; 
      }

    Block<T>& operator=(const T& v)
      {
        iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++){
          BMatrix<T> &op = (*iter);
          op = v;
        }
        return *this;
      }

    bool contains(const char *, int site, int internal_site = -1, size_t type = OP_SYSTEM) const;
    bool contains(const string&, int site, int internal_site = -1, size_t type = OP_SYSTEM) const;
    bool contains(const BasicOp<T> &op) const;
    bool contains(const BasicOp<T> *op) const;

    const BasicOp<T>* operator()(const char *, int site, int internal_site = 0, size_t type = OP_SYSTEM) const;
    BasicOp<T>* operator()(const char *, int site, int internal_site = 0, size_t type = OP_SYSTEM);
    const BasicOp<T>* operator()(const BasicOp<T> &op) const;
    BasicOp<T>* operator()(const BasicOp<T> &op);

    Block& reflect() { return *this; }

    Block& resize(const Block<T> &b1, const Block<T>& b2)
      { resize_basis(Basis(b1.basis,b2.basis)); _n_orbitals = b1._n_orbitals + b2._n_orbitals; return *this; }
    Block& resize(const Basis &b1, const Basis& b2)
      { resize_basis(Basis(b1,b2)); return *this; }
    Block& resize(const Basis &b)
      { resize_basis(b); return *this; }

    size_t dim() const { return _basis.dim(); }
    const Basis& basis() const { return _basis; }
    Basis& basis() { return _basis; }

    const Lattice& lattice() const { return _lattice; }
    Block& set_lattice(const Lattice& lattice) { _lattice = lattice; return *this; }
    size_t n_orbitals() const { return _n_orbitals; }
    Block& set_orbitals(size_t n) { _n_orbitals = n; return *this; }
    bool single_site() const { return _single_site; }
    Block& set_single_site(bool yn) { _single_site = yn; return *this; }

    friend class System<T>;
    friend class VectorState<T>;

    // Streams

    void write(std::ostream &s) const
    {
      size_t l = _V::size();
      s.write((const char *)&l, sizeof(size_t));
      s.write((const char *)&_n_orbitals, sizeof(size_t));
      s.write((const char *)&_single_site, sizeof(bool));

      _basis.write(s);
      _lattice.write(s);

      const_iterator iter;
      for(iter = _V::begin(); iter != _V::end(); iter++){
        (*iter).write(s);
      }
    }

    void read(std::istream &s)
    {
      size_t l;
      s.read((char *)&l, sizeof(size_t));
      s.read((char *)&_n_orbitals, sizeof(size_t));
      s.read((char *)&_single_site, sizeof(bool));
      list<BasicOp<T> >::resize(l);

      _basis.read(s);
      _lattice.read(s);

      iterator iter;
      for(iter = _V::begin(); iter != _V::end(); iter++) {
        (*iter).read(s);
      }
    }

};

template <class T>
inline BasicOp<T>*  
Block<T>::operator()(const char *name, int site, int internal_site, size_t type)
{
  iterator iter;
  for(iter = _V::begin(); iter != _V::end(); iter++){
    if((*iter).name() == string(name) && (*iter).site() == site &&
       (*iter).type() == type && (*iter).internal_site() == internal_site) {
      return &(*iter);
    }
  }
#ifdef WITH_WARNINGS
  cout << "*** WARNING: Operator 1: " << name << "(" << site << "," << internal_site << ") not found\n";
#endif // WITH_WARNINGS
  return NULL; 
}

template <class T>
inline const BasicOp<T>*
Block<T>::operator()(const char *name, int site, int internal_site, size_t type) const
{
  const_iterator iter;
  for(iter = _V::begin(); iter != _V::end(); iter++){
    if((*iter).name() == string(name) && (*iter).site() == site &&
       (*iter).type() == type && (*iter).internal_site() == internal_site) {
      return &(*iter);
    }
  }
#ifdef WITH_WARNINGS
  cout << "*** WARNING: Operator 2: " << name << "(" << site << "," << internal_site << ") not found\n";
#endif // WITH_WARNINGS
  return NULL; 
}

template <class T>
inline bool 
Block<T>::contains(const char *name, int site, int internal_site, size_t type) const
{
  const_iterator iter;
  for(iter = _V::begin(); iter != _V::end(); iter++){
    if((*iter).name() == string(name) && (*iter).site() == site &&
       (*iter).type() == type && (*iter).internal_site() == internal_site) {
      return true;
    }
  }
  return false; 
}

template <class T>
inline bool 
Block<T>::contains(const string &name, int site, int internal_site, size_t type) const
{
  const_iterator iter;
  for(iter = _V::begin(); iter != _V::end(); iter++){
    if((*iter).name() == name && (*iter).site() == site &&
       (*iter).type() == type && (*iter).internal_site() == internal_site) {
      return true;
    }
  }
  return false; 
}

template <class T>
inline bool
Block<T>::contains(const BasicOp<T> & op) const
{
  const_iterator iter;
  for(iter = _V::begin(); iter != _V::end(); iter++){
    if((*iter).name() == op.name() && (*iter).site() == op.site() &&
       (*iter).internal_site() == op.internal_site()) return true;
  }
  return false; 
}

template <class T>
inline bool
Block<T>::contains(const BasicOp<T> * op) const
{
  const_iterator iter;
  for(iter = _V::begin(); iter != _V::end(); iter++){
    if((*iter).name() == op->name() && (*iter).site() == op->site() &&
       (*iter).internal_site() == op->internal_site()) return true;
  }
  return false; 
}


template <class T>
inline const BasicOp<T>*
Block<T>::operator()(const BasicOp<T> & op) const
{
  const_iterator iter;
  for(iter = _V::begin(); iter != _V::end(); iter++){
    if((*iter).name() == op.name() && (*iter).site() == op.site() &&
       (*iter).internal_site() == op.internal_site()) return &(*iter);
  }
#ifdef WITH_WARNINGS
  cout << "*** WARNING: Operator 3: " << op.name() << "(" << op.site() << "," << op.internal_site() << ") not found\n";
//  for(iter = _V::begin(); iter != _V::end(); iter++)
//   cout << (*iter).name() << " " << (*iter).site() << endl;
#endif // WITH_WARNINGS
  return NULL; 
}

template <class T>
inline BasicOp<T>*
Block<T>::operator()(const BasicOp<T> & op) 
{
  iterator iter;
  for(iter = _V::begin(); iter != _V::end(); iter++){
    if((*iter).name() == op.name() && (*iter).site() == op.site() &&
       (*iter).internal_site() == op.internal_site()) return &(*iter);
  }
#ifdef WITH_WARNINGS
  cout << "*** WARNING: Operator 4: " << op.name() << "(" << op.site() << "," << op.internal_site() << ") not found\n";
#endif // WITH_WARNINGS
  return NULL; 
}

template<class T>
void
Block<T>::resize_basis(const Basis &b)
{
  // New basis

  _basis = b;
  _basis.reorder();

  // resize operators

  iterator iter;
  for(iter = _V::begin(); iter != _V::end(); iter++)
    (*iter).resize(_basis);
}

} // namespace dmtk

#endif // __DMTK_BLOCK_H__
