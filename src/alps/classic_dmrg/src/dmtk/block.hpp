// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/block.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef DMTK_BLOCK_H_
#define DMTK_BLOCK_H_

#include "dmtk/output.hpp"
#include <iosfwd>
#include <list>
#include "dmtk/enums.hpp"
#include "dmtk/vector.hpp"
#include "dmtk/matrix.hpp"
#include "dmtk/qn.hpp"
#include "dmtk/state.hpp"
#include "dmtk/subspace.hpp"
#include "dmtk/basis.hpp"
#include "dmtk/block_matrix.hpp"
#include "dmtk/operators.hpp"
#include "dmtk/lattice.hpp"

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

    // Shared linear scan backing the operator()/contains lookup overloads.
    // Each overload supplies its own predicate, preserving the deliberately
    // non-uniform match criteria (the by-name overloads test type(), the
    // by-op overloads do not).
    template <class Pred>
    BasicOp<T>* find_op(Pred pred)
    {
      for (BasicOp<T> & op : *this) if (pred(op)) return &op;
      return nullptr;
    }
    template <class Pred>
    const BasicOp<T>* find_op(Pred pred) const
    {
      for (BasicOp<T> const & op : *this) if (pred(op)) return &op;
      return nullptr;
    }
  public:
    using _V = typename std::list<BasicOp<T> >;
    using iterator = typename std::list<BasicOp<T> >::iterator;
    using const_iterator = typename std::list<BasicOp<T> >::const_iterator;

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
        for (BMatrix<T> & op : *this) op = v;
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
      s.write(reinterpret_cast<const char*>(&l), sizeof(size_t));
      s.write(reinterpret_cast<const char*>(&_n_orbitals), sizeof(size_t));
      s.write(reinterpret_cast<const char*>(&_single_site), sizeof(bool));

      _basis.write(s);
      _lattice.write(s);

      for (BasicOp<T> const & op : *this) op.write(s);
    }

    void read(std::istream &s)
    {
      size_t l;
      s.read(reinterpret_cast<char*>(&l), sizeof(size_t));
      s.read(reinterpret_cast<char*>(&_n_orbitals), sizeof(size_t));
      s.read(reinterpret_cast<char*>(&_single_site), sizeof(bool));
      list<BasicOp<T> >::resize(l);

      _basis.read(s);
      _lattice.read(s);

      for (BasicOp<T> & op : *this) op.read(s);
    }

};

// Predicate for the by-name lookups: matches name, site, type AND
// internal_site (the by-name overloads deliberately test type()).
template <class T>
inline BasicOp<T>*
Block<T>::operator()(const char *name, int site, int internal_site, size_t type)
{
  BasicOp<T>* r = find_op([&](BasicOp<T> const & op) {
    return op.name() == string(name) && op.site() == site &&
           op.type() == type && op.internal_site() == internal_site; });
#ifdef WITH_WARNINGS
  if (!r)
    dmtk::out() << "*** WARNING: Operator 1: " << name << "(" << site << "," << internal_site << ") not found\n";
#endif // WITH_WARNINGS
  return r;
}

template <class T>
inline const BasicOp<T>*
Block<T>::operator()(const char *name, int site, int internal_site, size_t type) const
{
  const BasicOp<T>* r = find_op([&](BasicOp<T> const & op) {
    return op.name() == string(name) && op.site() == site &&
           op.type() == type && op.internal_site() == internal_site; });
#ifdef WITH_WARNINGS
  if (!r)
    dmtk::out() << "*** WARNING: Operator 2: " << name << "(" << site << "," << internal_site << ") not found\n";
#endif // WITH_WARNINGS
  return r;
}

template <class T>
inline bool
Block<T>::contains(const char *name, int site, int internal_site, size_t type) const
{
  return find_op([&](BasicOp<T> const & op) {
    return op.name() == string(name) && op.site() == site &&
           op.type() == type && op.internal_site() == internal_site; }) != nullptr;
}

template <class T>
inline bool
Block<T>::contains(const string &name, int site, int internal_site, size_t type) const
{
  return find_op([&](BasicOp<T> const & op) {
    return op.name() == name && op.site() == site &&
           op.type() == type && op.internal_site() == internal_site; }) != nullptr;
}

// Predicate for the by-op lookups: matches name, site AND internal_site
// only (the by-op overloads deliberately do NOT test type()).
template <class T>
inline bool
Block<T>::contains(const BasicOp<T> & op) const
{
  return find_op([&](BasicOp<T> const & o) {
    return o.name() == op.name() && o.site() == op.site() &&
           o.internal_site() == op.internal_site(); }) != nullptr;
}

template <class T>
inline bool
Block<T>::contains(const BasicOp<T> * op) const
{
  return find_op([&](BasicOp<T> const & o) {
    return o.name() == op->name() && o.site() == op->site() &&
           o.internal_site() == op->internal_site(); }) != nullptr;
}


template <class T>
inline const BasicOp<T>*
Block<T>::operator()(const BasicOp<T> & op) const
{
  const BasicOp<T>* r = find_op([&](BasicOp<T> const & o) {
    return o.name() == op.name() && o.site() == op.site() &&
           o.internal_site() == op.internal_site(); });
#ifdef WITH_WARNINGS
  if (!r)
    dmtk::out() << "*** WARNING: Operator 3: " << op.name() << "(" << op.site() << "," << op.internal_site() << ") not found\n";
#endif // WITH_WARNINGS
  return r;
}

template <class T>
inline BasicOp<T>*
Block<T>::operator()(const BasicOp<T> & op)
{
  BasicOp<T>* r = find_op([&](BasicOp<T> const & o) {
    return o.name() == op.name() && o.site() == op.site() &&
           o.internal_site() == op.internal_site(); });
#ifdef WITH_WARNINGS
  if (!r)
    dmtk::out() << "*** WARNING: Operator 4: " << op.name() << "(" << op.site() << "," << op.internal_site() << ") not found\n";
#endif // WITH_WARNINGS
  return r;
}

template<class T>
void
Block<T>::resize_basis(const Basis &b)
{
  // New basis

  _basis = b;
  _basis.reorder();

  // resize operators
  for (BasicOp<T> & op : *this) op.resize(_basis);
}

} // namespace dmtk

#endif // DMTK_BLOCK_H_
