// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/hami.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef DMTK_HAMI_H_
#define DMTK_HAMI_H_

#include "dmtk/operators.hpp"
#include "dmtk/block.hpp"

namespace dmtk{

template <class T>
class Hami : public std::vector<Term<T> > 
{
  private:
    Lattice _lattice;
    std::string _name;
    T _val;
    bool _use_hc;
  public:
    using _V = typename std::vector<Term<T> >;
    using const_iterator = typename std::vector<Term<T> >::const_iterator;
    using iterator = typename std::vector<Term<T> >::iterator;

    Block<T> site;
    Vector<Block<T>* > sites;

    Hami(const char *__name = nullptr) : _val(0), _use_hc(false)
      { if(__name) _name = std::string(__name); else _name = std::string("H"); };
    Hami(const Lattice& lattice, const char *__name = nullptr): _lattice(lattice), _val(0), _use_hc(false) 
      { 
        sites.resize(lattice.size()); 
        if(__name) _name = std::string(__name); else _name = std::string("H"); 
      };
    Hami(const Hami& h): _V(h), _lattice(h._lattice), _name(h._name), _val(h._val), _use_hc(h._use_hc), site(h.site), sites(h.sites)
      { sites.resize(_lattice.size()); } 

    Hami& operator=(const Hami<T> &h)
      {
        _V::operator=(h);
        _lattice = h._lattice;
        site = h.site;
        sites = h.sites;
        _name = h._name;
        _val = h._val;
        _use_hc = h._use_hc;
        return *this;
      }
    Hami& operator+=(const Term<T>& t)
      { std::vector<Term<T> >::push_back(t); return *this; }
    Hami& operator-=(const Term<T>& t) 
      { t.coef() = - t.coef(); std::vector<Term<T> >::push_back(t); return *this; }
    Hami& operator+=(const BasicOp<T>& t) 
      { std::vector<Term<T> >::push_back(Term<T>(t)); return *this; }
    Hami& operator-=(const BasicOp<T>& t) 
      { std::vector<Term<T> >::push_back(Term<T>(t,T(-1))); return *this; }
    Hami& operator+=(const Hami& h)
      {
        for (Term<T> const & term : h)
          std::vector<Term<T> >::push_back(term);
        return *this;
      }
    Hami& operator*=(const T& v)
      {
        for (Term<T> & term : *this) term *= v;
        return *this;
      }
    Hami operator*(const Hami& h) const
      {
        Hami<T> res(*this);
        res.clear();
        for (Term<T> const & a : *this)
          for (Term<T> const & b : h)
            res += a * b;
        return res;
      }
    Hami operator*(const T& v) const
      {
        Hami<T> res(*this);
        res.clear();
        for (Term<T> const & term : *this)
          res += term * v;
        return res;
      }
    Hami operator+(const Hami& h) const
      {
        Hami<T> res(*this);
        for (Term<T> const & term : h)
          res += term;
        return res;
      }
 
    const Lattice &lattice() const { return _lattice; }
    Lattice &lattice() { return _lattice; }
    Block<T>& get_site(int pos) 
      { if(sites[pos]) return *sites[pos]; else { return site;} }
    const Block<T>& get_site(int pos) const
      { if(sites[pos]) return *sites[pos]; else { return site;} }

    std::string const& name() const { return _name; }; 
    Hami& set_name(const char *name) 
      { _name = std::string(name); return *this; }; 
    T value() const { return _val; }
    Hami& set_value(T val) { _val = val; return *this; }

    Hami& set_use_hc(bool use) { _use_hc = use; return *this; }
    bool use_hc() const { return _use_hc; }

    Hami& reorder_terms(bool use_sign = true); 

    Hami find_common_factors() const;
    std::string description() const;
};

template<class T>
Hami<T> &
Hami<T>::reorder_terms(bool use_sign)
{
  for (Term<T> & t : *this) {
    t = t.reorder(use_sign);
/*
    Vector<int> sites(t.size());
    Vector<size_t> indx(t.size());
    bool fermion = false;
    for(int i = 0; i < t.size(); i++){
      BasicOp<T> &op = t[i];
      sites(i) = op.site();
      if(op.fermion()) fermion = true;
    }
    int nex = indexx2<int, Vector<int> >(t.size(), sites, indx, false);
    if(nex > 0){
      Term<T> new_term;
      for(int i = 0; i < t.size(); i++){
        new_term *= t[indx(i)];
      }
      new_term.coef() = fermion && use_sign ? t.coef()*T(SGN(nex)) : t.coef();
      t.clear();
      t = new_term;
    }
*/
  }

  return *this;
}

template<class T>
std::string
Hami<T>::description() const
{
  ostringstream bf;

  if(_V::size() == 0) return std::string("");

  bf << _V::operator[](0).description();
  for (auto iter = std::next(_V::begin()); iter != _V::end(); ++iter) {
    bf << " + " << iter->description();
  }

  return bf.str();
}


} // namespace dmtk

#endif // DMTK_HAMI_H_
