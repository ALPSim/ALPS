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

#ifndef __DMTK_OPERATORS_H__
#define __DMTK_OPERATORS_H__

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <string>
#include <iosfwd>
#include "enums.h"
#include "conj.h"
#include "vector.h"
#include "qn.h"
#include "subspace.h"
#include "block_matrix.h"
#include "state.h"
#include "globals.h"

namespace dmtk
{

// Basic operators
template <class T> class Block;
template <class T> class Term;
template <class T> class Hami;

enum
{  
  OP_SYSTEM,
  OP_MEASURE,
  OP_ADDITIVE,
};

template <class T>
class BasicOp: public dmtk::BMatrix<T>
{
  protected:
    int _site;     // site where the operator is seated on
    int _label; // if we use more than one site per unit cell;
    int _internal_site;    // any arbitrary internal_site we want to use
                   // in general: relative position for a product operator
    bool _is_fermion;

    size_t _type;

    std::string _name;

    virtual void init() {};
    void init_blocks();
    void bmatrix_from_matrix(const Matrix<T> &);

    Term<T> _term; // if this is a composite operator originated in a product
    bool _is_term;
    Hami<T> _hami; // if this is a sum of terms
    bool _is_hami;

  public:
    typedef BMatrix<T> _BM;
    typedef std::list<SubMatrix<T> > _BMlist;

    QN dqn; // (Sz(after) - Sz(before))
            // (Nt(after) - Nt(before))

    BasicOp(): _site(0), _label(-1), _internal_site(0), _is_fermion(false), _type(OP_SYSTEM), _is_term(false), _is_hami(false) {}

    BasicOp(const char *name, int site, int internal_site = 0, int label = -1)
    : _site(site),  _label(label), _internal_site(internal_site), _is_fermion(false), _type(OP_SYSTEM), _name(name), _is_term(false), _is_hami(false) {}
    BasicOp(const char *name, const Basis& b, int site, int internal_site = 0, int label = -1): 
      BMatrix<T>(b), _site(site), _label(label), _internal_site(internal_site), _is_fermion(false), _type(OP_SYSTEM), _name(name),  _is_term(false), _is_hami(false) {}
    BasicOp(const char *name, const Matrix<T> &_m, const Basis& b, int site, int internal_site = 0, int label = -1):
       BMatrix<T>(b), _site(site),  _label(label), _internal_site(internal_site), _is_fermion(false), _type(OP_SYSTEM), _name(name), _is_term(false), _is_hami(false) {}
    BasicOp(const BasicOp<T> &o): BMatrix<T>(o), _site(o._site), _label(o._label), _internal_site(o._internal_site), _is_fermion(o._is_fermion), _type(o._type), _name(o._name), _term(o._term), _is_term(o._is_term), _hami(o._hami), _is_hami(o._is_hami), dqn(o.dqn) {}


    BasicOp(const Term<T>& t): _type(OP_SYSTEM)
      {
        _name = string(t.name());
        _site = t[0].site();
        _internal_site = t[0].internal_site();
        _term = t;
        _is_hami = false;
        _is_term = true;
        _is_fermion = false;
        typename vector<BasicOp<T> >::const_iterator iter;
        int nfermions = 0;
        for(iter = t.begin(); iter != t.end(); iter++){
          dqn += (*iter).dqn;
          if(iter->fermion()) nfermions++;
        }
        if(SGN(nfermions) == -1) _is_fermion = true;
      }

    BasicOp(const Hami<T> &h): _type(OP_ADDITIVE)
      {
        _name = string(h.name());
        _site = 0;
        _internal_site = 0;
        _hami = h;
        _is_hami = true;
        _is_fermion = false;
        _is_term = false;
        BasicOp<T> first_op = *(h.begin());
        dqn = first_op.dqn;
        _is_fermion = first_op.fermion();
      }

    BasicOp(const BasicOp<T> &op1, const BasicOp<T> &op2);
    virtual ~BasicOp() {}

//  Asignment

    BasicOp& operator=(const T& v)
      {BMatrix<T>::operator=(v); return *this;}
    BasicOp& operator=(const BMatrix<T> &m)
      {BMatrix<T>::operator=(m); return *this;}
    BasicOp& operator=(const BasicOp<T> &m)
      {
        BMatrix<T>::operator=(m); 
        dqn = m.dqn; 
        _type = m._type;
        _site = m._site;
        _internal_site = m._internal_site;
        _label = m._label;
        _name = m._name;
        _is_term = m._is_term;
        _is_fermion = m._is_fermion;
        _term = m._term;
        _is_hami = m._is_hami;
        _hami = m._hami;
        return *this;
      }
    BasicOp& operator=(const Matrix<T> &m)
      { bmatrix_from_matrix(m); return *this; }
    BasicOp& set_matrix(const Matrix<T> &m, const Basis &b, QN _dqn)
      {  dqn = _dqn; _BM::repack(b); bmatrix_from_matrix(m); return *this; }

    BasicOp& operator=(const Term<T>& t)
      {
        _type = OP_SYSTEM;
        _name = string(t.name());
        _site = t[0].site();
        _internal_site = t[0].internal_site();
        _term = t;
        _is_hami = false;
        _is_term = true;
        _is_fermion = false;
        typename vector<BasicOp<T> >::const_iterator iter;
        int nfermions = 0;
        for(iter = t.begin(); iter != t.end(); iter++){
          dqn += (*iter).dqn;
          if(iter->_is_fermion) nfermions++;
        }
        if(SGN(nfermions) == -1) _is_fermion = true;
        return *this;
      }

    Term<T> operator*(T v) const 
      { return Term<T>(*this,v); }
    Term<T> operator*(const BasicOp<T> &other) const 
      { return Term<T>(*this,other); }


    bool operator==(const BasicOp<T>& op) const
      { return (op._name == _name && op._site == _site && op._internal_site == _internal_site && op._label == _label); } 
    bool operator==(const Term<T>& t) const
      { return (name() == t.name()); }

//  Methods

    BasicOp& resize(const Basis &b) 
      { BMatrix<T>::resize(b); _BMlist::clear(); init_blocks(); return *this; }


    bool fermion() const { return (dqn.fermion_sign() == -1); }
    bool apply_sign() const { return fermion(); }
    bool is_diagonal() const { return (dqn.equal(QN(),QN::default_mask())); }

    std::string const& name() const { return _name; }
    int site() const { return _site; }
    int internal_site() const { return _internal_site; }
    int label() const { return _label; }
    size_t type() const { return _type; }
    bool is_term() const { return _is_term; }
    bool is_hami() const { return _is_hami; }

    const Hami<T> &hami() const { return _hami; }
    Hami<T> &hami() { return _hami; }
    const Term<T> &term() const { return _term; }
    Term<T> &term() { return _term; }

    BasicOp& set_type(size_t type) { _type = type; return *this; }
    BasicOp& set_internal_site(int internal_site) { _internal_site = internal_site; return *this; }
    BasicOp& set_fermion(bool fermion) { _is_fermion = fermion; return *this; }
    BasicOp& set_site(int site) { _site = site; return *this; }
    BasicOp& set_name(const char* name) { _name = std::string(name); return *this; }
    QN &get_dqn() { return dqn; }
    QN get_dqn() const { return dqn; }

    friend class Block<T>;
    friend class Hami<T>;
    friend class Term<T>;

    string description() const 
      {
        ostringstream bf;
        if(_is_term)
          return _term.description();
        else if(_is_hami)
          bf << "(" << _hami.description() << ")";
        else
          bf << name() << "(" << site() << "," << internal_site() << ")";
           
        return bf.str();
      }


    // Streams

    void write(std::ostream &s) const
    {
      s.write((const char *)&_site, sizeof(int));
      s.write((const char *)&_internal_site, sizeof(int));
      s.write((const char *)&_label, sizeof(int));
      s.write((const char *)&_type, sizeof(size_t));
      s.write((const char *)&_is_term, sizeof(bool));
      s.write((const char *)&_is_fermion, sizeof(bool));
      if(_is_term) _term.write(s);
      dqn.write(s);
      size_t n = _name.size()+1;
      char *t = new char[n];
      sprintf(t, "%s", _name.c_str());
      size_t l = strlen(t);
      s.write((const char *)&l, sizeof(size_t));
      s.write((const char *)t, l*sizeof(char));
      BMatrix<T>::write(s);
      delete(t);
    }

    void read(std::istream &s)
    {
      s.read((char *)&_site, sizeof(int));
      s.read((char *)&_internal_site, sizeof(int));
      s.read((char *)&_label, sizeof(int));
      s.read((char *)&_type, sizeof(size_t));
      s.read((char *)&_is_term, sizeof(bool));
      s.read((char *)&_is_fermion, sizeof(bool));
      if(_is_term) _term.read(s);
      dqn.read(s);
      size_t l;
      s.read((char *)&l, sizeof(size_t));
      int l1 = l+1;
      char *t = new char[l1];
      s.read((char *)t, l*sizeof(char));
      t[l] = '\0';
      _name = string(t);
      BMatrix<T>::read(s);
      delete(t);
    }

};
template <class T>
class NullOp: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init() 
      { _O::dqn = 0; }
  public:
    NullOp(): BasicOp<T>("NULL",0) 
      { init(); }
    NullOp(int site, int internal_site = 0): BasicOp<T>("NULL",site,internal_site)  
      { init(); }
    NullOp(const Basis& b, int site, int internal_site = 0): BasicOp<T>("NULL",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    NullOp(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("NULL",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
 
template <class T>
class Identity: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init() 
      { _O::dqn = 0; }
  public:
    Identity(): BasicOp<T>("I",0) 
      { init(); }
    Identity(int site, int internal_site = 0): BasicOp<T>("I",site,internal_site)  
      { init(); }
    Identity(const Basis& b, int site, int internal_site = 0): BasicOp<T>("I",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Identity(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("I",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
    
template <class T>
class H: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init() 
      { _O::_is_fermion = false; _O::dqn = 0; _O::_type = OP_ADDITIVE; }
  public:
    H(): BasicOp<T>("H",0) 
      { init(); }
    H(int site, int internal_site = 0): BasicOp<T>("H",site,internal_site)  
      { init(); }
    H(const Basis& b, int site, int internal_site = 0): BasicOp<T>("H",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    H(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("H",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};

/*    
template <class T>
class Sz: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = false; _O::dqn = 0; }
  public:
    Sz(): BasicOp<T>("Sz",0) 
      { init(); }
    Sz(int site, int internal_site = 0): BasicOp<T>("Sz",site,internal_site)  
      { init(); }
    Sz(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Sz",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Sz(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Sz",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
    
template <class T>
class N: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = false; _O::dqn = 0; }
  public:
    N(): BasicOp<T>("N",0) 
      { init(); }
    N(int site, int internal_site = 0): BasicOp<T>("N",site,internal_site)  
      { init(); }
    N(const Basis& b, int site, int internal_site = 0): BasicOp<T>("N",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    N(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("N",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
    
template <class T>
class Splus: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = false; _O::dqn.sz() = 2; }
  public:
    Splus(): BasicOp<T>("S+",0) 
      { init(); }
    Splus(int site, int internal_site = 0): BasicOp<T>("S+",site,internal_site)  
      { init(); }
    Splus(const Basis& b, int site, int internal_site = 0): BasicOp<T>("S+",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Splus(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("S+",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
    
template <class T>
class Sminus: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = false; _O::dqn.sz() = -2; }
  public:
    Sminus(): BasicOp<T>("S-",0) 
      { init(); }
    Sminus(int site, int internal_site = 0): BasicOp<T>("S-",site,internal_site)  
      { init(); }
    Sminus(const Basis& b, int site, int internal_site = 0): BasicOp<T>("S-",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Sminus(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("S-",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
    
template <class T>
class Cup: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = true; _O::dqn.sz() = -1; _O::dqn.n() = -1; }
  public:
    Cup(): BasicOp<T>("Cup",0) 
      { init(); }
    Cup(int site, int internal_site = 0): BasicOp<T>("Cup",site,internal_site)  
      { init(); }
    Cup(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cup",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Cup(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cup",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
    
template <class T>
class Cdup: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()  
      { _O::_is_fermion = true; _O::dqn.sz() = 1; _O::dqn.n() = 1; }
  public:
    Cdup(): BasicOp<T>("Cdup",0) 
      { init(); }
    Cdup(int site, int internal_site = 0): BasicOp<T>("Cdup",site,internal_site)  
      { init(); }
    Cdup(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cdup",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Cdup(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cdup",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
    
template <class T>
class Cdn: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = true; _O::dqn.sz() = 1; _O::dqn.n() = -1; }
  public:
    Cdn(): BasicOp<T>("Cdn",0) 
      { init(); }
    Cdn(int site, int internal_site = 0): BasicOp<T>("Cdn",site,internal_site)  
      { init(); }
    Cdn(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cdn",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Cdn(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cdn",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
    
template <class T>
class Cddn: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = true; _O::dqn.sz() = -1; _O::dqn.n() = 1; }
  public:
    Cddn(): BasicOp<T>("Cddn",0) 
      { init(); }
    Cddn(int site, int internal_site = 0): BasicOp<T>("Cddn",site,internal_site)  
      { init(); }
    Cddn(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cddn",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Cddn(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cddn",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
    
// Boson creation and anihilation operators    
template <class T>
class B: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = false; _O::dqn.sz() = 0; _O::dqn.n() = -1; }
  public:
    B(): BasicOp<T>("B",0) 
      { init(); }
    B(int site, int internal_site = 0): BasicOp<T>("B",site,internal_site)  
      { init(); }
    B(const Basis& b, int site, int internal_site = 0): BasicOp<T>("B",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    B(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("B",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
    
template <class T>
class Bd: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = false; _O::dqn.sz() = 0; _O::dqn.n() = 1; }
  public:
    Bd(): BasicOp<T>("Bd",0) 
      { init(); }
    Bd(int site, int internal_site = 0): BasicOp<T>("Bd",site,internal_site)  
      { init(); }
    Bd(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Bd",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Bd(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Bd",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};

*/
 
template<class T>
inline void
BasicOp<T>::init_blocks()
{
  typename PackedBasis::iterator iter;

  typedef BMatrix<T> bm;

  for(iter = bm::subspace_begin(); iter != bm::subspace_end(); iter++){
    SubSpace col_range = *iter;
    SubSpace row_range = bm::subspace(col_range.qn()+dqn);
    if(row_range.begin() == DMTK_ERROR) continue;
    SubMatrix<T> b(col_range.qn(),col_range,row_range);
    dmtk::BMatrix<T>::push_back(b);
  }
}

template<class T>
void
BasicOp<T>::bmatrix_from_matrix(const Matrix<T>& m)
{
  typename PackedBasis::iterator iter;

  typedef BMatrix<T> bm;
  _BMlist::clear();

  for(iter = bm::subspace_begin(); iter != bm::subspace_end(); iter++){
    SubSpace col_range = *iter;
    SubSpace row_range = bm::subspace(col_range.qn()+dqn);
    if(row_range.begin() == DMTK_ERROR) continue;
    SubMatrix<T> b(col_range.qn(),col_range,row_range);
    b = m(col_range,row_range);
    dmtk::BMatrix<T>::push_back(b);
  }
}

/////////////////////////////////////////////////////////////////////
// Product of two or more operators and a coefficient
//

typedef enum
{
  TERM_EMPTY,
  TERM_LOCAL,
  TERM_PRODUCT,
  TERM_EXTERN,
  TERM_MEASURE,
} TermType;

template <class T>
class Term: public std::vector<BasicOp<T> >
{
  private:
    T _coef;
    TermType _term_type;
    T _value;  // mean value, after measurement

  public:
    typedef typename std::vector<BasicOp<T> > _V;
    typedef typename _V::const_iterator const_iterator;
    typedef typename _V::iterator iterator;

    Term(): _coef(1), _term_type(TERM_EMPTY), _value(T(0)) {}

    Term(const BasicOp<T>& local_op): _coef(1), _term_type(TERM_PRODUCT), _value(T(0)) 
      { 
        if(local_op.is_term())
          operator=(local_op.term());
        else
          { _V::clear(); std::vector<BasicOp<T> >::push_back(local_op); }
      }

    Term(const BasicOp<T>& local_op, T coef): _coef(coef), _term_type(TERM_PRODUCT), _value(T(0)) 
      { _V::clear(); std::vector<BasicOp<T> >::push_back(local_op); }

    Term(const BasicOp<T>& op1, const BasicOp<T>& op2): _coef(1), _term_type(TERM_PRODUCT), _value(T(0)) 
      { _V::clear(); std::vector<BasicOp<T> >::push_back(op2); std::vector<BasicOp<T> >::push_back(op1); }

    Term(const BasicOp<T>& _op1, const BasicOp<T>& _op2, T coef): _coef(coef), _term_type(TERM_PRODUCT), _value(T(0))
          { _V::clear(); _V::push_back(_op2); _V::push_back(_op1); }

    Term(const Term& t): _V(t), _coef(t.coef()), _term_type(t.type()), _value(t.value()) {}

    Term& operator=(const Term& t) { _V::operator=(t); _coef = t.coef(); _term_type = t.type(); _value = t.value(); return *this; }

    Term operator*(const T& v) const 
      { Term aux(*this); aux._coef *= v; return aux; }

    Term operator*(const BasicOp<T>& op) const
      { 
        Term aux(op); 
        const_iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++) aux.push_back(*iter);
        aux._term_type = TERM_PRODUCT; 
        return aux; 
      }
      
    Term operator*(const Term& t) const
      {
        Term aux = t;
        const_iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++) aux.push_back(*iter);        aux._coef *= _coef;
        aux._term_type = TERM_PRODUCT;
        return aux;
      }

    Term& operator*=(const T& v) { _coef *= v; return *this; }
    Term& operator*=(const BasicOp<T>& op) { _V::push_back(op); _term_type = TERM_PRODUCT; return *this; }
    Term& operator*=(const Term& t) 
    { 
      const_iterator iter;
      for(iter = t.begin(); iter != t.end(); iter++) push_back(*iter);
      _coef *= t.coef();
      _term_type = TERM_PRODUCT; 
      return *this; 
    }

    QN dqn() const
    {
      QN _dqn;
      const_iterator iter;
      for(iter = _V::begin(); iter != _V::end(); iter++) _dqn += (*iter).dqn;
      return _dqn; 
    }

    Term reorder(bool use_sign = true) const
    {
      Term new_term;
      if(_V::size() == 1) {
        new_term = *this;
        BasicOp<T> &op = new_term[0];
        if(op.is_hami()){
          op.hami().reorder_terms();
          return new_term;
        } else if(op.is_term()){
          op.term() = op.term().reorder();
          return new_term;
        } else {
          return new_term;
        }
      }

      Term t = *this;
      Vector<int> sites(t.size());
      Vector<size_t> indx(t.size());
      Vector<char> is_fermion(t.size());
      bool fermion = false;
      for(int i = 0; i < t.size(); i++){
        const BasicOp<T> &op = t[i];
        sites(i) = op.site();
        is_fermion(i) = op.fermion();
        if(op.fermion()) fermion = true;
      }
      int nex = indexx2<int, Vector<int> >(t.size(), sites, indx, false);
      if(nex > 0){
        int fermion_sign = 1;
        for(int i = 0; i < t.size(); i++){
          if(is_fermion(indx(i))){
            for(int j = 0; j < indx(i); j++){
              if(is_fermion(indx(j))) fermion_sign = -fermion_sign;
            }
          }
          new_term *= t[indx(i)];
          is_fermion[indx(i)] = false;
        }
        new_term.coef() = fermion && use_sign ? t.coef()*T(fermion_sign) : t.coef();
      } else {
        new_term = *this;
      }
      new_term._term_type = _term_type;
      return new_term;
    }

    Term& set_type(TermType type) { _term_type = type; return *this; }

    TermType type() const { return _term_type; }
    T coef() const { return _coef; }
    T& coef() { return _coef; }

    string name(bool use_coef = false) const 
      {
        ostringstream bf;
        if(use_coef) bf << coef() << "*";
        const_iterator iter;
        iter = _V::end();
        iter--;
        bf << (*iter).name() << "(" << (*iter).site() << "," << (*iter).internal_site() << ")"; 
        if(_V::size() > 1){
          while(1){
            iter--;
            const BasicOp<T> &op = *iter;
            bf << "*" << op.name() << "(" << op.site() << "," << op.internal_site() << ")"; 
            if(iter == _V::begin()) break;
          }
        }
        return bf.str();
      }

    string description() const 
      {
        ostringstream bf;
        const_iterator iter;
if(_V::size() == 0) cout << "ERROR: size == 0\n";
        if(_V::size() == 1){
          iter = _V::begin();
          bf << coef() << "*" << (*iter).description();
          return bf.str();
        }
        iter = _V::end();
        iter--;
        bf << coef() << "*" << (*iter).description();
        if(_V::size() > 1){
          while(1){
            iter--;
//            bf << "*" << (*iter).description(); 
            bf << (*iter).description(); 
            if(iter == _V::begin()) break;
          }
        }
        return bf.str();
      }


    bool operator==(const BasicOp<T>& op)
      { return (name() == op.name()); }

    bool operator==(const Term<T>& t)
      {
        typename _V::iterator iter1;
        typename _V::const_iterator iter2;
        for(iter1 = _V::begin(), iter2 = t.begin(); iter1 != _V::end() && iter2 != t.end(); iter1++, iter2++){
          if(*iter1 != *iter2) return false;
        }

        return(_coef == t._coef && _V::size() == t.size());
      }

    T value() const { return _value; }
    Term<T>& set_value(T val) { _value = val; return *this; }

// Streams

    void write(std::ostream &s) const
    {
      size_t n = _V::size();
      s.write((const char *)&_coef, sizeof(T));
      s.write((const char *)&_term_type, sizeof(TermType));
      s.write((const char *)&_value, sizeof(T));
      s.write((const char *)&n, sizeof(size_t));
      const_iterator iter;
      for(iter = _V::begin(); iter != _V::end(); iter++) (*iter).write(s);
    }

    void read(std::istream &s)
    {
      size_t n;
      s.read((char *)&_coef, sizeof(T));
      s.read((char *)&_term_type, sizeof(TermType));
      s.read((char *)&_value, sizeof(T));
      s.read((char *)&n, sizeof(size_t));
      _V::clear();
      for(int i = 0; i < n; i++){
        BasicOp<T> op;
        op.read(s);
        std::vector<BasicOp<T> >::push_back(op);
      }
    }
};

////////////////////////////////////////////////////////////////////////////

//////////// NEW

template<class T>
void
new_operator(BasicOp<T> &new_op, const BasicOp<T>&op,
             const BMatrix<T>& rho, const Basis &new_basis, 
             int position, T coef = T(1), bool rotate = true)
{
  QN qn(-999);
  typename BMatrix<T>::iterator biter;
  const SubMatrix<T> *_sub;
  Matrix<T> m, aux;
  CTimer clock;
  clock.Start();

  for(biter = new_op.begin(); biter != new_op.end(); biter++){

    SubMatrix<T> &block(*biter);
    const SubMatrix<T> *_u1 = rho.block(block.qn());
    const SubMatrix<T> *_u2 = rho.block(block.qn()+new_op.dqn);

    if(!_u1 || !_u2) continue;

    const SubMatrix<T> &u1 = *_u1;
    const SubMatrix<T> &u2 = *_u2;

    SubSpace subspacei = new_basis(block.qn());
    SubSpace subspacej = new_basis(block.qn()+new_op.dqn);
    m.reshape(subspacei.dim(),subspacej.dim()); 
    m = T(0);
    int sign = 1;

    if(position == RIGHT){

      for(int col = subspacei.begin(); col <= subspacei.end(); col++){
        State si = new_basis(col);

        if(si.qn2() != qn || qn[0] == -999){
          qn = si.qn2();
          _sub = op.block(qn);
        }

        if(!_sub) continue;
        const SubMatrix<T> &sub = *_sub;

        if(op.fermion()) sign = si.qn1().fermion_sign();

        for(int row = subspacej.begin(); row <= subspacej.end(); row++){
          State sj = new_basis(row);

          if(sj.qn2() == qn+op.dqn && sj.i1 == si.i1){
            int i2 = si.i2 - sub.col_range().begin();
            int j2 = sj.i2 - sub.row_range().begin();

            int rcol = col - subspacei.begin();
            int rrow = row - subspacej.begin();
            m(rcol,rrow) = T(sign)*coef*sub(i2,j2);
          }

        }
      }

    } else if(position == LEFT) { 

      for(int col = subspacei.begin(); col <= subspacei.end(); col++){
        State si = new_basis(col);

        if(si.qn1() != qn || qn[0] == -999){
          qn = si.qn1();
          _sub = op.block(qn);
        }

        if(!_sub) continue;
        const SubMatrix<T> &sub = *_sub;

        for(int row = subspacej.begin(); row <= subspacej.end(); row++){
          State sj = new_basis(row);

          if(sj.qn1() == qn+op.dqn && sj.i2 == si.i2){ 
            int i1 = si.i1 - sub.col_range().begin();
            int j1 = sj.i1 - sub.row_range().begin();

            int rcol = col - subspacei.begin();
            int rrow = row - subspacej.begin();
            m(rcol,rrow) = coef*sub(i1,j1);
          }

        }
      }

    }


//    cout << "Lap: " << clock.LapTime() << endl;
//    clock.Lap();

    if(rotate){
      aux.reshape(u1.cols(),m.rows());
//      aux = product(m,u1);
//      Matrix<T> tu2(u2.rows(),u2.cols());
//      tu2 = ctranspose(u2);
//      block += product(tu2,aux);
      matrix_matrix_product('N','N',m,u1,aux);
      matrix_matrix_product('C','N',u2,aux,block,T(1),T(1));
    } else {
      block += m;
    }

/*
    for(int j = 0; j < block.col_range().dim(); j++){
      for(int i = 0; i < block.row_range().dim(); i++){
        T &val = block(j,i);
//        val = 0;
        for(int l = 0; l < u1.row_range().dim(); l++){
          for(int k = 0; k < u2.row_range().dim(); k++){
            val += conj(u2(i,k))*m(l,k)*u1(j,l);
          }
        }
      }
    }
*/

//    cout << "Lap: " << clock.LapTime() << endl;
  }
}


template<class T>
void
new_operator(BasicOp<T> &new_op, const BasicOp<T>&op1, const BasicOp<T>&op2,
             size_t m1, size_t m2,  
             const BMatrix<T>& rho, const Basis& new_basis, 
             T coef = T(1), bool use_hc = false, bool rotate = true)
{
  QN qn1(-999), qn2(-999);
  typename BMatrix<T>::iterator biter;
  const SubMatrix<T> *_sub1, *_sub2;
  Matrix<T> m, aux;
  int sign = 1;

  const BasicOp<T>* pop1 = &op1;
  const BasicOp<T>* pop2 = &op2;

  bool fermion1 = false;
  bool fermion2 = false;

  if(m2 < m1){
    pop2 = &op1;
    pop1 = &op2;
    if(op1.fermion()) fermion2 = true;        
  } else if(m2 > m1) {
    if(op2.fermion()) fermion1 = true;        
  }

//  cout << op1.name() << "(" << m1 << ")" << " " << op2.name() << "(" << m2 << ")" << " " << fermion1 << " " << fermion2 << endl;

  const BasicOp<T> &_op1 = *pop1;
  const BasicOp<T> &_op2 = *pop2;
//  cout << _op1.name() << " " <<  _op2.name() << endl;


  for(biter = new_op.begin(); biter != new_op.end(); biter++){

    SubMatrix<T> &block(*biter);
    const SubMatrix<T> *_u1 = rho.block(block.qn());
    const SubMatrix<T> *_u2 = rho.block(block.qn()+new_op.dqn);

    if(!_u1 || !_u2) continue;

    const dmtk::SubMatrix<T> &u1 = *_u1;
    const dmtk::SubMatrix<T> &u2 = *_u2;

    SubSpace subspacei = new_basis(block.qn());
    SubSpace subspacej = new_basis(block.qn()+new_op.dqn);
    m.reshape(subspacei.dim(),subspacej.dim()); 
    m = T(0);
////////////////////////////////////

    for(int col = subspacei.begin(); col <= subspacei.end(); col++){ 
      State si = new_basis(col);

      if(si.qn1() != qn1 || qn1[0] == -999){
        qn1 = si.qn1();
        _sub1 = _op1.block(qn1);
      }
      if(!_sub1) continue;

      if(si.qn2() != qn2 || qn2[0] == -999){
        qn2 = si.qn2();
        _sub2 = _op2.block(qn2);
      }
      if(!_sub2 || (m1 == m2 && qn1 != qn2 + op2.dqn)) continue;

      const dmtk::SubMatrix<T> &sub1 = *_sub1;
      const dmtk::SubMatrix<T> &sub2 = *_sub2;

      int i1 = si.i1 - sub1.col_range().begin();
      int i2 = si.i2 - sub2.col_range().begin();

      if(fermion1) sign = qn1.fermion_sign();
      else if(fermion2) sign = qn1.fermion_sign()*_op2.dqn.fermion_sign();

      for(int row = subspacej.begin(); row <= subspacej.end(); row++){ 
        State sj = new_basis(row);

        if(sj.qn2() == qn2+_op2.dqn && sj.qn1() == qn1+_op1.dqn){ 
          int j1 = sj.i1 - sub1.row_range().begin();
          int j2 = sj.i2 - sub2.row_range().begin();

          int rcol = col - subspacei.begin();
          int rrow = row - subspacej.begin();
          m(rcol,rrow) += T(sign)*coef*sub1(i1,j1)*sub2(i2,j2);
        }

      }

      if(use_hc){

        for(int row = subspacej.begin(); row <= subspacej.end(); row++){ 
          State sj = new_basis(row);

          if(sj.qn2() == qn2+_op2.dqn && sj.qn1() == qn1+_op1.dqn){ 
            int j1 = sj.i1 - sub1.row_range().begin();
            int j2 = sj.i2 - sub2.row_range().begin();
      
            int rcol = col - subspacei.begin();
            int rrow = row - subspacej.begin();
            m(rrow,rcol) += T(sign)*conj(coef)*conj(sub1(i1,j1))*conj(sub2(i2,j2));
          }

        }
      }

    }

    if(rotate){   
      aux.reshape(u1.cols(),m.rows());
/*
      aux = product(m,u1);
      Matrix<T> tu2(u2.rows(),u2.cols());
      tu2 = ctranspose(u2);
      block += product(tu2,aux);
*/
      dmtk::matrix_matrix_product('N','N',m,u1,aux);
      dmtk::matrix_matrix_product('C','N',u2,aux,block,T(1),T(1));
    } else {
      block += m;
    }


/*
    for(int j = 0; j < block.col_range().dim(); j++){
      for(int i = 0; i < block.row_range().dim(); i++){
        T &val = block(j,i);
//        val = 0;
        for(int l = 0; l < u1.row_range().dim(); l++){
          for(int k = 0; k < u2.row_range().dim(); k++){
            val += conj(u2(i,k))*m(l,k)*u1(j,l);
          }
        }
      }
    }
*/

  }


}

/////////////// OLD

template<class T>
void
new_operator2(BasicOp<T> &new_op, const BasicOp<T>&op,
             const BMatrix<T>& rho, const Basis &new_basis, 
             int position, T coef = T(1))
{
  const SubMatrix<T> *_u1;
  const SubMatrix<T> *_u2;
  QN qnref(-999);
 
  typename BMatrix<T>::const_iterator iter;
  for(iter = op.begin(); iter != op.end(); iter++){

    const SubMatrix<T> &sub = *iter;
    QN qn = sub.qn();

    typename BMatrix<T>::iterator biter = new_op.begin();
    while(biter != new_op.end()){

      SubMatrix<T> &block(*biter);

      _u1 = rho.block(block.qn());
      _u2 = rho.block(block.qn()+new_op.dqn);

      if(!_u1 || !_u2) continue;

      const SubMatrix<T>& u1 = *_u1;
      const SubMatrix<T>& u2 = *_u2;

      SubSpace row_range = block.row_range();
      SubSpace col_range = block.col_range();

      int sign = 1;

      for(int i = row_range.begin(); i <= row_range.end(); i++){
        int ri = i - row_range.begin();
        for(int j = col_range.begin(); j <= col_range.end(); j++){
          int rj = j - col_range.begin();

          T &val = block(rj,ri);
//          val = 0;

          for(int l = u1.row_range().begin(); l <= u1.row_range().end(); l++){
            int l1 = new_basis(l).i1;
            int l2 = new_basis(l).i2;
            int rl = l - u1.row_range().begin();

            if(position == LEFT && new_basis(l).qn1() != qn) continue;
            if(position == RIGHT && new_basis(l).qn2() != qn) continue;

            for(int k = u2.row_range().begin(); k <= u2.row_range().end(); k++){
              int k1 = new_basis(k).i1;
              int k2 = new_basis(k).i2;
              int rk = k - u2.row_range().begin();

              if(position == RIGHT){
                if(new_basis(k).qn2() == qn+op.dqn && k1 == l1){
                  if(op.fermion()) sign = new_basis(l).qn1().fermion_sign();
                  val += T(sign)*coef*conj(u2(ri,rk))*sub(l2-sub.col_range().begin(),k2-sub.row_range().begin())*u1(rj,rl);
                }
              } else if(position == LEFT) { 
                if(new_basis(k).qn1() == qn+op.dqn && k2 == l2) 
                  val += T(sign)*coef*conj(u2(ri,rk))*sub(l1-sub.col_range().begin(),k1-sub.row_range().begin())*u1(rj,rl);
              }
            }
          }

        }
      }

      biter++;
    }
  }
}

template<class T>
void
new_operator2(BasicOp<T> &new_op, const BasicOp<T>&op1, const BasicOp<T>&op2,
             size_t m1, size_t m2,  
             const BMatrix<T>& rho, const Basis& new_basis, T coef = T(1))
{
  const BasicOp<T>* pop1 = &op1;
  const BasicOp<T>* pop2 = &op2;

  bool fermion1 = false;
  bool fermion2 = false;

  if(mask(m1,m2) == (MASK_BLOCK1|MASK_BLOCK2)){
    if(m2 < m1){
      pop2 = &op1;
      pop1 = &op2;
      if(op2.fermion()) fermion2 = true;        
    } else {
      if(op1.fermion()) fermion1 = true;        
    }
  } else if(mask(m1,m2) == (MASK_BLOCK3|MASK_BLOCK4)) {
    if(m1 < m2){
      pop2 = &op1;
      pop1 = &op2;
      if(op1.fermion()) fermion2 = true;        
    } else {
      if(op2.fermion()) fermion1 = true;        
    } 
  }

//  cout << op1.name() << "(" << m1 << ")" << " " << op2.name() << "(" << m2 << ")" << " " << fermion1 << " " << fermion2 << endl;

  const BasicOp<T> &_op1 = *pop1;
  const BasicOp<T> &_op2 = *pop2;
//  cout << _op1.name() << " " <<  _op2.name() << endl;


  QN qnref(-999);
  const SubMatrix<T> *_u1;
  const SubMatrix<T> *_u2;
  SubMatrix<T> *_block;

  typename BMatrix<T>::const_iterator iter1, iter2;

  for(iter1 = _op1.begin(); iter1 != _op1.end(); iter1++){
    for(iter2 = _op2.begin(); iter2 != _op2.end(); iter2++){

      const SubMatrix<T> &sub1 = *iter1;
      const SubMatrix<T> &sub2 = *iter2;

      QN qn1 = sub1.qn(); 
      QN qn2 = sub2.qn(); 
      QN qn = qn1 + qn2; 

      if(qn != qnref){
        
        _block = new_op.block(qn);
        if(!_block) continue;

        _u1 = rho.block(qn);
        _u2 = rho.block(qn+new_op.dqn);

        if(!_u1 || !_u2) continue;

        qnref = qn;
      }

      SubMatrix<T> &block = *_block;
      const SubMatrix<T> &u1 = *_u1;
      const SubMatrix<T> &u2 = *_u2;

      SubSpace row_range = block.row_range();
      SubSpace col_range = block.col_range();

      int sign = 1; 

      for(int i = row_range.begin(); i <= row_range.end(); i++){
        int ri = i - row_range.begin();
        for(int j = col_range.begin(); j <= col_range.end(); j++){
          int rj = j - col_range.begin();

          T &val = block(rj,ri);
//          val = 0;

          for(int l = u1.row_range().begin(); l <= u1.row_range().end(); l++){
            int l1 = new_basis(l).i1;
            int l2 = new_basis(l).i2;
            int rl = l - u1.row_range().begin();

            if(new_basis(l).qn1() != qn1) continue;
            if(new_basis(l).qn2() != qn2) continue;

            if(fermion1) sign = new_basis(l).qn1().fermion_sign();

            for(int k = u2.row_range().begin(); k <= u2.row_range().end(); k++){
              int k1 = new_basis(k).i1;
              int k2 = new_basis(k).i2;
              int rk = k - u2.row_range().begin();

              if(new_basis(k).qn1() == qn1+_op1.dqn &&
                 new_basis(k).qn2() == qn2+_op2.dqn){ 

                if(fermion2) sign = new_basis(k).qn1().fermion_sign();

                T x1 = sub1(l1-sub1.col_range().begin(), k1-sub1.row_range().begin());
                T x2 = sub2(l2-sub2.col_range().begin(), k2-sub2.row_range().begin());
                val += conj(u2(ri,rk))*coef*sign*x1*x2*u1(rj,rl);
                
              }
            }

          }


        }
      }

    }
  }

}


////////////////////////////////////////////////////////////////////////////

template<class T>
/*static*/ const SubMatrix<T>* 
get_block(const vector<const SubMatrix<T> *> &op, const QN &qn)
{
  int origin = 0;
  int end = op.size() - 1;
  int index = -1;

  while(origin <= end){
    int index_old = index;
    index = (origin + end) / 2;

    const SubMatrix<T> *_s = op[index];

    if(_s->qn() == qn) return _s; // { cout << index << endl; return _s; }
    if(_s->qn() > qn)
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

// cout <<"NOT FOUND\n";

  return 0;
}

////////////////////////////////////////////////////////////////////
// product:
////////////////////////////////////////////////////////////////////
template<class T>
void
product(const BasicOp<T> &op,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m, T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL)
{
  Matrix<T> *_maux1, *_maux2;
  if(globals){
    _maux1 = &globals->m1;
    _maux2 = &globals->m2;
  } else {
    _maux1 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
    _maux2 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
  }
  Matrix<T> &maux1 = *_maux1;
  Matrix<T> &maux2 = *_maux2;

  Vector<QN> dqn(5);
  if(!hc) 
    dqn(m) += op.dqn; 
  else 
    dqn(m) -= op.dqn;

  char do_hc = (!hc) ? 'N' : 'C';
  T real_coef = (!hc) ? coef : conj(coef);

  vector<const SubMatrix<T> *> auxv;
  typename BasicOp<T>::const_iterator citer;
  for(citer = op.begin(); citer != op.end(); citer++)
    auxv.push_back(&(*citer));

  typename VectorState<T>::const_iterator siter;
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;

    const SubMatrix<T> *_block;

    if(!hc)
      _block = get_block(auxv, ss[m].qn());
    else
      _block = get_block(auxv, ss[m].qn()-op.dqn);

    if(!_block) continue;
    const SubMatrix<T> &block(*_block);

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    int sign = 1;
    if(op.fermion()){
      for(int ib = m-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }

    switch(m){
      case BLOCK1:
        for(int i2 = 0; i2 < ss[2].dim(); i2++){
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,slice(0,res_slice.size4(),1));
            maux2.reshape(subres.size2(),subres.size1());
            maux1 = subv;
            matrix_matrix_product(do_hc,'T',static_cast<Matrix<T> >(block),maux1,maux2,T(sign)*real_coef);
            subres.transpose() += maux2.array();
          }
        }
        break;
      case BLOCK2:
        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),i3,slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),i3,slice(0,res_slice.size4(),1));
            maux1 = subv;
            maux2.reshape(subres.size2(),subres.size1());
            matrix_matrix_product(do_hc,'T',static_cast<Matrix<T> >(block),maux1,maux2,T(sign)*real_coef);
            subres.transpose() += maux2.array();
          }
        }
        break;
      case BLOCK3:
        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),slice(0,v_slice.size3(),1),i4);
            gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),slice(0,res_slice.size3(),1),i4);
            maux1 = subv;
            maux2.reshape(subres.size1(),subres.size2());
            matrix_matrix_product(do_hc,'N',static_cast<Matrix<T> >(block),maux1,maux2,T(sign)*real_coef);
            subres += maux2.array();
          }
        }
        break;
      case BLOCK4:
        for(int i2 = 0; i2 < ss[2].dim(); i2++){
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,slice(0,res_slice.size4(),1));
            maux1 = subv;
            maux2.reshape(subres.size1(),subres.size2());
            matrix_matrix_product(do_hc,'N',static_cast<Matrix<T> >(block),maux1,maux2,T(sign)*real_coef);
            subres += maux2.array();
          }
        }
        break;
    }
  }
  if(!globals){
    delete(_maux1);
    delete(_maux2);
  }
}

template<class T>
void
product(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL)
{
  int _mask = mask(m1,m2);
  int _m1, _m2;
  Vector<QN> dqn(5);
  dqn(m2) += op2.dqn;
  dqn(m1) += op1.dqn;
  Vector<QN> dqn2(5);
  dqn2(m2) += op2.dqn;
  Vector<T> *_vaux1, *_vaux2;
  Matrix<T> *_maux1, *_maux2, *_maux3;
  if(globals){
    _vaux1 = &globals->v1;
    _vaux2 = &globals->v2;
    _maux1 = &globals->m1;
    _maux2 = &globals->m2;
    _maux3 = &globals->m3;
  } else {
    _vaux1 = new Vector<T>(MIN_VECTOR_SIZE);
    _vaux2 = new Vector<T>(MIN_VECTOR_SIZE);
    _maux1 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
    _maux2 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
    _maux3 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
  }
  Vector<T> &vaux1 = *_vaux1;
  Vector<T> &vaux2 = *_vaux2;
  Matrix<T> &maux1 = *_maux1;
  Matrix<T> &maux2 = *_maux2;
  Matrix<T> &maux3 = *_maux3;

  const BasicOp<T> *_op1;
  const BasicOp<T> *_op2;
  int sign0 = 1;
  if(m1 <= m2){
    _op1 = &op1;
    _op2 = &op2;
    _m1 = m1;
    _m2 = m2;
  } else {
    _op1 = &op2;
    _op2 = &op1;
    _m1 = m2;
    _m2 = m1;
    if(_op1->fermion() && _op2->fermion()) sign0 = -1;
  }

/* Test 
  if(m1 == m2){
    VectorState<T> aux(v.b1(),v.b2(),v.b3(),v.b4(),v.qn()+op2.dqn, v.qn_mask());
    product(op2,v,aux,m2,T(1),false);
    product(op1,aux,res,m1,coef,false);
    if(hc){
      aux.set_qn_mask(v.qn()-op1.dqn, v.qn_mask());
      aux = T(0);
      product(op1,v,aux,m1,T(1),true);
      product(op2,aux,res,m2,coef,true);
    }
    return;
  }
*/
        
  vector<const SubMatrix<T> *> auxv1, auxv2;
  typename BasicOp<T>::const_iterator citer;
  for(citer = _op1->begin(); citer != _op1->end(); citer++) 
    auxv1.push_back(&(*citer));
  for(citer = _op2->begin(); citer != _op2->end(); citer++) 
    auxv2.push_back(&(*citer));

  typename VectorState<T>::const_iterator siter;
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;

    const SubMatrix<T> *_block1;
    const SubMatrix<T> *_block2;

    _block2 = get_block(auxv2, ss[_m2].qn());
    if(m1 != m2)
      _block1 = get_block(auxv1, ss[_m1].qn());
    else
      _block1 = get_block(auxv1, ss[_m1].qn() + op2.dqn);
    

    if(!_block1 || !_block2) continue;

    const SubMatrix<T> &block1(*_block1);
    const SubMatrix<T> &block2(*_block2);

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    cstate_slice<T> v_slice_hc = v(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    state_slice<T> res_slice_hc = res(ss[1].qn(),ss[2].qn(),ss[3].qn(),ss[4].qn());
    if(v_slice_hc.size() == 0 || res_slice_hc.size() == 0) continue; 

    int sign = sign0;
    if(op2.fermion()){
      for(int ib = m2-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }
    if(op1.fermion()){
      for(int ib = m1-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn2[ib].fermion_sign();
    }

    switch(_mask){
      case (MASK_BLOCK1|MASK_BLOCK2):
        for(int i3 = 0; i3 < ss[3].dim(); i3++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),
                                         slice(0,v_slice.size2(),1),i3,i4);
          gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),
                                            slice(0,res_slice.size2(),1),i3,i4);
          maux1 = subv;
          maux2 = subres;
          product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
          subres = maux2.array();
        }
        if(hc){
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),
                                              slice(0,v_slice_hc.size2(),1),i3,i4);
            gslice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),
                                                 slice(0,res_slice_hc.size2(),1),i3,i4);
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
            subres = maux2.array();
          }
        }
        break;
      case (MASK_BLOCK1|MASK_BLOCK3):
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,
                                         slice(0,v_slice.size3(),1),i4);
          gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,
                                            slice(0,res_slice.size3(),1),i4);
          maux1 = subv;
          maux2 = subres;
          product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
          subres = maux2.array();
        }
        if(hc){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),i2,
                                              slice(0,v_slice_hc.size3(),1),i4);
            gslice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),i2,
                                                 slice(0,res_slice_hc.size3(),1),i4);
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
            subres = maux2.array();
          }
        }
        break;
      case (MASK_BLOCK1|MASK_BLOCK4):
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++){
          cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,
                                         slice(0,v_slice.size4(),1));
          gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,
                                            slice(0,res_slice.size4(),1));
          maux1 = subv;
          maux2 = subres;
          product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
          subres = maux2.array();
        }
        if(hc){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),i2,i3,
                                              slice(0,v_slice_hc.size4(),1));
            gslice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),i2,i3,
                                                 slice(0,res_slice_hc.size4(),1));
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
            subres = maux2.array();
          }
        }
        break;
      case (MASK_BLOCK2|MASK_BLOCK3):
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),
                                         slice(0,v_slice.size3(),1),i4);
          gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),
                                            slice(0,res_slice.size3(),1),i4);
          maux1 = subv;
          maux2 = subres;
          product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
          subres = maux2.array();
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice_hc(i1,slice(0,v_slice_hc.size2(),1),
                                           slice(0,v_slice_hc.size3(),1),i4);
            gslice_iter<T> subres = res_slice_hc(i1,slice(0,res_slice_hc.size2(),1),
                                              slice(0,res_slice_hc.size3(),1),i4);
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
            subres = maux2.array();
          }
        }
        break;
      case (MASK_BLOCK2|MASK_BLOCK4):
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++){
          cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),i3,
                                         slice(0,v_slice.size4(),1));
          gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),i3,
                                            slice(0,res_slice.size4(),1));
          maux1 = subv;
          maux2 = subres;
          product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
          subres = maux2.array();
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice_hc(i1,slice(0,v_slice_hc.size2(),1),i3,
                                              slice(0,v_slice_hc.size4(),1));
            gslice_iter<T> subres = res_slice_hc(i1,slice(0,res_slice_hc.size2(),1),i3,
                                                 slice(0,res_slice_hc.size4(),1));
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
            subres = maux2.array();
          }
        }
        break;
      case (MASK_BLOCK3|MASK_BLOCK4):
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++){
          cgslice_iter<T> subv = v_slice(i1,i2,slice(0,v_slice.size3(),1),
                                         slice(0,v_slice.size4(),1));
          gslice_iter<T> subres = res_slice(i1,i2,slice(0,res_slice.size3(),1),
                                            slice(0,res_slice.size4(),1));
          maux1 = subv;
          maux2 = subres;
          product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
          subres = maux2.array();
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++){
            cgslice_iter<T> subv = v_slice_hc(i1,i2,slice(0,v_slice_hc.size3(),1),
                                              slice(0,v_slice_hc.size4(),1));
            gslice_iter<T> subres = res_slice_hc(i1,i2,slice(0,res_slice_hc.size3(),1),
                                                 slice(0,res_slice_hc.size4(),1));
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
            subres = maux2.array();
          }
        }
        break;
      case MASK_BLOCK1:
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,i4);
          slice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,i4);
          vaux1 = subv;
          vaux2 = subres;
          product(block1,block2,vaux1,vaux2,coef,T(1));
          subres = vaux2.array();
        }
        if(hc){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),i2,i3,i4);
            slice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),i2,i3,i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
      case MASK_BLOCK2:
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),i3,i4);
          slice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),i3,i4);
          vaux1 = subv;
          vaux2 = subres;
          product(block1,block2,vaux1,vaux2,coef,T(1));
          subres = vaux2.array();
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice_hc(i1,slice(0,v_slice_hc.size2(),1),i3,i4);
            slice_iter<T> subres = res_slice_hc(i1,slice(0,res_slice_hc.size2(),1),i3,i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
      case MASK_BLOCK3:
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cslice_iter<T> subv = v_slice(i1,i2,slice(0,v_slice.size3(),1),i4);
          slice_iter<T> subres = res_slice(i1,i2,slice(0,res_slice.size3(),1),i4);
          vaux1 = subv;
          vaux2 = subres;
          product(block1,block2,vaux1,vaux2,coef,T(1));
          subres = vaux2.array();
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice_hc(i1,i2,slice(0,v_slice_hc.size3(),1),i4);
            slice_iter<T> subres = res_slice_hc(i1,i2,slice(0,res_slice_hc.size3(),1),i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
      case MASK_BLOCK4:
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++){
          cslice_iter<T> subv = v_slice(i1,i2,i3,slice(0,v_slice.size4(),1));
          slice_iter<T> subres = res_slice(i1,i2,i3,slice(0,res_slice.size4(),1));
          vaux1 = subv;
          vaux2 = subres;
          product(block1,block2,vaux1,vaux2,coef,T(1));
          subres = vaux2.array();
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cslice_iter<T> subv = v_slice_hc(i1,i2,i3,slice(0,v_slice_hc.size4(),1));
            slice_iter<T> subres = res_slice_hc(i1,i2,i3,slice(0,res_slice_hc.size4(),1));
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
    }
  }
  if(!globals){
    delete(_vaux1);
    delete(_vaux2);
    delete(_maux1);
    delete(_maux2);
    delete(_maux3);
  }
}

template<class T>
void
product(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const BasicOp<T> &op3, const BasicOp<T>& op4,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, size_t m3, size_t m4,
        T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL)
{
  Vector<QN> dqn(5);
  size_t _m1, _m2, _m3, _m4;
  dqn(m4) += op4.dqn;
  dqn(m3) += op3.dqn;
  dqn(m2) += op2.dqn;
  dqn(m1) += op1.dqn;
  Vector<QN> dqn2(5);
  dqn2(m2) += op2.dqn;
  Vector<QN> dqn3(5);
  dqn3(m3) += op3.dqn;
  Vector<QN> dqn4(5);
  dqn4(m4) += op4.dqn;

  const BasicOp<T> *_op1;
  const BasicOp<T> *_op2;
  const BasicOp<T> *_op3;
  const BasicOp<T> *_op4;

  Vector<size_t> m(4), indx(4);
  Vector<const BasicOp<T>* > ops(4);
  ops(0) = &op1;
  ops(1) = &op2;
  ops(2) = &op3;
  ops(3) = &op4;
  m(0) = int(m1);
  m(1) = int(m2);
  m(2) = int(m3);
  m(3) = int(m4);
  indexx<size_t,Vector<size_t> >(4, m, indx);
  _m1 = (size_t)m(indx(0));
  _m2 = (size_t)m(indx(1));
  _m3 = (size_t)m(indx(2));
  _m4 = (size_t)m(indx(3));
  _op1 = ops(indx(0));
  _op2 = ops(indx(1));
  _op3 = ops(indx(2));
  _op4 = ops(indx(3));

  vector<const SubMatrix<T> *> auxv1, auxv2, auxv3, auxv4;
  typename BasicOp<T>::const_iterator citer;
  for(citer = _op1->begin(); citer != _op1->end(); citer++) 
    auxv1.push_back(&(*citer));
  for(citer = _op2->begin(); citer != _op2->end(); citer++) 
    auxv2.push_back(&(*citer));
  for(citer = _op3->begin(); citer != _op3->end(); citer++) 
    auxv3.push_back(&(*citer));
  for(citer = _op4->begin(); citer != _op4->end(); citer++) 
    auxv4.push_back(&(*citer));

  typename VectorState<T>::const_iterator siter;
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;

    const SubMatrix<T> *_block1;
    const SubMatrix<T> *_block2;
    const SubMatrix<T> *_block3;
    const SubMatrix<T> *_block4;

    _block4 = get_block(auxv4, ss[_m4].qn());
    _block3 = get_block(auxv3, ss[_m3].qn());
    _block2 = get_block(auxv2, ss[_m2].qn());
    _block1 = get_block(auxv1, ss[_m1].qn());

    if(!_block1 || !_block2 || !_block3 || !_block4) continue;

    const SubMatrix<T> &block1(*_block1);
    const SubMatrix<T> &block2(*_block2);
    const SubMatrix<T> &block3(*_block3);
    const SubMatrix<T> &block4(*_block4);

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    int sign = 1;
    if(op4.fermion()){
      for(int ib = m4-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }
    if(op3.fermion()){
      for(int ib = m3-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn4[ib].fermion_sign();
    }
    if(op2.fermion()){
      for(int ib = m2-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn3[ib].fermion_sign()*dqn4[ib].fermion_sign();
    }
    if(op1.fermion()){
      for(int ib = m1-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn2[ib].fermion_sign()*dqn3[ib].fermion_sign()*dqn4[ib].fermion_sign();
    }

    for(int i1 = 0; i1 < ss[1].dim(); i1++)
      for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++)
            for(int j1 = 0; j1 < res_slice.size1(); j1++)
              for(int j2 = 0; j2 < res_slice.size2(); j2++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(j1,j2,j3,j4) += coef*T(sign)*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i2,j2)*block3(i3,j3)*block4(i4,j4);

  }

}

template<class T>
void
product(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const BasicOp<T> &op3, 
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, size_t m3,
        T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL)
{
  Vector<QN> dqn(5);
  size_t _m1, _m2, _m3;
  dqn(m3) += op3.dqn;
  dqn(m2) += op2.dqn;
  dqn(m1) += op1.dqn;
  Vector<QN> dqn2(5);
  dqn2(m2) += op2.dqn;
  Vector<QN> dqn3(5);
  dqn3(m3) += op3.dqn;

  const BasicOp<T> *_op1;
  const BasicOp<T> *_op2;
  const BasicOp<T> *_op3;

  Vector<size_t> m(3), indx(3);
  Vector<const BasicOp<T>* > ops(3);
  ops(0) = &op1;
  ops(1) = &op2;
  ops(2) = &op3;
  m(0) = int(m1);
  m(1) = int(m2);
  m(2) = int(m3);
  indexx<size_t,Vector<size_t> >(3, m, indx);
  _m1 = (size_t)m(indx(0));
  _m2 = (size_t)m(indx(1));
  _m3 = (size_t)m(indx(2));
  _op1 = ops(indx(0));
  _op2 = ops(indx(1));
  _op3 = ops(indx(2));

  vector<const SubMatrix<T> *> auxv1, auxv2, auxv3;
  typename BasicOp<T>::const_iterator citer;
  for(citer = _op1->begin(); citer != _op1->end(); citer++) 
    auxv1.push_back(&(*citer));
  for(citer = _op2->begin(); citer != _op2->end(); citer++) 
    auxv2.push_back(&(*citer));
  for(citer = _op3->begin(); citer != _op3->end(); citer++) 
    auxv3.push_back(&(*citer));

  typename VectorState<T>::const_iterator siter;
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;

    const SubMatrix<T> *_block1;
    const SubMatrix<T> *_block2;
    const SubMatrix<T> *_block3;

    _block3 = get_block(auxv3, ss[_m3].qn());
    _block2 = get_block(auxv2, ss[_m2].qn());
    _block1 = get_block(auxv1, ss[_m1].qn());

    if(!_block1 || !_block2 || !_block3) continue;

    const SubMatrix<T> &block1(*_block1);
    const SubMatrix<T> &block2(*_block2);
    const SubMatrix<T> &block3(*_block3);

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    int sign = 1;
    if(op3.fermion()){
      for(int ib = m3-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }
    if(op2.fermion()){
      for(int ib = m2-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn3[ib].fermion_sign();
    }
    if(op1.fermion()){
      for(int ib = m1-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn2[ib].fermion_sign()*dqn3[ib].fermion_sign();
    }

    if(_m1 == BLOCK1 && _m2 == BLOCK2 && _m3 == BLOCK3){ 
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j2 = 0; j2 < res_slice.size2(); j2++)
                  for(int j3 = 0; j3 < res_slice.size3(); j3++)
                      res_slice(j1,j2,j3,i4) += coef*T(sign)*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i2,j2)*block3(i3,j3);
    } else if(_m1 == BLOCK1 && _m2 == BLOCK2 && _m3 == BLOCK4){ 
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j2 = 0; j2 < res_slice.size2(); j2++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(j1,j2,i3,j4) += coef*T(sign)*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i2,j2)*block3(i4,j4);
    } else if(_m1 == BLOCK1 && _m2 == BLOCK3 && _m3 == BLOCK4){ 
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(j1,i2,j3,j4) += coef*T(sign)*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i3,j3)*block3(i4,j4);
    } else if(_m1 == BLOCK2 && _m2 == BLOCK3 && _m3 == BLOCK4){ 
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j2 = 0; j2 < res_slice.size2(); j2++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(i1,j2,j3,j4) += coef*T(sign)*v_slice(i1,i2,i3,i4)*block1(i2,j2)*block2(i3,j3)*block3(i4,j4);
    }
  }

}

//////////////////////////////////////////////////////////////////////////

} // dmtk

#endif // __DMTK_OPERATORS_H__
