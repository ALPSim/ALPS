#ifndef __DMTK_HAMI_H__
#define __DMTK_HAMI_H__

#include "operators.h"
#include "block.h"

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
    typedef typename std::vector<Term<T> > _V;
    typedef typename std::vector<Term<T> >::const_iterator const_iterator;
    typedef typename std::vector<Term<T> >::iterator iterator;

    Block<T> site;
    Vector<Block<T>* > sites;

    Hami(const char *__name = NULL) : _val(0), _use_hc(false)
      { if(__name) _name = std::string(__name); else _name = std::string("H"); };
    Hami(const Lattice& lattice, const char *__name = NULL): _lattice(lattice), _val(0), _use_hc(false) 
      { 
        sites.resize(lattice.size()); 
        if(__name) _name = std::string(__name); else _name = std::string("H"); 
      };
    Hami(const Hami& h): _V(h), _lattice(h._lattice), site(h.site), sites(h.sites), _name(h._name), _val(h._val), _use_hc(h._use_hc)
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
/*
    ~Hami() 
      { 
        for(int i = 0; i < _lattice.size(); i++)
          { if(sites[i]) delete sites[i]; sites[i] = NULL; }
      }
*/
    Hami& operator+=(const Term<T>& t) 
      { push_back(t); return *this; }
    Hami& operator-=(const Term<T>& t) 
      { t.coef() = - t.coef(); push_back(t); return *this; }
    Hami& operator+=(const BasicOp<T>& t) 
      { push_back(Term<T>(t)); return *this; }
    Hami& operator-=(const BasicOp<T>& t) 
      { push_back(Term<T>(t,T(-1))); return *this; }
    Hami& operator+=(const Hami& h) 
      {
        const_iterator iter;
        for(iter = h.begin(); iter != h.end(); iter++)
          push_back(*iter);
        return *this;
      }
    Hami& operator*=(const T& v) 
      {
        iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++) (*iter) *= v;
        return *this;
      }
    Hami operator*(const Hami& h) const
      {
        Hami<T> res(*this);
        res.clear();
        const_iterator iter1;
        const_iterator iter2;
        for(iter1 = _V::begin(); iter1 != _V::end(); iter1++)
          for(iter2 = h.begin(); iter2 != h.end(); iter2++)
            res += ((*iter1)*(*iter2));
        return res;
      }
    Hami operator*(const T& v) const
      {
        Hami<T> res(*this);
        res.clear();
        const_iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++)
          res += (*iter)*v;
        return res;
      }
    Hami operator+(const Hami& h) const
      {
        Hami<T> res(*this);
        const_iterator iter;
        for(iter = h.begin(); iter != h.end(); iter++)
          res += (*iter);
        return res;
      }
 
    const Lattice &lattice() const { return _lattice; }
    Lattice &lattice() { return _lattice; }
    Block<T>& get_site(int pos) 
      { if(sites[pos]) return *sites[pos]; else { return site;} }
    const Block<T>& get_site(int pos) const
      { if(sites[pos]) return *sites[pos]; else { return site;} }

    std::string name() const { return _name; }; 
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
  typename Hami<T>::iterator iter;
  for(iter = _V::begin(); iter != _V::end(); iter++){
    Term<T> &t = *iter;
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
  const_iterator iter;
  iter = _V::begin();
  iter++;
  for(; iter != _V::end(); iter++){
    bf << " + " << (*iter).description();
  }

  return bf.str();
}


} // namespace dmtk

#endif // __DMTK_HAMI_H__
