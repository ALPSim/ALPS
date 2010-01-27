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

#ifndef __DMTK_LATTICE_H__
#define __DMTK_LATTICE_H__

#include <iosfwd>
#include <vector>
#include <valarray>
#include <stdio.h>
#include "bits.h" 

namespace dmtk{

enum
{
  LATTICE_NULL,
  LATTICE_1D,
  LATTICE_2D,
};

enum
{
  NN0,
  NN1_RIGHT,
  NN1_TOP,
  NN1_LEFT,
  NN1_BOTTOM,
  NN2_TOP_RIGHT,
  NN2_TOP_LEFT,
  NN2_BOTTOM_LEFT,
  NN2_BOTTOM_RIGHT,
  NN3_RIGHT,
  NN3_TOP,
  NN3_LEFT,
  NN3_BOTTOM,
  NN_LAST,
};

enum 
{
  OBC, // Open Boundary Conditions
  PBC, // Periodic Boundary Conditions
  CBC, // Cylindric Boundary Conditions
};


class Lattice;

class Site
{
  private:
    std::valarray<int> _nn;
    size_t _index;
    int ix, iy;
  public:
    Site(): _nn(std::valarray<int>(NN_LAST)) {} 
    Site(int x, int y = 0): _nn(std::valarray<int>(NN_LAST)),ix(x),iy(y) {} 
    Site(int x, const std::valarray<int> &nn): _nn(nn),ix(x), iy(0)  {} 
    Site(int x, int y, const std::valarray<int> &nn): _nn(nn),ix(x), iy(y) {} 

    int & operator()(size_t n) { return _nn[n]; }
    int operator()(size_t n) const { return _nn[n]; }
    
    int & operator[](size_t n) { return _nn[n]; }
    int operator[](size_t n) const { return _nn[n]; }

    int & nn(size_t n) { return _nn[n]; }
    int nn(size_t n) const { return _nn[n]; }

    int x() const { return ix; } 
    int y() const { return iy; } 
    size_t index() const { return _index; } 

    friend class Lattice;

};

class Lattice : public std::vector<Site> 
{
  private:
    size_t _type; // 1d / 2d
    int _nx, _ny;
    size_t _bc;  // boundary conditions
    bool _snake; // topolgy of the lattice for 2d

    int curr;    // current position of the pointer

    void init();
    int map2d(int x, int y) const 
      { 
        if(_snake){
          int sign = SGN(x);
          int index = x*_ny; 
          index += (1-sign)/2*(_ny-1) + sign*y; 
          return index;
        } else {
          return x*_ny+y;
        }
      }
    int iy(int idx) const 
      { 
        if(_snake){
          int x = idx/_ny;
          int sign = SGN(x);
          int i = idx - x*_ny;
          int y = (1-sign)/2*(_ny-1) + sign*i; 
          return y; 
        } else {
          return idx-(idx/_ny)*_ny;
        }
      }
    int ix(int idx) const { return idx/_ny; }

  public:
    typedef std::vector<Site>::const_iterator const_iterator;
    typedef std::vector<Site>::iterator iterator;

    Lattice(): _type(LATTICE_NULL), _bc(OBC), _snake(false), curr(-1) {};
    Lattice(int nx, size_t bc): 
      _type(LATTICE_1D), _nx(nx), _ny(1), _bc(bc), _snake(false), curr(0) { init(); }
    Lattice(int nx, int ny, size_t bc, bool snake = false): 
      _type(LATTICE_2D), _nx(nx), _ny(ny), _bc(bc), _snake(snake), curr(0){ init(); }
    Lattice(const Lattice& l): std::vector<Site>(l),
      _type(l._type), _nx(l._nx), _ny(l._ny), _bc(l._bc), _snake(l._snake), curr(l.curr)  {}

    Lattice& operator=(const Lattice& l)
      {
        std::vector<Site>::operator=(l);
        _type = l._type;
        _nx = l._nx;
        _ny = l._ny;
        _bc = l._bc;
        curr = l.curr;
        _snake = l._snake;
        return *this;
      }

    Lattice& operator+=(const Site& s) { push_back(s); return *this; }

    Lattice& push_back(const Site &s)
      {
        Site aux(s);
        switch(_type){
          case LATTICE_1D:
            aux._index = size();
            break;
          case LATTICE_2D:
            aux._index = map2d(s.x(),s.y());
            break;
        }
        std::vector<Site>::push_back(aux);
        return *this;
      }

    Lattice& operator()(size_t x) 
      { 
        std::vector<Site>::iterator iter;
        int n = 0;
        for(iter = begin(); iter != end(); iter++, n++){
          if((*iter).x() == x) { curr = n; return *this; }
        }
        std::cout << "*** WARNING: Lattice 1: Site " << x << " not found\n";
        return(*this);
      }

    Lattice& operator()(size_t x, size_t y) 
      { 
        std::vector<Site>::iterator iter;
        int n = 0;
        for(iter = begin(); iter != end(); iter++, n++){
          if((*iter).x() == x && (*iter).y() == y) 
            { curr = n; return (*this); }
        }
        std::cout << "*** WARNING: Lattice 2: Site " << x << " " << y << " not found\n";
        return(*this);
      }

    Lattice& operator[](size_t nn)
      { 
        Site s = std::vector<Site>::operator[](curr);
        size_t index = s[nn];
        int n = 0;
        std::vector<Site>::iterator iter;
        for(iter = begin(); iter != end(); iter++, n++){
          if((*iter)._index == index) { curr = n; return *this; }
        }
        return(*this);
      }

    Lattice operator()(size_t x) const 
      { 
        Lattice aux(*this);
        std::vector<Site>::iterator iter;
        int n = 0;
        for(iter = aux.begin(); iter != aux.end(); iter++, n++){
          if((*iter).x() == x) { aux.curr = n; return aux; }
        }
        std::cout << "*** WARNING: Lattice 3: Site " << x << " not found\n";
        return(aux);
      }

    Lattice operator()(size_t x, size_t y) const
      { 
        Lattice aux(*this);
        std::vector<Site>::iterator iter;
        int n = 0;
        for(iter = aux.begin(); iter != aux.end(); iter++, n++){
          if((*iter).x() == x && (*iter).y() == y) 
            { aux.curr = n; return (aux); }
        }
        std::cout << "*** WARNING: Lattice 4: Site " << x << " " << y << " not found\n";
        return(aux);
      }

    Lattice operator[](size_t nn) const
      { 
        Lattice aux(*this);
        Site s = site();
        size_t index = s[nn];
        std::vector<Site>::iterator iter;
        int n = 0;
        for(iter = aux.begin(); iter != aux.end(); iter++, n++){
          if((*iter)._index == index) { aux.curr = n; return aux; }
        }
        return(aux);
      }

    Site site() const 
      { 
        Site s;
        int i;
        for(i = 0; i < size(); i++) {
          s = std::vector<Site>::operator[](i);
          if(s.index() == curr) break;
        }
        return std::vector<Site>::operator[](i);
      }

    Site& site() 
      { 
        Site s;
        int i;
        for(i = 0; i < size(); i++) {
          s = std::vector<Site>::operator[](i);
          if(s.index() == curr) break;;
        }
        return std::vector<Site>::operator[](i);
      }

    Site site(int x) const { return operator()(x).site(); }
    Site& site(int x) { return operator()(x).site(); }

    Site site(int x, int y) const { return operator()(x,y).site(); }
    Site& site(int x, int y) { return operator()(x,y).site(); }


    int pos() const { return curr; }

    size_t type() const { return _type; }
    int nx() const { return _nx; }
    int ny() const { return _ny; }
    int& nx() { return _nx; }
    int& ny() { return _ny; }
    int x(int pos) const { return ix(pos); }
    int y(int pos) const { return iy(pos); }
    int x(const Site &s) const { return ix(s.index()); }
    int y(const Site &s) const { return iy(s.index()); }
    int ls() const { return size(); }
    size_t bc() const { return _bc; }
    bool snake() const { return _snake; }
    void use_snake(bool snake) { _snake = snake; }    

    // Streams

    void write(std::ostream &s) const
    {
      s.write((const char *)&_type, sizeof(size_t));
      s.write((const char *)&_nx, sizeof(size_t));
      s.write((const char *)&_ny, sizeof(size_t));
      size_t l = size();
      s.write((const char *)&l, sizeof(size_t));
      s.write((const char *)&_bc, sizeof(size_t));
    }

    void read(std::istream &s)
    {
      s.read((char *)&_type, sizeof(size_t));
      s.read((char *)&_nx, sizeof(size_t));
      s.read((char *)&_ny, sizeof(size_t));
      size_t l;
      s.read((char *)&l, sizeof(size_t));
      s.read((char *)&_bc, sizeof(size_t));

      init();
    }

};



void
Lattice::init()
{
  size_t _ls = _nx * _ny;
  clear();
  switch(_type){
    case LATTICE_1D:
      for(int is = 0; is < _ls; is++){
        std::valarray<int> nn(NN_LAST);
        nn[NN0] = is;
        nn[NN1_RIGHT] = is + 1; 
        nn[NN1_LEFT] = is - 1; 
        nn[NN1_TOP] = is; 
        nn[NN1_BOTTOM] = is; 
        nn[NN2_TOP_RIGHT] = is; 
        nn[NN2_TOP_LEFT] = is; 
        nn[NN2_BOTTOM_LEFT] = is; 
        nn[NN2_BOTTOM_RIGHT] = is; 
        nn[NN3_RIGHT] = is + 2; 
        nn[NN3_LEFT] = is - 2; 
        nn[NN3_TOP] = is; 
        nn[NN3_BOTTOM] = is; 
        push_back(Site(is,nn));
      }
      _ls = size();
      if(_bc == OBC){
        site(0)[NN1_LEFT] = 0;
        site(0)[NN3_LEFT] = 0;
        if(_nx > 1) site(1)[NN3_LEFT] = 1;
        site(_ls-1)[NN1_RIGHT] = _ls-1;
        site(_ls-1)[NN3_RIGHT] = _ls-1;
        if(_ls > 1) site(_ls-2)[NN3_RIGHT] = _ls-2;
      }else if(_bc == PBC){
        site(0)[NN1_LEFT] = _ls-1;
        site(0)[NN3_LEFT] = _ls-2;
        if(_nx > 1)site(1)[NN3_LEFT] = _ls-1;
        site(_ls-1)[NN1_RIGHT] = 0;
        site(_ls-1)[NN3_RIGHT] = 1;
        if(_ls > 1)site(_ls-2)[NN3_RIGHT] = 0;
      }
      break;
    case LATTICE_2D:
      for(int i = 0; i < _nx*_ny; i++){
          int x = ix(i); 
          int y = iy(i); 
          std::valarray<int> nn(NN_LAST);
          nn[NN0] = map2d(x,y);
          nn[NN1_RIGHT] = map2d(x+1,y); 
          nn[NN1_LEFT] = map2d(x-1,y); 
          nn[NN1_TOP] = map2d(x,y+1); 
          nn[NN1_BOTTOM] = map2d(x,y-1); 
          nn[NN2_TOP_RIGHT] = map2d(x+1,y+1); 
          nn[NN2_TOP_LEFT] = map2d(x-1,y+1); 
          nn[NN2_BOTTOM_LEFT] = map2d(x-1,y-1); 
          nn[NN2_BOTTOM_RIGHT] = map2d(x+1,y-1); 
          nn[NN3_RIGHT] = map2d(x+2,y); 
          nn[NN3_LEFT] = map2d(x-2,y); 
          nn[NN3_TOP] = map2d(x,y+2); 
          nn[NN3_BOTTOM] = map2d(x,y-2); 
          push_back(Site(x,y,nn));
      }
      for(int x = 0; x < _nx; x++){
        if(_bc == OBC){
          site(x,0)[NN1_BOTTOM] = map2d(x,0);
          site(x,0)[NN2_BOTTOM_RIGHT] = map2d(x,0);
          site(x,0)[NN2_BOTTOM_LEFT] = map2d(x,0);
          site(x,0)[NN3_BOTTOM] = map2d(x,0);
          site(x,1)[NN3_BOTTOM] = map2d(x,1);
          site(x,_ny-1)[NN1_TOP] = map2d(x,_ny-1);
          site(x,_ny-1)[NN2_TOP_RIGHT] = map2d(x,_ny-1);
          site(x,_ny-1)[NN2_TOP_LEFT] = map2d(x,_ny-1);
          site(x,_ny-1)[NN3_TOP] = map2d(x,_ny-1);
          site(x,_ny-2)[NN3_TOP] = map2d(x,_ny-2);
        }else if(_bc == PBC || _bc == CBC){
          site(x,0)[NN1_BOTTOM] = map2d(x,_ny-1);
          if(x < _nx-1){
            site(x,0)[NN2_BOTTOM_RIGHT] = map2d(x+1,_ny-1);
            site(x,_ny-1)[NN2_TOP_RIGHT] = map2d(x+1,0);
          }else{
            site(x,0)[NN2_BOTTOM_RIGHT] = map2d(0,_ny-1);
            site(x,_ny-1)[NN2_TOP_RIGHT] = map2d(0,0);
          }
          if(x > 0){
            site(x,0)[NN2_BOTTOM_LEFT] = map2d(x-1,_ny-1);
            site(x,_ny-1)[NN2_TOP_LEFT] = map2d(x-1,0);
          }else{
            site(x,0)[NN2_BOTTOM_LEFT] = map2d(_nx-1,_ny-1);
            site(x,_ny-1)[NN2_TOP_LEFT] = map2d(_nx-1,0);
          }
          site(x,0)[NN3_BOTTOM] = map2d(x,_ny-2);
          site(x,1)[NN3_BOTTOM] = map2d(x,_ny-1);
          site(x,_ny-1)[NN1_TOP] = map2d(x,0);
          site(x,_ny-1)[NN3_TOP] = map2d(x,1);
          site(x,_ny-2)[NN3_TOP] = map2d(x,0);
        }
      }
      if(_bc == OBC || _bc == CBC){
        for(int y = 0; y < _ny; y++){
          site(0,y)[NN1_LEFT] = map2d(0,y);
          site(0,y)[NN3_LEFT] = map2d(0,y);
          site(1,y)[NN3_LEFT] = map2d(1,y);
          site(0,y)[NN2_BOTTOM_LEFT] = map2d(0,y);
          site(0,y)[NN2_TOP_LEFT] = map2d(0,y);
          site(_nx-1,y)[NN1_RIGHT] = map2d(_nx-1,y);
          site(_nx-1,y)[NN3_RIGHT] = map2d(_nx-1,y);
          site(_nx-2,y)[NN3_RIGHT] = map2d(_nx-2,y);
          site(_nx-1,y)[NN2_BOTTOM_RIGHT] = map2d(_nx-1,y);
          site(_nx-1,y)[NN2_TOP_RIGHT] = map2d(_nx-1,y);
        }
      }else if(_bc == PBC){
        for(int y = 0; y < _ny; y++){
          site(0,y)[NN1_LEFT] = map2d(_nx-1,y);
          site(0,y)[NN3_LEFT] = map2d(_nx-2,y);
          site(1,y)[NN3_LEFT] = map2d(_nx-1,y);
          site(0,y)[NN2_BOTTOM_LEFT] = map2d(_nx-1,y-1);
          site(0,y)[NN2_TOP_LEFT] = map2d(_nx-1,y+1);
          site(_nx-1,y)[NN1_RIGHT] = map2d(0,y);
          site(_nx-1,y)[NN3_RIGHT] = map2d(1,y);
          site(_nx-2,y)[NN3_RIGHT] = map2d(1,y);
          site(_nx-1,y)[NN2_BOTTOM_RIGHT] = map2d(0,y-1);
          site(_nx-1,y)[NN2_TOP_RIGHT] = map2d(0,y+1);
        }
        site(0,0)[NN2_BOTTOM_LEFT] = map2d(_nx-1,_ny-1);
        site(0,_ny-1)[NN2_TOP_LEFT] = map2d(_nx-1,0);
        site(_nx-1,0)[NN2_BOTTOM_RIGHT] = map2d(0,_ny-1);
        site(_nx-1,_ny-1)[NN2_TOP_RIGHT] = map2d(0,0);
      }
    case LATTICE_NULL:
    default:
      break;
  }
}

} // namespace dmtk

#endif // __DMTK_LATTICE_H__
