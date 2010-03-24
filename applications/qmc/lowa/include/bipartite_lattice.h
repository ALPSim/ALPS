/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Applications
*
* Copyright (C) 2006-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Lode Pollet <lpollet@physics.harvard.edu>
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>
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

/* $Id: lattice.h 3520 2010-03-03 16:55:00Z tamama $ */


/*
 *
 * 1) Code modification      -- done
 * 2) Replacing raw pointers -- not done yet
 *
 */


#ifndef ALPS_APPLICATION_LOWA_LATTICE_H
#define ALPS_APPLICATION_LOWA_LATTICE_H


#include <iostream>
#include <cmath>

namespace alps {
namespace applications {
namespace lowa {


template<class S, class R>
class bipartite_lattice
{
public:
  // typedefs
  typedef uint8_t dim_type;
  typedef S       site_type;
  typedef S       index_type;
  typedef R       real_coordinate_type;

  // get_functions
  inline dim_type  dim() const  { return _dim; }
  inline site_type Nx()  const  { return _Nx;  }
  inline site_type Ny()  const  { return _Ny;  }
  inline site_type Nz()  const  { return _Nz;  }
  inline site_type Nxy() const  { return _Nxy; }
  inline site_type N()   const  { return _N;   }

  inline site_type i(index_type index) const { return _xintcoord[index]; }
  inline site_type j(index_type index) const { return _yintcoord[index]; }
  inline site_type k(index_type index) const { return _zintcoord[index]; }

  inline real_coordinate_type distsq(index_type index, dim_type which_dimension)
  {
    if (which_dimension == 0)   {  return _xdistsq[index];  }
    if (which_dimension == 1)   {  return _ydistsq[index];  }
    if (which_dimension == 2)   {  return _zdistsq[index];  }
#ifdef DEBUGMODE
    std::cout << "\nCheck <which dimension> in <inline real_coordinate_type distsq(index_type index, int which_dimension)>...\n";
#endif
  }

  inline real_coordinate_type mid(dim_type which_dimension)
  {
    if (which_dimension == 0)   {  return _xmid;  }
    if (which_dimension == 1)   {  return _ymid;  }
    if (which_dimension == 2)   {  return _zmid;  }
  } 

  inline site_type nb(index_type index, dim_type which_direction)
  {
    if (_dim == 3) 
    {
      if (which_direction == 0)   {  return _east[index];    }
      if (which_direction == 1)   {  return _south[index];   }
      if (which_direction == 2)   {  return _bottom[index];  }
      if (which_direction == 3)   {  return _west[index];    }
      if (which_direction == 4)   {  return _north[index];   }
      if (which_direction == 5)   {  return _top[index];     }
    }
    else if (_dim == 2)
    {
      if (which_direction == 0)   {  return _east[index];    }
      if (which_direction == 1)   {  return _south[index];   }
      if (which_direction == 2)   {  return _west[index];    }
      if (which_direction == 3)   {  return _north[index];   }
    }
    else if (_dim == 1)
    {
      if (which_direction == 0)   {  return _east[index];    }
      if (which_direction == 1)   {  return _west[index];    }
    }
    std::cout << "Something is wrong with <inline site_type nb(index_type index, dim_type which_direction)> of <bipartite_lattice.h>\n";
  }

  inline bool is_border_site(index_type index)  {  return _is_border_site[index];  }


private:
  dim_type _dim;

  site_type _Nx;
  site_type _Ny;
  site_type _Nz;
  site_type _Nxy;
  site_type _N;

  site_type* _xintcoord;
  site_type* _yintcoord;
  site_type* _zintcoord;

  real_coordinate_type _xmid;
  real_coordinate_type _ymid;
  real_coordinate_type _zmid;

  real_coordinate_type* _xdistsq;
  real_coordinate_type* _ydistsq;
  real_coordinate_type* _zdistsq;

  site_type* _east;
  site_type* _west;
  site_type* _south;
  site_type* _north;
  site_type* _bottom;
  site_type* _top;

  bool*  _is_border_site;

public:
  // destructor
  ~bipartite_lattice()
  {
    if (_dim == 3)
    {
      delete [] _xintcoord;
      delete [] _yintcoord;
      delete [] _zintcoord;

      delete [] _xdistsq;
      delete [] _ydistsq;
      delete [] _zdistsq;

      delete [] _east;
      delete [] _west;
      delete [] _south;
      delete [] _north;
      delete [] _bottom;
      delete [] _top;
    }
    else if (_dim == 2)
    {
      delete [] _xintcoord;
      delete [] _yintcoord;

      delete [] _xdistsq;
      delete [] _ydistsq;

      delete [] _east;
      delete [] _west;
      delete [] _south;
      delete [] _north;
    }
    else if (_dim == 1)
    {
      delete [] _xintcoord;

      delete [] _xdistsq;

      delete [] _east;
      delete [] _west;
    }
    delete [] _is_border_site;
  }

  // constructors
  bipartite_lattice()  {}
  bipartite_lattice(site_type const Nx)  { init(Nx); }
  bipartite_lattice(site_type const Nx, site_type const Ny)  { init(Nx,Ny); }
  bipartite_lattice(site_type const Nx, site_type const Ny, site_type const Nz)  { init(Nx,Ny,Nz); }


  void init(site_type const Nx)
  {
    _dim = 1;

    _xmid = 0.5*(Nx-1);

    _Nx = Nx;

    _N = Nx;

    _xintcoord = new site_type [Nx];

    _xdistsq = new real_coordinate_type [Nx];

    _east   = new site_type [Nx];
    _west   = new site_type [Nx];

    _is_border_site = new bool [Nx];

    for (index_type index=0, xindex=0; xindex < Nx; ++xindex, ++index) {
      _xintcoord[index] = xindex;

      _xdistsq[index] = (xindex - _xmid) * (xindex - _xmid);

      _east[index]   =  ((xindex == (_Nx-1)) ? (index+1   -_Nx)  : (index+1));
      _west[index]   =  ((xindex == 0) ? (index-1   +_Nx)  : (index-1));

      _is_border_site[index] = false;
      if (xindex == 0)       _is_border_site[index] = true;
      if (xindex == (_Nx-1)) _is_border_site[index] = true;
    }
  }


  void init(site_type const Nx, site_type const Ny)
  {
    _dim = 2;

    _xmid = 0.5*(Nx-1);
    _ymid = 0.5*(Ny-1);

    _Nx = Nx;
    _Ny = Ny;

    _Nxy = Nx*Ny;
    _N   = Nx*Ny;

    _xintcoord = new site_type [Nx*Ny];
    _yintcoord = new site_type [Nx*Ny];

    _xdistsq = new real_coordinate_type [Nx*Ny];
    _ydistsq = new real_coordinate_type [Nx*Ny];

    _east   = new site_type [Nx*Ny];
    _west   = new site_type [Nx*Ny];
    _south  = new site_type [Nx*Ny];
    _north  = new site_type [Nx*Ny];

    _is_border_site = new bool [Nx*Ny];

    for (index_type index=0, yindex=0; yindex < Ny; ++yindex) {
      for (index_type xindex=0; xindex < Nx; ++xindex, ++index) {
        _xintcoord[index] = xindex;
        _yintcoord[index] = yindex;

        _xdistsq[index] = (xindex - _xmid) * (xindex - _xmid);
        _ydistsq[index] = (yindex - _ymid) * (yindex - _ymid);

        _east[index]   =  ((xindex == (_Nx-1)) ? (index+1   -_Nx)  : (index+1));
        _south[index]  =  ((yindex == (_Ny-1)) ? (index+_Nx -_Nxy) : (index+_Nx));
        _west[index]   =  ((xindex == 0) ? (index-1   +_Nx)  : (index-1));
        _north[index]  =  ((yindex == 0) ? (index-_Nx +_Nxy) : (index-_Nx));

        _is_border_site[index] = false;
        if (xindex == 0)       _is_border_site[index] = true;
        if (xindex == (_Nx-1)) _is_border_site[index] = true;
        if (yindex == 0)       _is_border_site[index] = true;
        if (yindex == (_Ny-1)) _is_border_site[index] = true;
      }
    }
  }


  void init(site_type const Nx, site_type const Ny, site_type const Nz)
  {
    _dim = 3;

    _xmid = 0.5*(Nx-1);
    _ymid = 0.5*(Ny-1);
    _zmid = 0.5*(Nz-1);

    _Nx = Nx;
    _Ny = Ny;
    _Nz = Nz;

    _Nxy = Nx*Ny;
    _N   = Nx*Ny*Nz;

    _xintcoord = new site_type [Nx*Ny*Nz];
    _yintcoord = new site_type [Nx*Ny*Nz];
    _zintcoord = new site_type [Nx*Ny*Nz];

    _xdistsq = new real_coordinate_type [Nx*Ny*Nz];
    _ydistsq = new real_coordinate_type [Nx*Ny*Nz];
    _zdistsq = new real_coordinate_type [Nx*Ny*Nz];

    _east   = new site_type [Nx*Ny*Nz];
    _west   = new site_type [Nx*Ny*Nz];
    _south  = new site_type [Nx*Ny*Nz];
    _north  = new site_type [Nx*Ny*Nz];
    _bottom = new site_type [Nx*Ny*Nz];
    _top    = new site_type [Nx*Ny*Nz];
    
    _is_border_site = new bool [Nx*Ny*Nz];

    for (index_type index=0, zindex=0; zindex < Nz; ++zindex) {
      for (index_type yindex=0; yindex < Ny; ++yindex) {
        for (index_type xindex=0; xindex < Nx; ++xindex, ++index) {
          _xintcoord[index] = xindex;
          _yintcoord[index] = yindex;
          _zintcoord[index] = zindex;

          _xdistsq[index] = (xindex - _xmid) * (xindex - _xmid);
          _ydistsq[index] = (yindex - _ymid) * (yindex - _ymid);
          _zdistsq[index] = (zindex - _zmid) * (zindex - _zmid);

          _east[index]   =  ((xindex == (_Nx-1)) ? (index+1   -_Nx)  : (index+1));
          _south[index]  =  ((yindex == (_Ny-1)) ? (index+_Nx -_Nxy) : (index+_Nx));
          _bottom[index] =  ((zindex == (_Nz-1)) ? (index+_Nxy-_N)   : (index+_Nxy));
          _west[index]   =  ((xindex == 0) ? (index-1   +_Nx)  : (index-1));
          _north[index]  =  ((yindex == 0) ? (index-_Nx +_Nxy) : (index-_Nx));
          _top[index]    =  ((zindex == 0) ? (index-_Nxy+_N)   : (index-_Nxy));

          _is_border_site[index] = false;
          if (xindex == 0)       _is_border_site[index] = true;
          if (xindex == (_Nx-1)) _is_border_site[index] = true;
          if (yindex == 0)       _is_border_site[index] = true;
          if (yindex == (_Ny-1)) _is_border_site[index] = true;
          if (zindex == 0)       _is_border_site[index] = true;
          if (zindex == (_Nz-1)) _is_border_site[index] = true;
        }
      }
    }
  }


/* Lode's comment as follows:
 *  -LATTICE requirements.
 *  -Different for open boundary condition (OBC) and periodic boundary conditions (PBC)
 *  -returns the site that behaves to 0 as s2 behaves to s1
 *
 * Tama's comment as follows:
 *  -very frequently used function, modified to speed up calculations
 */


  inline site_type dist(site_type const s1, site_type const s2)
  {
    if (_dim == 3) 
    {
#ifndef TRAPPEDSYSTEM
      site_type i0 = _xintcoord[s2] - _xintcoord[s1];
      site_type j0 = _yintcoord[s2] - _yintcoord[s1];
      site_type k0 = _zintcoord[s2] - _zintcoord[s1];

      if (i0 >= 0)  {  i0 += _Nx; }
      if (j0 >= 0)  {  j0 += _Ny; }
      if (k0 >= 0)  {  k0 += _Nz; }

      site_type s0 = i0 + j0*_Nx + k0*_Nxy;
      return s0;
#else
      site_type i0 = (_xintcoord[s2] - _xintcoord[s1] > 0 ? _xintcoord[s2] - _xintcoord[s1] : -(_xintcoord[s2] - _xintcoord[s1]));
      site_type j0 = (_yintcoord[s2] - _yintcoord[s1] > 0 ? _yintcoord[s2] - _yintcoord[s1] : -(_yintcoord[s2] - _yintcoord[s1]));
      site_type k0 = (_zintcoord[s2] - _zintcoord[s1] > 0 ? _zintcoord[s2] - _zintcoord[s1] : -(_zintcoord[s2] - _zintcoord[s1]));

      site_type s0 = i0 + j0*_Nx + k0*_Nxy;
      return s0;
#endif
    }




  }

};


} // ending namespace lowa
} // ending namespace applications
} // ending namespace alps


#endif
