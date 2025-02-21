/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_MATRIX_H
#define LOOPER_MATRIX_H

#include "power.h"
#include <alps/model/half_integer.h>
#include <boost/multi_array.hpp>


// IBM workaround
#ifdef hz
#undef hz
#endif

namespace looper {

//
// site matrix
//

class site_matrix {
public:
  typedef boost::multi_array<double, 2> matrix_type;
  typedef matrix_type::size_type size_type;

  site_matrix() : mat_() {}
  site_matrix(site_matrix const& m) : mat_(m.mat_) {}
  template<typename SP>
  site_matrix(SP const& sp) : mat_() { build(sp); }

  double operator()(size_type i, size_type j) const { return mat_[i][j]; }
  matrix_type const& matrix() const { return mat_; }

protected:
  template<typename SP>
  void build(SP const& sp);

private:
  boost::multi_array<double, 2> mat_;
};

//
// bond matrices
//

class bond_matrix_xxz {
public:
  typedef boost::multi_array<double, 4> matrix_type;
  typedef matrix_type::size_type size_type;

  bond_matrix_xxz() : mat_() {}
  bond_matrix_xxz(const bond_matrix_xxz& m) : mat_(m.mat_) {}
  template<typename I, typename BP>
  bond_matrix_xxz(alps::half_integer<I> const& s0, alps::half_integer<I> const& s1, BP const& bp) :
    mat_() {
    build(s0, s1, bp);
  }
  template<typename SP, typename BP>
  bond_matrix_xxz(SP const& sp0, SP const& sp1, BP const& bp) : mat_() {
    build(sp0.s, sp1.s, bp);
  }

  double operator()(size_type i, size_type j, size_type k, size_type l) const {
    return mat_[i][j][k][l];
  }
  const matrix_type& matrix() const { return mat_; }

protected:
  template<typename I, typename BP>
  void build(alps::half_integer<I> const& s0, alps::half_integer<I> const& s1, BP const& bp);

private:
  boost::multi_array<double, 4> mat_;
};

class bond_matrix_xyz {
public:
  typedef boost::multi_array<double, 4> matrix_type;
  typedef matrix_type::size_type size_type;

  bond_matrix_xyz() : mat_() {}
  bond_matrix_xyz(const bond_matrix_xyz& m) : mat_(m.mat_) {}
  template<typename I, typename BP>
  bond_matrix_xyz(alps::half_integer<I> const& s0, alps::half_integer<I> const& s1, BP const& bp) :
    mat_() {
    build(s0, s1, bp);
  }
  template<typename SP, typename BP>
  bond_matrix_xyz(SP const& sp0, SP const& sp1, BP const& bp) : mat_() {
    build(sp0.s, sp1.s, bp);
  }

  double operator()(size_type i, size_type j, size_type k, size_type l) const {
    return mat_[i][j][k][l];
  }
  const matrix_type& matrix() const { return mat_; }

protected:
  template<typename I, typename BP>
  void build(alps::half_integer<I> const& s0, alps::half_integer<I> const& s1, BP const& bp);

private:
  boost::multi_array<double, 4> mat_;
};

//
// implementation of site_matrix
//

template<typename SP>
void site_matrix::build(SP const& sp) {
  using std::sqrt;
  typedef typename SP::spin_type spin_type;

  spin_type s = sp.s;

  // set matrix dimension
  int dim = s.get_twice()+1;
  mat_.resize(boost::extents[dim][dim]);
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      mat_[i][j] = 0;

  // diagonal elements: c - hz sz + d sz^2
  for (spin_type sz = -s; sz <= s; ++sz)
    mat_[sz.distance(-s)][sz.distance(-s)] =
      sp.c - sp.hz * to_double(sz) +  sp.d * power2(to_double(sz));

  // off-diagonal elements: - hx s+ / 2
  for (spin_type sz = -s; sz <= s-1; ++sz)
    mat_[sz.distance(-s)+1][sz.distance(-s)] =
      - 0.5 * sp.hx * sqrt(to_double(s-sz) * to_double(s+sz+1));

  // off-diagonal elements: - hx s- / 2
  for (spin_type sz = -s+1; sz <= s; ++sz)
    mat_[sz.distance(-s)-1][sz.distance(-s)] =
      - 0.5 * sp.hx * sqrt(to_double(s+sz) * to_double(s-sz+1));
}

//
// implementation of bond matxices
//

template<typename I, typename BP>
void bond_matrix_xxz::build(alps::half_integer<I> const& s0, alps::half_integer<I> const& s1,
  BP const& bp) {
  typedef alps::half_integer<I> spin_type;
  using std::sqrt;

  // set matrix dimension
  int d0 = s0.get_twice()+1;
  int d1 = s1.get_twice()+1;
  mat_.resize(boost::extents[d0][d1][d0][d1]);
  for (int i0 = 0; i0 < d0; ++i0)
    for (int i1 = 0; i1 < d1; ++i1)
      for (int j0 = 0; j0 < d0; ++j0)
        for (int j1 = 0; j1 < d1; ++j1)
          mat_[i0][i1][j0][j1] = 0;

  // diagonal elements: c + jz sz0 sz1
  for (spin_type sz0 = -s0; sz0 <= s0; ++sz0) {
    for (spin_type sz1 = -s1; sz1 <= s1; ++sz1) {
      mat_[sz0.distance(-s0)][sz1.distance(-s1)]
        [sz0.distance(-s0)][sz1.distance(-s1)] =
        bp.c + bp.jz * to_double(sz0) * to_double(sz1);
    }
  }

  // off-diagonal elements: jxy s0+ s1- / 2
  for (spin_type sz0 = -s0; sz0 <= s0-1; ++sz0) {
    for (spin_type sz1 = -s1+1; sz1 <= s1; ++sz1) {
      mat_[sz0.distance(-s0)+1][sz1.distance(-s1)-1]
        [sz0.distance(-s0)][sz1.distance(-s1)] =
        0.5 * bp.jxy *
        sqrt(to_double(s0-sz0) * to_double(s0+sz0+1)) *
        sqrt(to_double(s1+sz1) * to_double(s1-sz1+1));
    }
  }

  // off-diagonal elements: jxy s0- s1+ / 2
  for (spin_type sz0 = -s0+1; sz0 <= s0; ++sz0) {
    for (spin_type sz1 = -s1; sz1 <= s1-1; ++sz1) {
      mat_[sz0.distance(-s0)-1][sz1.distance(-s1)+1]
        [sz0.distance(-s0)][sz1.distance(-s1)] =
        0.5 * bp.jxy *
        sqrt(to_double(s0+sz0) * to_double(s0-sz0+1)) *
        sqrt(to_double(s1-sz1) * to_double(s1+sz1+1));
    }
  }
}

template<typename I, typename BP>
void bond_matrix_xyz::build(alps::half_integer<I> const& s0, alps::half_integer<I> const& s1,
  BP const& bp) {
  typedef alps::half_integer<I> spin_type;
  using std::sqrt;

  // set matrix dimension
  int d0 = s0.get_twice()+1;
  int d1 = s1.get_twice()+1;
  mat_.resize(boost::extents[d0][d1][d0][d1]);
  for (int i0 = 0; i0 < d0; ++i0)
    for (int i1 = 0; i1 < d1; ++i1)
      for (int j0 = 0; j0 < d0; ++j0)
        for (int j1 = 0; j1 < d1; ++j1)
          mat_[i0][i1][j0][j1] = 0;

  // diagonal elements: c + jz sz0 sz1
  for (spin_type sz0 = -s0; sz0 <= s0; ++sz0) {
    for (spin_type sz1 = -s1; sz1 <= s1; ++sz1) {
      mat_[sz0.distance(-s0)][sz1.distance(-s1)]
        [sz0.distance(-s0)][sz1.distance(-s1)] =
        bp.c + bp.jz * to_double(sz0) * to_double(sz1);
    }
  }

  // off-diagonal elements: (jx+jy) s0+ s1- / 4
  for (spin_type sz0 = -s0; sz0 <= s0-1; ++sz0) {
    for (spin_type sz1 = -s1+1; sz1 <= s1; ++sz1) {
      mat_[sz0.distance(-s0)+1][sz1.distance(-s1)-1]
        [sz0.distance(-s0)][sz1.distance(-s1)] =
        0.25 * (bp.jx + bp.jy) *
        sqrt(to_double(s0-sz0) * to_double(s0+sz0+1)) *
        sqrt(to_double(s1+sz1) * to_double(s1-sz1+1));
    }
  }

  // off-diagonal elements: (jx+jy) s0- s1+ / 4
  for (spin_type sz0 = -s0+1; sz0 <= s0; ++sz0) {
    for (spin_type sz1 = -s1; sz1 <= s1-1; ++sz1) {
      mat_[sz0.distance(-s0)-1][sz1.distance(-s1)+1]
        [sz0.distance(-s0)][sz1.distance(-s1)] =
        0.25 * (bp.jx + bp.jy) *
        sqrt(to_double(s0+sz0) * to_double(s0-sz0+1)) *
        sqrt(to_double(s1-sz1) * to_double(s1+sz1+1));
    }
  }

  // off-diagonal elements: (jx-jy) s0+ s1+ / 4
  for (spin_type sz0 = -s0; sz0 <= s0-1; ++sz0) {
    for (spin_type sz1 = -s1; sz1 <= s1-1; ++sz1) {
      mat_[sz0.distance(-s0)+1][sz1.distance(-s1)+1]
        [sz0.distance(-s0)][sz1.distance(-s1)] =
        0.25 * (bp.jx - bp.jy) *
        sqrt(to_double(s0-sz0) * to_double(s0+sz0+1)) *
        sqrt(to_double(s1-sz1) * to_double(s1+sz1+1));
    }
  }

  // off-diagonal elements: (jx-jy) s0- s1- / 4
  for (spin_type sz0 = -s0+1; sz0 <= s0; ++sz0) {
    for (spin_type sz1 = -s1+1; sz1 <= s1; ++sz1) {
      mat_[sz0.distance(-s0)-1][sz1.distance(-s1)-1]
        [sz0.distance(-s0)][sz1.distance(-s1)] =
        0.25 * (bp.jx - bp.jy) *
        sqrt(to_double(s0+sz0) * to_double(s0-sz0+1)) *
        sqrt(to_double(s1+sz1) * to_double(s1-sz1+1));
    }
  }
}

} // end namespace looper

#endif // LOOPER_MATRIX_H
