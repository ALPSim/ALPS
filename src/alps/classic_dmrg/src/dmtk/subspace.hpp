// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/subspace.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef __DMTK_SUBSPACE_H__
#define __DMTK_SUBSPACE_H__

#include <iostream>
#include <iosfwd>
#include "dmtk/qn.hpp"
#include "dmtk/range.hpp"

namespace dmtk
{

class SubSpace: public Range
{
  private:
    QN _qn;
  public:

    SubSpace(): Range(0,0), _qn(0) {};
    SubSpace(const QN &qn, int _begin, int _end): 
      Range(_begin,_end), _qn(qn)  {}
    SubSpace(const Range &s): Range(s) {}

    size_t dim() const { return std::slice::size(); }
    QN& qn() { return _qn; }
    QN qn() const { return _qn; }

    // Streams

    void write(std::ostream &s) const
    {
      _qn.write(s);
      int _begin = begin(), _end = end();
      s.write(reinterpret_cast<const char*>(&_begin), sizeof(int));
      s.write(reinterpret_cast<const char*>(&_end), sizeof(int));
    }

    void read(std::istream &s)
    {
      QN qn;
      int _begin, _end;

      qn.read(s);
      s.read(reinterpret_cast<char*>(&_begin), sizeof(int));
      s.read(reinterpret_cast<char*>(&_end), sizeof(int));

      *this = SubSpace(qn, _begin, _end);
    }

};

} // namespace dmtk

#endif // __DMTK_SUBSPACE_H__
