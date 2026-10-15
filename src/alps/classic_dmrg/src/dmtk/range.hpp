// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/range.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef __RANGE_H__
#define __RANGE_H__

namespace dmtk
{

class Range: public std::slice
{
  public:
    Range(size_t _o, size_t _d, size_t _s): std::slice(_o, (_d-_o+1)/_s, _s) {}
    Range(size_t _o, size_t _d): std::slice(_o, (_d-_o+1), 1) {}
    Range(const Range &r): std::slice(r) {}

    bool operator==(const Range &r)
     { 
       if(begin() == r.begin() && end() == r.end() && stride() == r.stride())
         return true;
       return false;
     }

    size_t begin() const { return std::slice::start(); }
    size_t end() const { return std::slice::start()+std::slice::size()-1; }
};

} // namespace dmtk

#endif // __RANGE_H__
