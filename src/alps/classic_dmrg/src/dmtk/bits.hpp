// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/bits.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef __DMTK_BITS_H__
#define __DMTK_BITS_H__

// Only IS_EVEN, SGN, IBITS and mask() have live call sites; IEOR, IAND,
// IOR, IBSET, IBCLR, ISHFTC and mod2() had zero call sites.
namespace dmtk
{
#define IS_EVEN(n) (((n) & 1) == 0)
#define SGN(n) (1 - (((n) & 1) << 1))
#define IBITS(n,i) (((n) & (1 << (i))) >> (i))

inline size_t mask(int b1, int b2) { return (1 << (b1-1)) | (1 << (b2-1)); }
inline size_t mask(int b) { return (1 << (b-1)); }
} // namespace dmtk

#endif // __DMTK_BITS_H__
