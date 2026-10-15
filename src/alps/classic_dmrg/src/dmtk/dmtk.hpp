// Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/dmtk.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#ifndef __DMTK_H__
#define __DMTK_H__

namespace dmtk
{

enum
{
  LEFT,
  RIGHT,
};

} //namespace dmtk

#include <cmath>
#include <valarray>
#include "dmtk/conj.hpp"
#include "dmtk/support.hpp"
#include "dmtk/bits.hpp"
#include "dmtk/enums.hpp"
#include "dmtk/basis.hpp"
#include "dmtk/vector.hpp"
#include "dmtk/matrix.hpp"
#include "dmtk/qn.hpp"
#include "dmtk/subspace.hpp"
#include "dmtk/operators.hpp"
#include "dmtk/block_matrix.hpp"
#include "dmtk/block.hpp"
#include "dmtk/state.hpp"
#include "dmtk/state_slice.hpp"
#include "dmtk/lattice.hpp"
#include "dmtk/system.hpp"
#include "dmtk/util.hpp"
#include "dmtk/hami.hpp"

#endif // __DMTK_H__
