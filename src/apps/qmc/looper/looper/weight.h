/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef LOOPER_WEIGHT_H
#define LOOPER_WEIGHT_H

namespace looper {

//
// default graph weights for path-integral and SSE loop algorithms
//

struct site_weight_helper;

struct xxz_bond_weight_helper;

struct xyz_bond_weight_helper;

template<typename SITE_WEIGHT_HELPER = site_weight_helper,
         typename BOND_WEIGHT_HELPER = xyz_bond_weight_helper>
class weight_helper;

} // namespace looper

#endif // LOOPER_WEIGHT_H
