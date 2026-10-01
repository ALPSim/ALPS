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

#ifndef LOOPER_LOCATION_H
#define LOOPER_LOCATION_H

namespace looper {

class location;

// optimized version for models with bond terms only
class location_bond;

// for long-range models
class location_longrange;

} // end namespace looper

#endif // LOOPER_LOCATION_H
