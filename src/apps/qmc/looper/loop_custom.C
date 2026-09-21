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

#include "loop_config.h"
#include <looper/custom_impl.h>

//
// explicit instantiation of custom_measurement_initializer
//

template struct
looper::custom_measurement_initializer<looper::lattice_helper<loop_config::lattice_graph_t> >;
