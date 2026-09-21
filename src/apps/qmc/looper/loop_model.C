/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "loop_config.h"
#include <looper/model_impl.h>

//
// explicit instantiation of spinmodel_helper::init() member function
//

template
void looper::spinmodel_helper<loop_config::lattice_graph_t, loop_config::loop_graph_t>::
  init(alps::Parameters const& p, looper::lattice_helper<loop_config::lattice_graph_t>& lat,
  bool is_path_integral);
