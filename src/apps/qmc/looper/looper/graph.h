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

#ifndef LOOPER_GRAPH_H
#define LOOPER_GRAPH_H

#include "location.h"

namespace looper {

//
// site graph type
//

struct site_graph_type;

//
// bond graph types
//

// optimized bond_graph_type for Ising model
struct ising_bond_graph_type;

// optimized bond_graph_type for Heisenberg Antiferromagnet
struct haf_bond_graph_type;

// optimized bond_graph_type for Heisenberg Ferromagnet
struct hf_bond_graph_type;

// bond_graph_type for XXZ interaction
struct xxz_bond_graph_type;

// bond_graph_type for XYZ interaction
struct xyz_bond_graph_type;

template<typename SITE = site_graph_type, typename BOND = xxz_bond_graph_type,
  typename LOC = location>
class local_graph;

} // end namespace looper

#endif // LOOPER_GARPH_H
