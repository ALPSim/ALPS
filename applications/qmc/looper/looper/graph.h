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
