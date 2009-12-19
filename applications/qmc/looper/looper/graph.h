/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
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
