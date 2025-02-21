/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_ALPS_LATTICE_H
#define APP_ALPS_LATTICE_H

#include "dmrg/models/lattice.h"

#include <vector>
#include <map>
#include <algorithm>
#include <sstream>

#include <alps/parameter.h>
#include <alps/lattice.h>


class alps_lattice : public lattice_impl
{

public:
    typedef lattice_impl::pos_t pos_t;
    typedef alps::graph_helper<> graph_type;
    typedef graph_type::site_descriptor site_descriptor;
    typedef graph_type::site_iterator site_iterator;
    
    alps_lattice (const alps::Parameters& p)
    : parms(p)
    , graph(parms)
    {
        // storing lattice informations
        forward_.resize(size());
        backward_.resize(size());
                
        for (graph_type::bond_iterator it=graph.bonds().first; it!=graph.bonds().second; ++it) {
            graph_type::size_type s, t;
            s = graph.vertex_index(graph.source(*it));
            t = graph.vertex_index(graph.target(*it));
            
            forward_[s].push_back(t);
            backward_[t].push_back(s);
            
            bond_index_map[s][t] = graph.edge_index(*it);
            bond_index_map[t][s] = graph.edge_index(*it);
        }
    }
    
    virtual std::vector<pos_t> forward(pos_t p) const
    {
        return forward_[p];
    }
    virtual std::vector<pos_t> all(pos_t p) const
    {
        std::vector<pos_t> ret = forward_[p];
        std::copy(backward_[p].begin(), backward_[p].end(), std::back_inserter(ret));
        return ret; 
    }
    
    virtual pos_t size() const
    {
        return graph.num_sites();
    }
    
    virtual int maximum_vertex_type() const
    {
        return alps::maximum_vertex_type(graph.graph());
    }
    
    virtual boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "label" && pos.size() == 1)
            return boost::any( alps::site_label(graph.graph(), graph.site(pos[0])) );
        else if (property == "label" && pos.size() == 2)
            return boost::any( alps::bond_label(graph.graph(), graph.bond(bond_index_map[pos[0]][pos[1]])) );
        else if (property == "type" && pos.size() == 1)
            return boost::any( static_cast<int>(graph.site_type(graph.site(pos[0]))) );
        else if (property == "type" && pos.size() == 2)
            return boost::any( static_cast<int>(graph.bond_type(graph.bond(bond_index_map[pos[0]][pos[1]]))) );
        else if (property == "wraps_pbc" && pos.size() == 2)
            return boost::any( static_cast<bool>(boost::get(alps::boundary_crossing_t(),
                                                            graph.graph(),
                                                            graph.bond(bond_index_map[pos[0]][pos[1]]))) );
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented."; 
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }

    virtual const graph_type& alps_graph() const
    {
        return graph;
    }
    
protected:
    alps::Parameters parms;
    graph_type graph;
    mutable std::vector<std::vector<pos_t> > forward_;
    mutable std::vector<std::vector<pos_t> > backward_;
    // TODO: find better solution instead of this (operator[] is non-const)
    mutable std::map<pos_t, std::map<pos_t, pos_t> > bond_index_map;
};


#endif
