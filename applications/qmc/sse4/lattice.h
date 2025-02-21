/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2010 by Sergei Isakov <isakov@itp.phys.ethz.ch>
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

#ifndef __LATTICE_H__
#define __LATTICE_H__

#include <map>
#include <set>
#include <vector>
#include <boost/tuple/tuple.hpp>

#include <alps/lattice.h>

const unsigned UNIT_SIZE = 2;

template <typename L>
class Lattice {
public:
    typedef unsigned lat_unit_sites_type[UNIT_SIZE];
    typedef unsigned bond_sites_type[2];
    typedef unsigned lat_unit_bonds_type[1];
private:
    template<unsigned UNIT_SIZE>
    struct Lat_unit {
        unsigned type;
        lat_unit_sites_type sites;
    };    
public:
    typedef L graph_type;

    typedef Lat_unit<UNIT_SIZE> lat_unit_type;
    typedef std::vector<lat_unit_type> lat_units_type;
    
    typedef typename graph_type::sites_size_type sites_size_type;
    typedef typename graph_type::site_descriptor site_descriptor;
    typedef typename graph_type::bonds_size_type bonds_size_type;
    typedef typename graph_type::bond_descriptor bond_descriptor;
    typedef typename graph_type::vector_type vector_type;
    
    typedef boost::tuple<unsigned, unsigned, unsigned,
                unsigned, unsigned> bond_tuple_type;
    typedef std::map<bond_tuple_type, unsigned> bond_types_type;
        
    Lattice(alps::Parameters const& params, graph_type& graph) :
        graph(graph),
        params(params)
    {
        _max_site_type = 0;
        _lat_units.resize(graph.num_bonds());
        for (bonds_size_type i = 0; i < graph.num_bonds(); ++i) {
            bond_descriptor bond = graph.bond(i);
            site_descriptor source = graph.source(bond);
            site_descriptor target = graph.target(bond);
            
            bond_tuple_type bond_tuple(
                graph.inhomogeneous_bond_type(bond),
                graph.inhomogeneous_site_type(source),
                graph.inhomogeneous_site_type(target),
                graph.num_neighbors(source), graph.num_neighbors(target));
            
            typename bond_types_type::const_iterator ite = bond_types.find(bond_tuple);
            if (ite != bond_types.end())
                _lat_units[i].type = ite->second;
            else {
                bond_type2bondi.push_back(i);
                
                _lat_units[i].type = bond_types.size();
                bond_types[bond_tuple] = _lat_units[i].type;
            }
            
            _lat_units[i].sites[0] = source;
            _lat_units[i].sites[1] = target;
            
            unsigned stype1 = sitei2alps_type(source);
            unsigned stype2 = sitei2alps_type(target);
            
            if (stype1 > _max_site_type) _max_site_type = stype1;
            if (stype2 > _max_site_type) _max_site_type = stype2;
        }
    }
    
    unsigned dimension() const
    {
        return graph.dimension();
    }
    
    unsigned nsites() const
    {
        return graph.num_sites();
    }
    
    unsigned nbonds() const
    {
        return graph.num_bonds();
    }
    
    unsigned ndistances() const
    {
        return graph.num_distances();
    }
    
    unsigned distance(site_descriptor x, site_descriptor y) const
    {
        return graph.distance(x, y);
    }
    
    bool inhomogeneous() const
    {
        return graph.inhomogeneous();
    }
    
    unsigned nlat_unit_types() const
    {
        return bond_types.size();
    }
    
    unsigned max_site_type()
    {
        return _max_site_type;
    }
    
    lat_units_type const& lat_units() const
    {
        return _lat_units;
    }
    
    unsigned sitei2alps_type(sites_size_type i) const
    {
        // fix me: we assume that graph.site(i) == i
        return graph.site_type(i);
    }
    
    unsigned bondi2alps_type(bonds_size_type i) const
    {
        return graph.bond_type(graph.bond(i));
    }
    
    template<typename P>
    void check_homogenity(P const& params, bonds_size_type i)
    {
        if (graph.inhomogeneous())
            graph.throw_if_xyz_defined(params, graph.bond(i));
    }
    
    template<typename P>
    void do_inhomogeneous_source(P& params, bonds_size_type i)
    {
        if (graph.inhomogeneous_sites())
            params << graph.coordinate_as_parameter(graph.source(graph.bond(i)));
    }
    
    template<typename P>
    void do_inhomogeneous_target(P& params, bonds_size_type i)
    {
        if (graph.inhomogeneous_sites())
            params << graph.coordinate_as_parameter(graph.target(graph.bond(i)));
    }
    
    template<typename P>
    void do_inhomogeneous_bond(P& params, bonds_size_type i)
    {
        if (graph.inhomogeneous_bonds())
            params << graph.coordinate_as_parameter(graph.bond(i));
    }
    
    void bondi2sitei(bonds_size_type i, bond_sites_type& sites) const
    {
        bond_descriptor bond = graph.bond(i);
        
        sites[0] = graph.source(bond);
        sites[1] = graph.target(bond);
    }
    
    void lat_unit_type2sitei(unsigned type, lat_unit_sites_type& sites) const
    {
        unsigned i = bond_type2bondi[type];
        
        sites[0] = _lat_units[i].sites[0];
        sites[1] = _lat_units[i].sites[1];
    }
    
    void lat_unit_type2bondi(unsigned type, lat_unit_bonds_type& bonds) const
    {
        bonds[0] = bond_type2bondi[type];
    }
    
    unsigned nneighbors(unsigned i) const
    {
        return graph.num_neighbors(i);
    }
    
    vector_type const& bond_vector_relative(bonds_size_type i) const
    {
        return graph.bond_vector_relative(graph.bond(i));
    }
private:
    lat_units_type _lat_units;
        
    graph_type& graph;
    alps::Parameters const& params;
    
    bond_types_type bond_types;
    unsigned _max_site_type;

    std::vector<unsigned> bond_type2bondi;
};

#endif
