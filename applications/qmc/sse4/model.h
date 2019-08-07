/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2010 by Sergei Isakov <isakov@itp.phys.ethz.ch>
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

#ifndef __MODEL_H__
#define __MODEL_H__

#include <set>
#include <cmath>
#include <vector>
#include <limits>

#include <alps/lattice.h>
#include <alps/model.h>

#include "lattice.h"

template <typename M, typename L>
class Model {
private:
    typedef unsigned state_type;
    typedef state_type Vertex_state[2 * UNIT_SIZE];
    typedef state_type Lat_unit_state[UNIT_SIZE];

    struct Vertex {
        double me;
        bool diagonal;
        unsigned unit_type;
        unsigned num_state;        // numerical represenation of the state
        Vertex_state state;
        
        bool operator==(const Vertex& v2) const
        {
            return unit_type == v2.unit_type && num_state == v2.num_state;
        }
    };
public:
    static const unsigned INVALID_VERTEX = 0xffffffff;
    
    typedef M model_type;
    typedef L lattice_type;

    typedef Vertex vertex_type;
    typedef Vertex_state vertex_state_type;
    
    typedef Lat_unit_state lat_unit_state_type;
    typedef typename lattice_type::lat_unit_sites_type lat_unit_sites_type;
    typedef typename lattice_type::lat_unit_bonds_type lat_unit_bonds_type;
    
    typedef boost::multi_array<double, 2> site_tensor_type;
    typedef boost::multi_array<double, 4> bond_tensor_type;
    
    Model(alps::Parameters const& params,
            model_type& model, lattice_type& lattice) :
        model(model),
        lattice(lattice),
        params(params)
    {
        epsilon = params.value_or_default("EPSILON", 0.0);
        if (epsilon <= 0.0 && model.is_signed())
            std::cout << "Warning: Hamiltonian has a sign problem and EPSILON=0; "
                "make sure that EPSILON is ergodic.\n";
        
        _nbstates.resize(lattice.max_site_type() + 1);
        for (unsigned i = 0; i < lattice.nsites(); ++i) {
            unsigned stype = lattice.sitei2alps_type(i);
            _nbstates[stype] = model.site_basis(stype).num_states();
        }
        
        _raising_matrix_elements.resize(lattice.max_site_type() + 1);
        _lowering_matrix_elements.resize(lattice.max_site_type() + 1);
        
        construct_vertices();
    }
    
    std::vector<unsigned> const& nbstates() const
    {
        return _nbstates;
    }
    
    std::vector<std::vector<double> > const& raising_matrix_elements() const
    {
        return _raising_matrix_elements;
    }
    
    std::vector<std::vector<double> > const& lowering_matrix_elements() const
    {
        return _lowering_matrix_elements;
    }
    
    unsigned nvertices() const
    {
        return vertices.size();
    }

    const vertex_type& vertex(unsigned i) const
    {
        return vertices[i];
    }

    double me(unsigned i) const
    {
        return vertices[i].me;
    }

    const vertex_state_type& vertex_state(unsigned i) const
    {
        return vertices[i].state;
    }
    
    bool diagonal(unsigned i) const
    {
        return vertices[i].diagonal;
    }
    
    unsigned diag_vertex_index(lat_unit_state_type& state, unsigned utype,
            lat_unit_sites_type const& sites) const
    {
        return diag_vertex_indices[e_index(state, utype, sites)];
    }
    
    // double max_diag_me() const
    // {
    //     return _max_diag_me;
    // }
    
    double c(unsigned unit_type) const
    {
        return epsilon + _max_diag_me[unit_type];
    }
    
    unsigned find_vertex(
            vertex_state_type const& vstate, unsigned utype) const
    {
        // this function is slow; it's okay
        
        vertex_type vertex;

        vertex.unit_type = utype;
        vertex.num_state = 0;
        for (unsigned k = 2 * UNIT_SIZE - 1; k > 0; --k) {
            lat_unit_sites_type sites;
            lattice.lat_unit_type2sitei(utype, sites);
            
            vertex.num_state = (vertex.num_state + vstate[k])
                * _nbstates[lattice.sitei2alps_type(sites[(k - 1) % UNIT_SIZE])];
        }
        vertex.num_state += vstate[0];
        
        unsigned pos;
        for (pos = 0; pos < vertices.size(); ++pos)
            if (vertex == vertices[pos])
                return pos;

        return INVALID_VERTEX;
    }
private:
    Model();
    
    model_type& model;
    lattice_type& lattice;
    alps::Parameters const& params;

    double epsilon;
    std::vector<double> _max_diag_me;
    
    std::vector<vertex_type> vertices;
    std::vector<unsigned> diag_vertex_indices;
    std::vector<unsigned> diag_vertex_indices_offsets;
    
    std::vector<unsigned> _nbstates;
    
    std::vector<site_tensor_type> site_tensors;
    std::vector<bond_tensor_type> bond_tensors;
    
    std::vector<std::vector<double> > _raising_matrix_elements;
    std::vector<std::vector<double> > _lowering_matrix_elements;
    
    unsigned e_index(state_type const* state, unsigned utype,
        lat_unit_sites_type const& sites) const
    {
        unsigned index = state[0];
        
        for (unsigned k = 1; k < UNIT_SIZE; ++k)
            index = index * _nbstates[lattice.sitei2alps_type(sites[k])] + state[k];
            
        return diag_vertex_indices_offsets[utype] + index;
    }
        
    void construct_vertices()
    {
        unsigned ndstates = 0;
        std::vector<unsigned> nstates(lattice.nlat_unit_types());
        diag_vertex_indices_offsets.resize(lattice.nlat_unit_types());
        
        for (unsigned i = 0; i < lattice.nlat_unit_types(); ++i) {
            diag_vertex_indices_offsets[i] = ndstates;
            
            lat_unit_sites_type sites;
            lattice.lat_unit_type2sitei(i, sites);
            
            unsigned ndstates2 = 1;
            for (unsigned k = 0; k < UNIT_SIZE; ++k)
                ndstates2 *= _nbstates[lattice.sitei2alps_type(sites[k])];
            
            ndstates += ndstates2;
            nstates[i] = ndstates2 * ndstates2;
        }
        
        unsigned invalid_vertex = Model::INVALID_VERTEX;
        diag_vertex_indices.resize(ndstates, invalid_vertex);
        
        unsigned vertex_index = 0;
        std::set<unsigned> site_types;
        std::set<std::string> allops;
        
        bool have_diagonal = false;
        _max_diag_me.resize(lattice.nlat_unit_types());
        
        for (unsigned i = 0; i < lattice.nlat_unit_types(); ++i) {
            lat_unit_sites_type sites;
            lattice.lat_unit_type2sitei(i, sites);
            lat_unit_bonds_type bonds;
            lattice.lat_unit_type2bondi(i, bonds);
            
            // fix me: the following line might not be good; should we have btype = i?
            unsigned btype = lattice.bondi2alps_type(bonds[0]);
            unsigned stype0 = lattice.sitei2alps_type(sites[0]);
            unsigned stype1 = lattice.sitei2alps_type(sites[1]);
            
            lattice.check_homogenity(params, bonds[0]);
            
            alps::Parameters p(params);
            
            lattice.do_inhomogeneous_source(p, bonds[0]);
            site_tensor_type site_tensor0 = alps::get_matrix(double(),
                    model.site_term(stype0), model.site_basis(stype0), p);
            check_site_tensor(site_tensor0, sites[0]);
            
            lattice.do_inhomogeneous_target(p, bonds[0]);
            site_tensor_type site_tensor1 = alps::get_matrix(double(),
                    model.site_term(stype1), model.site_basis(stype1), p);
            check_site_tensor(site_tensor1, sites[1]);
            
            lattice.do_inhomogeneous_bond(p, bonds[0]);
            bond_tensor_type bond_tensor = alps::get_matrix(double(),
                    model.bond_term(btype),
                    model.site_basis(stype0),
                    model.site_basis(stype1), p);
            check_bond_tensor(bond_tensor, sites[0], sites[1]);
            
            std::set<std::string> ops = model.bond_term(btype).operator_names(params);
            allops.insert(ops.begin(), ops.end());
            site_types.insert(stype0);
            site_types.insert(stype1);
            
            unsigned nneighbors0 = lattice.nneighbors(sites[0]);
            unsigned nneighbors1 = lattice.nneighbors(sites[1]);
            
            _max_diag_me[i] = std::numeric_limits<double>::min();
            
            for (unsigned l = 0; l < nstates[i]; ++l) {
                vertex_type vertex;
                
                unsigned state = l;
                for (unsigned k = 0; k < 2 * UNIT_SIZE; ++k) {
                    unsigned nbstates =
                        _nbstates[lattice.sitei2alps_type(sites[k % UNIT_SIZE])];
                    vertex.state[k] = state % nbstates;
                    state /= nbstates;
                }
                
                vertex.unit_type = i;
                vertex.me = bond_matrix_element(vertex, bond_tensor,
                    site_tensor0, site_tensor1, nneighbors0, nneighbors1);
                    
                vertex.diagonal = true;
                for (unsigned k = 0; k < UNIT_SIZE; ++k)
                    if (vertex.state[k] != vertex.state[UNIT_SIZE + k]) {
                        vertex.diagonal = false;
                        break;
                    }
                    
                if (!vertex.diagonal && vertex.me == 0.0)
                    continue;
                    
                vertex.num_state = l;
                
                vertices.push_back(vertex);
                
                if (vertex.diagonal) {
                    if (vertex.me != 0.0)
                        have_diagonal = true;
                    
                    unsigned k = e_index(vertex.state, vertex.unit_type, sites);
                    diag_vertex_indices[k] = vertex_index;
                    
                    if (vertex.me > _max_diag_me[i])
                        _max_diag_me[i] = vertex.me;
                }
                                    
                ++vertex_index;
            }
        }
        
        if (epsilon == 0.0 && !have_diagonal)
            throw std::runtime_error("Hamiltonian looks purely off-diagonal. "
                "Parameter EPSILON has to be non zero for SSE to work.");
                
        if (epsilon <= 0.0)
            epsilon = 1e-6;
        
        std::set<unsigned>::const_iterator sti = site_types.begin();
        for (; sti != site_types.end(); ++sti)
                lowering_and_rasing_operators(allops, *sti);
    }
    
    double bond_matrix_element(vertex_type const& vertex,
        bond_tensor_type const& bond_tensor,
        site_tensor_type const& site_tensor0,
        site_tensor_type const& site_tensor1,
        unsigned nneighbors0, unsigned nneighbors1) const
    {        
        vertex_state_type const& state = vertex.state;
        double me = bond_tensor[state[0]][state[1]][state[2]][state[3]];
        
        if (state[1] == state[3])
            me += site_tensor0[state[0]][state[2]] / nneighbors0;
        if (state[0] == state[2])
            me += site_tensor1[state[1]][state[3]] / nneighbors1;
        
        return me;
    }
    
    void check_site_tensor(site_tensor_type const& site_tensor, unsigned site)
    {
        unsigned dim = _nbstates[lattice.sitei2alps_type(site)];
        
        for (unsigned i = 0; i < dim; ++i)
        for (unsigned j = 0; j < dim; ++j)
            if (i != j && site_tensor[i][j] != 0.0)
                throw std::runtime_error("Offdiagonal site terms are "
                        "not implemented in this SSE code.");
    }
    
    void check_bond_tensor(bond_tensor_type const& bond_tensor,
            unsigned site0, unsigned site1)
    {
        int dim0 = _nbstates[lattice.sitei2alps_type(site0)];
        int dim1 = _nbstates[lattice.sitei2alps_type(site1)];
        
        for (int i = 0; i < dim0; ++i)
        for (int j = 0; j < dim1; ++j)
        for (int k = 0; k < dim0; ++k)
        for (int l = 0; l < dim1; ++l)
            if (bond_tensor[i][j][k][l] != 0.0) {
                if (i + j != k + l)
                    throw std::runtime_error("Offdiagonal bond terms that "
                        "do not conserve quantum numbers are not "
                        "implemented in this SSE code.");
                if (std::abs(i - k) > 1 || std::abs(j - l) > 1)
                    throw std::runtime_error("Offdiagonal bond terms that "
                        "change quantum numbers by 2 or more are not "
                        "implemented in this SSE code.");
            }
    }
    
    void lowering_and_rasing_operators(
        std::set<std::string> const& ops, unsigned stype)
    {
        _raising_matrix_elements[stype].resize(_nbstates[stype]);
        _lowering_matrix_elements[stype].resize(_nbstates[stype]);
        
        unsigned nraising = 0;
        unsigned nlowering = 0;
        
        unsigned nbstates = _nbstates[stype];
        std::set<std::string>::const_iterator it = ops.begin();
        for (; it != ops.end(); ++it) {
            site_tensor_type opmatrix = alps::get_matrix(
                double(), alps::SiteOperator(*it),
                model.basis().site_basis(stype), params);
            
            bool have_raising = false;
            for (unsigned i = 0; i < nbstates - 1; ++i)
                if (opmatrix[i][i + 1] != 0.0) {
                    _raising_matrix_elements[stype][i] = opmatrix[i][i + 1];
                    have_raising = true;
                }
            if (have_raising) ++nraising;
            
            bool have_lowering = false;
            for (unsigned i = 1; i < nbstates; ++i)
                if (opmatrix[i][i - 1] != 0.0) {
                    _lowering_matrix_elements[stype][i] = opmatrix[i][i - 1];
                    have_lowering = true;
                }
            if (have_lowering) ++nlowering;
        }
        
        if (nraising > 1)
            throw std::runtime_error("Found two raising operators.");
        if (nlowering > 1)
            throw std::runtime_error("Found two lowering operators.");
    }
};

#endif
