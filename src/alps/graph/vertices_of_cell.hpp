#ifndef ALPS_GRAPH_VERTICES_OF_CELL_HPP
#define ALPS_GRAPH_VERTICES_OF_CELL_HPP

#include <vector>
#include <boost/graph/graph_traits.hpp>
#include <alps/lattice/lattice.h>

namespace alps {
namespace graph {

    template <typename Lattice>
    std::vector<typename boost::graph_traits<typename alps::graph_traits<Lattice>::graph_type>::vertex_descriptor>
    vertices_of_cell(typename alps::lattice_traits<Lattice>::cell_descriptor c, Lattice const& L)
    {
        typedef typename alps::graph_traits<Lattice>::graph_type graph_type;
        typename alps::lattice_traits<Lattice>::size_type const cell_id = index(c, L);
        std::size_t unit_cell_size = num_vertices(alps::graph::graph(unit_cell(L)));
        std::vector<typename boost::graph_traits<graph_type>::vertex_descriptor> vertices;
        for(unsigned int v = 0; v < unit_cell_size; ++v)
            vertices.push_back(cell_id * unit_cell_size + v);
        return vertices;
    }

} // end namespace graph
} // end namespace alps

#endif // ALPS_GRAPH_VERTICES_OF_CELL_HPP
