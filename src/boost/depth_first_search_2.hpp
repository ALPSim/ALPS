#ifndef BOOST_GRAPH_DFS2_HPP
#define BOOST_GRAPH_DFS2_HPP

// Another implementation of the depth first search algorithm without
// using recursive function calls

#include <boost/graph/depth_first_search.hpp>
#include <stack>

namespace boost {

  namespace detail {

    template <class IncidenceGraph, class DFSVisitor, class ColorMap>
    void depth_first_visit_2_impl
      (const IncidenceGraph& g,
       typename graph_traits<IncidenceGraph>::vertex_descriptor u, 
       DFSVisitor vis,
       ColorMap color)
    {
      function_requires<IncidenceGraphConcept<IncidenceGraph> >();
      function_requires<DFSVisitorConcept<DFSVisitor, IncidenceGraph> >();
      typedef typename graph_traits<IncidenceGraph>::vertex_descriptor Vertex;
      function_requires< ReadWritePropertyMapConcept<ColorMap, Vertex> >();
      typedef typename property_traits<ColorMap>::value_type ColorValue;
      function_requires< ColorValueConcept<ColorValue> >();
      typedef color_traits<ColorValue> Color;
      typename graph_traits<IncidenceGraph>::out_edge_iterator ei, ei_end;

      std::stack<
        std::pair<Vertex,
	typename graph_traits<IncidenceGraph>::out_edge_iterator> > vs;

      vs.push(std::make_pair(u, out_edges(u, g).first));
      put(color, u, Color::gray());        vis.discover_vertex(u, g);
      while (!vs.empty()) {
	tie(u, ei) = vs.top();
	bool discover = false;
	for (ei_end = out_edges(u, g).second; ei != ei_end; ++ei) {
	  Vertex v = boost::target(*ei, g);       vis.examine_edge(*ei, g);
	  ColorValue v_color = get(color, v);
	  if (v_color == Color::white()) { vis.tree_edge(*ei, g);
	    vs.top().second = ++ei; // save NEXT out edge
	    vs.push(std::make_pair(v, out_edges(v, g).first));
	    put(color, v, Color::gray());  vis.discover_vertex(v, g);
	    discover = true;
	    break;
	  } else if (v_color == Color::gray()) {
	                                   vis.back_edge(*ei, g);
	  } else                           vis.forward_or_cross_edge(*ei, g);
	}
	if (!discover) {
	  put(color, u, Color::black());   vis.finish_vertex(u, g);
	  vs.pop();
	}
      }
    }
  } // namespace detail

  template <class VertexListGraph, class DFSVisitor, class ColorMap, 
            class Vertex>
  void
  depth_first_search_2(const VertexListGraph& g, DFSVisitor vis,
		       ColorMap color, Vertex start_vertex)
  {
    function_requires<DFSVisitorConcept<DFSVisitor, VertexListGraph> >();
    typedef typename property_traits<ColorMap>::value_type ColorValue;
    typedef color_traits<ColorValue> Color;

    typename graph_traits<VertexListGraph>::vertex_iterator ui, ui_end;
    for (tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      put(color, *ui, Color::white());       vis.initialize_vertex(*ui, g);
    }

    if (start_vertex != *vertices(g).first){ vis.start_vertex(start_vertex, g);
      detail::depth_first_visit_2_impl(g, start_vertex, vis, color);
    }

    for (tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui) {
      ColorValue u_color = get(color, *ui);
      if (u_color == Color::white()) {       vis.start_vertex(*ui, g);
        detail::depth_first_visit_2_impl(g, *ui, vis, color);
      }
    }
  }

  template <class VertexListGraph, class DFSVisitor, class ColorMap>
  void
  depth_first_search_2(const VertexListGraph& g, DFSVisitor vis,
		       ColorMap color)
  {
    depth_first_search_2(g, vis, color, *vertices(g).first);
  }

  namespace detail {
    template <class ColorMap>
    struct dfs2_dispatch {

      template <class VertexListGraph, class Vertex, class DFSVisitor, 
                class P, class T, class R>
      static void
      apply(const VertexListGraph& g, DFSVisitor vis, Vertex start_vertex,
            const bgl_named_params<P, T, R>&,
            ColorMap color)
      {
        depth_first_search_2(g, vis, color, start_vertex);
      }
    };

    template <>
    struct dfs2_dispatch<detail::error_property_not_found> {
      template <class VertexListGraph, class Vertex, class DFSVisitor,
                class P, class T, class R>
      static void
      apply(const VertexListGraph& g, DFSVisitor vis, Vertex start_vertex,
            const bgl_named_params<P, T, R>& params,
            detail::error_property_not_found)
      {
        std::vector<default_color_type> color_vec(num_vertices(g));
        default_color_type c = white_color; // avoid warning about un-init
        depth_first_search_2
          (g, vis, make_iterator_property_map
           (color_vec.begin(),
            choose_const_pmap(get_param(params, vertex_index),
                              g, vertex_index), c), 
           start_vertex);
      }
    };

  } // namespace detail

  template <class VertexListGraph, class P, class T, class R>
  void
  depth_first_search_2(const VertexListGraph& g, 
		       const bgl_named_params<P, T, R>& params)
  {
    typedef typename property_value< bgl_named_params<P, T, R>, 
      vertex_color_t>::type C;
    detail::dfs2_dispatch<C>::apply
      (g,
       choose_param(get_param(params, graph_visitor),
                    make_dfs_visitor(null_visitor())),
       choose_param(get_param(params, root_vertex_t()),
                    *vertices(g).first),
       params,
       get_param(params, vertex_color)
       );
  }
  
  template <class IncidenceGraph, class DFSVisitor, class ColorMap>
  void depth_first_visit_2
    (const IncidenceGraph& g,
     typename graph_traits<IncidenceGraph>::vertex_descriptor u, 
     DFSVisitor vis, ColorMap color)
  {
    detail::depth_first_visit_2_impl(g, u, vis, color);
  }

} // namespace boost

#endif // BOOST_GRAPH_DFS2_HPP
