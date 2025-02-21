/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef ALPS_GRAPH_CANONICAL_PROPERTIES_TRAITS_HPP
#define ALPS_GRAPH_CANONICAL_PROPERTIES_TRAITS_HPP

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/container/container_fwd.hpp>
#include <alps/lattice/graphproperties.h>
#include <alps/lattice/propertymap.h>

namespace alps {
namespace graph {
namespace detail {

    class graph_label_matrix_type : public boost::dynamic_bitset<> {
        public:
            bool operator< (graph_label_matrix_type const & rhs) const {
                return size() < rhs.size() || (
                       !(rhs.size() < size())
                    && static_cast<boost::dynamic_bitset<> const &>(*this) < static_cast<boost::dynamic_bitset<> const &> (rhs)
                );
            }
            bool operator == (graph_label_matrix_type const & rhs) const {
                return size() == rhs.size()
                    && static_cast<boost::dynamic_bitset<> const &>(*this) == static_cast<boost::dynamic_bitset<> const &>(rhs);
            }
    };

    // printable color list
    template<typename Value> class graph_label_color_vector : public std::vector<Value>
    {
      public:
        graph_label_color_vector() {}
        template<typename A> graph_label_color_vector(std::vector<Value,A> const& v) : std::vector<Value>(v) {}
    };

    template <typename T>
    void swap(graph_label_color_vector<T>& v1, graph_label_color_vector<T>& v2) {
        swap(static_cast<std::vector<T>&>(v1), static_cast<std::vector<T>&>(v2));
    }

    // ostream operator for color list
    template<typename Stream, typename Value> Stream & operator<< (Stream & os, graph_label_color_vector<Value> const & vec) {
        os << "(";
        for (typename graph_label_color_vector<Value>::const_iterator it = vec.begin(); it != vec.end(); ++it)
            os << (it == vec.begin() ? "" : " ") << *it;
        os << ")";
        return os;
    }

    // istream operator for color list
    template<typename Stream, typename Value> Stream & operator>> (Stream & is, graph_label_color_vector<Value> & vec) {
        typename Stream::char_type c = 0;
        is >> c;
        if (!is || c != '(')
            throw std::invalid_argument("could not read color list (start not found)");
        for (Value v; is >> v; )
            vec.push_back(v);
        if (!is.bad() && !is.eof())
        {
            is.clear();
            if (is >> c && c == ')')
                return is;
        }
        throw std::invalid_argument("could not read color list (end not found)");
    }

    // no coloring
    template<typename Graph, bool, bool> struct graph_label_helper {
        typedef boost::tuple<
              // #vertices * (#vertices + 1) / 2 bits: triangular adjacency matrix
              graph_label_matrix_type
        > type;
    };

    // vertex coloring
    template<typename Graph> struct graph_label_helper<Graph, true, false> {
        typedef boost::tuple<
              // #vertices * (#vertices + 1) / 2 bits: triangular adjacency matrix
              graph_label_matrix_type
              // #vertices * (#vertex colors) bits: vertex vs color matrix
            , graph_label_matrix_type
              // vertex color list
            , graph_label_color_vector<typename has_property<alps::vertex_type_t,Graph>::vertex_property_type>
        > type;
    };

    // edge coloring
    template<typename Graph> struct graph_label_helper<Graph, false, true> {
        typedef boost::tuple<
              // #vertices * (#vertices + 1) / 2 bits: triangular adjacency matrix
              graph_label_matrix_type
              // #edge * (#edge colors) bits: edge vs color matrix
            , graph_label_matrix_type
              // edge color list
            , graph_label_color_vector<typename has_property<alps::edge_type_t,Graph>::edge_property_type>
        > type;
    };

    // vertex and edge coloring
    template<typename Graph> struct graph_label_helper<Graph, true, true> {
        typedef boost::tuple<
              // #vertices * (#vertices + 1) / 2 bits: triangular adjacency matrix
              graph_label_matrix_type
              // #vertices * (#vertex colors) bits: vertex vs color matrix
            , graph_label_matrix_type
              // vertex color list
            , graph_label_color_vector<typename has_property<alps::vertex_type_t,Graph>::vertex_property_type>
              // #edge * (#edge colors) bits: edge vs color matrix
            , graph_label_matrix_type
              // edge color list
            , graph_label_color_vector<typename has_property<alps::edge_type_t,Graph>::edge_property_type>
        > type;
    };

    template<typename Graph, bool> struct color_partition_helper {
    };

    template<typename Graph> struct color_partition_helper<Graph, true> {
        typedef boost::container::flat_map<typename boost::property_map<Graph,alps::edge_type_t>::type::value_type, unsigned int> type;
    };

} // end namespace detail

// pi = (V1, V2, ..., Vr), Vi = (n1, n2, ..., nk), ni element of G
template<typename Graph> struct partition_type {
    typedef std::vector<std::vector<typename boost::graph_traits<Graph>::vertex_descriptor> > type;
};

// comparable graph label
// vertex coloring tag: alps::vertex_type_t
// edge coloring tag: alps::edge_type_t
template<typename Graph> struct graph_label {
    typedef typename detail::graph_label_helper<
          Graph
        , has_property<alps::vertex_type_t,Graph>::vertex_property
        , has_property<alps::edge_type_t,Graph>::edge_property
    >::type type;
};

template<typename Graph> struct color_partition : detail::color_partition_helper<Graph, has_property<alps::edge_type_t,Graph>::edge_property> {
};

struct no_color_symmetry {};

template<typename Graph> struct canonical_properties_type {
    typedef boost::tuple<
      // canonical ordering
      std::vector<typename boost::graph_traits<Graph>::vertex_descriptor>
      // canonical label
    , typename graph_label<Graph>::type
      // orbit partition
    , typename partition_type<Graph>::type
    > type;
};

// Enum to select the canonical properties with boost::get<enum>
enum { ordering, label, partition };

} // end namespace graph
} // end namespace alps

#endif // ALPS_GRAPH_CANONICAL_PROPERTIES_TRAITS_HPP
