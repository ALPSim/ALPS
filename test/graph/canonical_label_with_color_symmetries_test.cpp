/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>                          *
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

#include <alps/graph/canonical_properties.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <iostream>

typedef boost::property<alps::edge_type_t,unsigned int> edge_props;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_props> graph_type;
typedef boost::graph_traits<graph_type>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<graph_type>::edge_iterator edge_iterator;
typedef boost::property_map<graph_type,alps::edge_type_t>::type edge_color_map_type;


template <typename Map>
void print_colormap(Map const& m)
{
    for(typename Map::const_iterator it(m.begin()), end(m.end()); it != end; ++it)
        std::cout << it->first << "->" << it->second << " ";
    std::cout << std::endl;
}

bool colored_edges_with_color_symmetry_test1() {
    std::cout << "colored_edges_with_color_symmetry_test1()" << std::endl;
    typedef alps::graph::graph_label<graph_type>::type label_type;
    using alps::graph::canonical_properties;

    // g: O---O+++O+++O
    graph_type g(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(1, 2, g).first;
        edge_color[e] = 1;
        e = add_edge(2, 3, g).first;
        edge_color[e] = 1;
    }

    // h: O+++O---O---O
    graph_type h(g);
    {
        unsigned int const map[] = {1,0};
        edge_iterator it,end;
        for(boost::tie(it,end) = edges(h); it != end; ++it)
            get(alps::edge_type_t(),h)[*it] = map[get(alps::edge_type_t(),h)[*it]];
    }

    // i: O---O+++O---O
    graph_type i(g);
    {
        edge_descriptor e = edge(2,3,i).first;
        get(alps::edge_type_t(),i)[e] = 0;
    }

    // j: O+++O---O---O
    graph_type j(i);
    {
        unsigned int const map[] = {1,0};
        edge_iterator it,end;
        for(boost::tie(it,end) = edges(j); it != end; ++it)
            get(alps::edge_type_t(),j)[*it] = map[get(alps::edge_type_t(),j)[*it]];
    }


    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 0;

    label_type lg(get<1>(canonical_properties(g)));
    label_type lh(get<1>(canonical_properties(h)));
    label_type li(get<1>(canonical_properties(i)));
    label_type lj(get<1>(canonical_properties(j)));


    label_type lg_with_sym(get<1>(canonical_properties(g,color_symmetry)));
    label_type lh_with_sym(get<1>(canonical_properties(h,color_symmetry)));
    label_type li_with_sym(get<1>(canonical_properties(i,color_symmetry)));
    label_type lj_with_sym(get<1>(canonical_properties(j,color_symmetry)));

    std::cout << lg << std::endl;
    std::cout << lh << std::endl;
    std::cout << li << std::endl;
    std::cout << lj << std::endl;
    std::cout << lg_with_sym << std::endl;
    std::cout << lh_with_sym << std::endl;
    std::cout << li_with_sym << std::endl;
    std::cout << lj_with_sym << std::endl;

    // False statements
    std::cout << std::boolalpha
        << (lg == lh) << (lg == li) << (lg == lj)
        << (lh == li) << (lh == lj)
        << (li == lj) << std::endl;

    // True statements
    std::cout << std::boolalpha
        << (lg_with_sym == lh_with_sym)
        << (li_with_sym == lj_with_sym) << std::endl;

    // False statements
    std::cout << std::boolalpha
        << (lg_with_sym == li_with_sym) << std::endl;

    return true;
}

bool colored_edges_with_color_symmetry_test2() {
    std::cout << "colored_edges_with_color_symmetry_test2()" << std::endl;
    typedef alps::graph::graph_label<graph_type>::type label_type;
    using alps::graph::canonical_properties;

    // g:
    //       4              c0 ---
    //       +              c1 +++
    //   1---0+++2          c2 ===
    //   \\ / 
    //     3
    graph_type g(5);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0,1,g).first;
        edge_color[e] = 0;
        e = add_edge(0,2,g).first;
        edge_color[e] = 1;
        e = add_edge(0,3,g).first;
        edge_color[e] = 0;
        e = add_edge(0,4,g).first;
        edge_color[e] = 1;
        e = add_edge(1,3,g).first;
        edge_color[e] = 2;
    }

    // h:
    //       4              c0 ---
    //       +              c1 +++
    //   1===0+++2          c2 ===
    //    \ // 
    //     3
    graph_type h(g);
    {
        unsigned int const map[] = {2,1,0};
        edge_iterator it,end;
        for(boost::tie(it,end) = edges(h); it != end; ++it)
            get(alps::edge_type_t(),h)[*it] = map[get(alps::edge_type_t(),h)[*it]];
    }


    // i:
    //       4              c0 ---
    //       ||             c1 +++
    //   1---0===2          c2 ===
    //    + / 
    //     3
    graph_type i(g);
    {
        unsigned int const map[] = {0,2,1};
        edge_iterator it,end;
        for(boost::tie(it,end) = edges(i); it != end; ++it)
            get(alps::edge_type_t(),i)[*it] = map[get(alps::edge_type_t(),i)[*it]];
    }


    // color   grp0  grp1
    // colors (0,2)  (1)
    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 1;
    color_symmetry[2] = 0;

    label_type lg(get<1>(canonical_properties(g)));
    label_type lh(get<1>(canonical_properties(h)));
    label_type li(get<1>(canonical_properties(i)));

    label_type lg_with_sym(get<1>(canonical_properties(g,color_symmetry)));
    label_type lh_with_sym(get<1>(canonical_properties(h,color_symmetry)));
    label_type li_with_sym(get<1>(canonical_properties(i,color_symmetry)));

    std::cout << lg << std::endl;
    std::cout << lh << std::endl;
    std::cout << li << std::endl;
    std::cout << lg_with_sym << std::endl;
    std::cout << lh_with_sym << std::endl;
    std::cout << li_with_sym << std::endl;

    // False statements
    std::cout << std::boolalpha
        << (lg == lh) << (lg == li)
        << (lh == li) << std::endl;

    // True statements
    std::cout << std::boolalpha
        << (lg_with_sym == lh_with_sym) << std::endl;

    // False statements
    std::cout << std::boolalpha
        << (lg_with_sym == li_with_sym) << std::endl;

    return true;
}

bool colored_edges_with_color_symmetry_test3() {
    std::cout << "colored_edges_with_color_symmetry_test3()" << std::endl;
    typedef alps::graph::graph_label<graph_type>::type label_type;
    using alps::graph::canonical_properties;

    // g: O---O---O---O
    graph_type g(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(1, 2, g).first;
        edge_color[e] = 0;
        e = add_edge(2, 3, g).first;
        edge_color[e] = 0;
    }

    // h: O+++O+++O+++O
    graph_type h(g);
    {
        unsigned int const map[] = {1,0};
        edge_iterator it,end;
        for(boost::tie(it,end) = edges(h); it != end; ++it)
            get(alps::edge_type_t(),h)[*it] = map[get(alps::edge_type_t(),h)[*it]];
    }

    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 0;

    label_type lg(get<1>(canonical_properties(g)));
    label_type lh(get<1>(canonical_properties(h)));

    label_type lg_with_sym(get<1>(canonical_properties(g,color_symmetry)));
    label_type lh_with_sym(get<1>(canonical_properties(h,color_symmetry)));

    std::cout << lg << std::endl;
    std::cout << lh << std::endl;
    std::cout << lg_with_sym << std::endl;
    std::cout << lh_with_sym << std::endl;

    // False statements
    std::cout << std::boolalpha
        << (lg == lh) << std::endl;

    // True statements
    std::cout << std::boolalpha
        << (lg_with_sym == lh_with_sym) << std::endl;
    return true;
}

void colored_edges_with_color_symmetry_test4() {
    std::cout << "colored_edges_with_color_symmetry_test4()" << std::endl;
    typedef alps::graph::graph_label<graph_type>::type label_type;
    using alps::graph::canonical_properties;

    // g:
    //   2++++4     c0 +++
    //  /    /      c1 ---
    // 0++++1       c2 ...
    //  . 
    //   3
    graph_type g(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 2, g).first;
        edge_color[e] = 1;
        e = add_edge(0, 3, g).first;
        edge_color[e] = 2;
        e = add_edge(1, 4, g).first;
        edge_color[e] = 1;
        e = add_edge(2, 4, g).first;
        edge_color[e] = 0;
    }

    // h:
    //        2++++4     c0 +++
    //       /    /      c1 ---
    //      0++++1       c2 ...
    //     /
    //    3
    graph_type h(g);
    {
        unsigned int const map[] = {0,1,0};
        edge_iterator it,end;
        for(boost::tie(it,end) = edges(h); it != end; ++it)
            get(alps::edge_type_t(),h)[*it] = map[get(alps::edge_type_t(),h)[*it]];
    }

    // i:
    //        2++++4     c0 +++
    //       /    /      c1 ---
    // 3++++0++++1       c2 ...
    //
    graph_type i(g);
    {
        unsigned int const map[] = {0,1,1};
        edge_iterator it,end;
        for(boost::tie(it,end) = edges(i); it != end; ++it)
            get(alps::edge_type_t(),i)[*it] = map[get(alps::edge_type_t(),i)[*it]];
    }
    // color   grp0
    // colors (0,1,2)
    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 0;
    color_symmetry[2] = 0;

    label_type lg(get<1>(canonical_properties(g)));
    label_type lh(get<1>(canonical_properties(h)));
    label_type li(get<1>(canonical_properties(i)));

    label_type lg_with_sym(get<1>(canonical_properties(g,color_symmetry)));
    label_type lh_with_sym(get<1>(canonical_properties(h,color_symmetry)));
    label_type li_with_sym(get<1>(canonical_properties(i,color_symmetry)));

    std::cout << lg << std::endl;
    std::cout << lh << std::endl;
    std::cout << li << std::endl;
    std::cout << lg_with_sym << std::endl;
    std::cout << lh_with_sym << std::endl;
    std::cout << li_with_sym << std::endl;

    // False statements
    std::cout << std::boolalpha
        << (lg == lh) << (lg == li)
        << (lh == li) << std::endl;

    // True statements
    std::cout << std::boolalpha
        << (lh_with_sym == li_with_sym) << std::endl;

    // False statements
    std::cout << std::boolalpha
        << (lg_with_sym == lh_with_sym) << std::endl;
}

void colored_edges_with_color_symmetry_test5() {
    std::cout << "colored_edges_with_color_symmetry_test5()" << std::endl;
    typedef alps::graph::graph_label<graph_type>::type label_type;
    typedef boost::container::flat_map<alps::type_type, alps::type_type> color_map_type;
    using alps::graph::canonical_properties;

    // g:
    //     3        c0 +++
    //     |        c1 ---
    // 2+++1+++0
    //
    graph_type g(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(1, 2, g).first;
        edge_color[e] = 0;
        e = add_edge(1, 3, g).first;
        edge_color[e] = 1;
    }

    // h:
    //     1        c0 +++
    //     +        c1 ---
    // 2---0---3
    //
    graph_type h(4);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),h);
        edge_descriptor e;
        e = add_edge(0, 1, h).first;
        edge_color[e] = 0;
        e = add_edge(0, 2, h).first;
        edge_color[e] = 1;
        e = add_edge(0, 3, h).first;
        edge_color[e] = 1;
    }

    // color   grp0
    // colors (0,1)
    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0;
    color_symmetry[1] = 0;

    label_type lg(get<1>(canonical_properties(g)));
    label_type lh(get<1>(canonical_properties(h)));

    color_map_type lg_colormap;
    color_map_type lh_colormap;
    label_type lg_with_sym(get<1>(canonical_properties(g, color_symmetry, &lg_colormap)));
    label_type lh_with_sym(get<1>(canonical_properties(h, color_symmetry, &lh_colormap)));

    std::cout << lg << std::endl;
    std::cout << lh << std::endl;
    std::cout << lg_with_sym << std::endl;
    print_colormap(lg_colormap);
    std::cout << lh_with_sym << std::endl;
    print_colormap(lh_colormap);

    // False statements
    std::cout << std::boolalpha
        << (lg == lh) << std::endl;

    // True statements
    std::cout << std::boolalpha
        << (lg_with_sym == lh_with_sym) << std::endl;
}

void colored_edges_with_color_symmetry_test6() {
    std::cout << "colored_edges_with_color_symmetry_test6()" << std::endl;
    typedef alps::graph::graph_label<graph_type>::type label_type;
    typedef boost::container::flat_map<alps::type_type, alps::type_type> color_map_type;
    using alps::graph::canonical_properties;

    // g:
    //              c0 +++
    //              c2 ...
    // 1+++0...2
    //
    graph_type g(3);
    {
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        edge_descriptor e;
        e = add_edge(0, 1, g).first;
        edge_color[e] = 0;
        e = add_edge(0, 2, g).first;
        edge_color[e] = 2;
    }

    // color   grp0
    // colors (0,1,2)
    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0; 
    color_symmetry[1] = 0;
    color_symmetry[2] = 0;

    label_type lg(get<1>(canonical_properties(g)));

    color_map_type lg_colormap;
    label_type lg_with_sym(get<1>(canonical_properties(g, color_symmetry, &lg_colormap)));

    std::cout << lg << std::endl;
    std::cout << lg_with_sym << std::endl;
    print_colormap(lg_colormap);
}

void colored_edges_with_color_symmetry_test7() {
    std::cout << "colored_edges_with_color_symmetry_test7()" << std::endl;
    typedef alps::graph::graph_label<graph_type>::type label_type;
    typedef boost::container::flat_map<alps::type_type, alps::type_type> color_map_type;
    using alps::graph::canonical_properties;
    graph_type g(8);
    {
        edge_descriptor e;
        edge_color_map_type edge_color = get(alps::edge_type_t(),g);
        e = add_edge(0, 1, g).first; edge_color[e] = 1;
        e = add_edge(2, 1, g).first; edge_color[e] = 1;
        e = add_edge(3, 1, g).first; edge_color[e] = 0;
        e = add_edge(4, 0, g).first; edge_color[e] = 0;
        e = add_edge(2, 0, g).first; edge_color[e] = 0;
        e = add_edge(5, 3, g).first; edge_color[e] = 0;
        e = add_edge(2, 3, g).first; edge_color[e] = 1;
        e = add_edge(6, 0, g).first; edge_color[e] = 1;
        e = add_edge(7, 3, g).first; edge_color[e] = 1;
    }
    graph_type h(8);
    {
        edge_descriptor e;
        edge_color_map_type edge_color = get(alps::edge_type_t(), h);
        e = add_edge(4, 7, h).first; edge_color[e] = 1;
        e = add_edge(4, 5, h).first; edge_color[e] = 1;
        e = add_edge(4, 6, h).first; edge_color[e] = 0;
        e = add_edge(0, 7, h).first; edge_color[e] = 0;
        e = add_edge(5, 7, h).first; edge_color[e] = 0;
        e = add_edge(3, 6, h).first; edge_color[e] = 0;
        e = add_edge(5, 6, h).first; edge_color[e] = 1;
        e = add_edge(1, 7, h).first; edge_color[e] = 1;
        e = add_edge(2, 6, h).first; edge_color[e] = 1;
    }

    alps::graph::color_partition<graph_type>::type color_symmetry;
    color_symmetry[0] = 0; 
    color_symmetry[1] = 0;

    color_map_type lg_colormap;
    label_type lg(get<1>(canonical_properties(g)));
    label_type lg_with_sym(get<1>(canonical_properties(g, color_symmetry, &lg_colormap)));

    color_map_type lh_colormap;
    label_type lh(get<1>(canonical_properties(h)));
    label_type lh_with_sym(get<1>(canonical_properties(h, color_symmetry, &lh_colormap)));
    // True statements
    std::cout << std::boolalpha
        << (lg == lh)
        << (lg_with_sym == lh_with_sym) << std::endl;
    std::cout << lg << std::endl;
    std::cout << lg_with_sym << std::endl;
    print_colormap(lg_colormap);
    std::cout << lh << std::endl;
    std::cout << lh_with_sym << std::endl;
    print_colormap(lh_colormap);
}

int main() {
    colored_edges_with_color_symmetry_test1();
    colored_edges_with_color_symmetry_test2();
    colored_edges_with_color_symmetry_test3();
    colored_edges_with_color_symmetry_test4();
    colored_edges_with_color_symmetry_test5();
    colored_edges_with_color_symmetry_test6();
    colored_edges_with_color_symmetry_test7();
    return 0;
}
