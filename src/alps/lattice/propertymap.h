/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2004 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

/* $Id$ */

#ifndef ALPS_LATTICE_PROPERTYMAP_H
#define ALPS_LATTICE_PROPERTYMAP_H

#include <alps/config.h>
#include <boost/limits.hpp>
#include <boost/pending/property.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/any.hpp>

#include <string>
#include <vector>

namespace alps {

template <class V, class K=boost::any> // singleton_property_map and constant_property_map
class singleton_property_map {
public:
  typedef K key_type;
  typedef V value_type;
  typedef V& reference;
  typedef boost::lvalue_property_map_tag category;

  singleton_property_map(V v=V()) : v_(v) {}

  operator V () const { return v_;}
  V value() const { return v_;}

  const singleton_property_map<V,K>& operator=(const V& v) { v_=v; return *this;}

  template <class T>
  V& operator[](T ) { return v_;}

  template <class T>
  const V& operator[](T ) const { return v_;}
private:
  V v_;
};

} // end namespace

namespace boost {
template <class V, class K>
V get(const alps::singleton_property_map<V,K>& m, const K&) { return m.value();}

template <class V, class K>
void put(alps::singleton_property_map<V,K>& m, const K& k, const V& v) { m[k]=v;}
}

namespace alps {
namespace detail {

// helper functions to probe for graph properties

template <class T, class DEFAULT=int>
struct existing_property {
  BOOST_STATIC_CONSTANT(bool, result=true);
  typedef T type;
};

template <class DEFAULT>
struct existing_property<boost::detail::error_property_not_found,DEFAULT> {
  BOOST_STATIC_CONSTANT(bool, result=false);
  typedef DEFAULT type;
};

template <bool FLAG, class T1, class T2>
struct choose {
  typedef T1 type;
};

template <class T1, class T2>
struct choose<false,T1,T2> {
  typedef T2 type;
};

} // end namespace detail



template <class Property, class Graph, class Default=int>
struct has_property {
  BOOST_STATIC_CONSTANT(bool, edge_property=false);
  BOOST_STATIC_CONSTANT(bool, vertex_property=false);
  BOOST_STATIC_CONSTANT(bool, graph_property=false);
  BOOST_STATIC_CONSTANT(bool, any_property=false);
  typedef Default vertex_property_type;
  typedef Default edge_property_type;
  typedef Default graph_property_type;
  typedef Default property_type;
  typedef property_type type;
};

template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
struct has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>, D>
{
  typedef boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4> Graph;
  BOOST_STATIC_CONSTANT(bool, edge_property = (detail::existing_property<typename boost::property_value<EP,P>::type,D>::result));
  BOOST_STATIC_CONSTANT(bool, vertex_property = (detail::existing_property<typename boost::property_value<VP,P>::type,D>::result));
  BOOST_STATIC_CONSTANT(bool, graph_property = (detail::existing_property<typename boost::property_value<GP,P>::type,D>::result));
  BOOST_STATIC_CONSTANT(bool, any_property = (edge_property || vertex_property || graph_property));
  typedef typename detail::existing_property<
    typename boost::property_value<EP,P>::type,D>::type edge_property_type;
  typedef typename detail::existing_property<
    typename boost::property_value<VP,P>::type,D>::type vertex_property_type;
  typedef typename detail::existing_property<
    typename boost::property_value<GP,P>::type,D>::type graph_property_type;
  typedef typename detail::choose<edge_property,edge_property_type,
    typename detail::choose<vertex_property,vertex_property_type,
    graph_property_type>::type>::type property_type;
  typedef property_type type;
};

#ifndef BOOST_NO_INCLASS_MEMBER_INITIALIZATION
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::edge_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::vertex_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::graph_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::any_property;
#endif

#if  !defined(__IBMCPP__) || defined (__APPLE__)

template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
struct has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>, D>
{
  typedef boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>  Graph;
  static const bool edge_property = detail::existing_property<
    typename boost::property_value<EP,P>::type,D>::result;
  static const bool vertex_property = detail::existing_property<
    typename boost::property_value<VP,P>::type,D>::result;
  static const bool graph_property = detail::existing_property<
    typename boost::property_value<GP,P>::type,D>::result;
  static const bool any_property =
    edge_property || vertex_property || graph_property;
  typedef typename detail::existing_property<
    typename boost::property_value<EP,P>::type,D>::type edge_property_type;
  typedef typename detail::existing_property<
    typename boost::property_value<VP,P>::type,D>::type vertex_property_type;
  typedef typename detail::existing_property<
    typename boost::property_value<GP,P>::type,D>::type graph_property_type;
  typedef typename detail::choose<edge_property,edge_property_type,
    typename detail::choose<vertex_property,vertex_property_type,
    graph_property_type>::type>::type property_type;
  typedef property_type type;
};

#ifndef BOOST_NO_INCLASS_MEMBER_INITIALIZATION
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::edge_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::vertex_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::graph_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::any_property;
#endif

#endif 


template <class P, class G, class Default>
struct property_map
{
  typedef 
    typename detail::choose<has_property<P,G>::graph_property,
      typename has_property<P,G>::graph_property_type&,
      typename detail::choose<has_property<P,G>::any_property,
        typename boost::property_map<G,P>::type,
        singleton_property_map<Default> 
      >::type
    >::type type;

  typedef
    typename detail::choose<has_property<P,G>::graph_property,
      const typename has_property<P,G>::graph_property_type&,
      typename detail::choose<has_property<P,G>::any_property,
        typename boost::property_map<G,P>::const_type,
        singleton_property_map<Default>
      >::type
    >::type const_type;
};

template <class P, class G, class Default>
struct property_map<P, const G, Default>
{
  typedef
    typename detail::choose<has_property<P,G>::graph_property,
      const typename has_property<P,G>::graph_property_type&,
      typename detail::choose<has_property<P,G>::any_property,
        typename boost::property_map<G,P>::const_type,
        singleton_property_map<Default>
      >::type
    >::type type;

  typedef type const_type;
};


namespace detail {

template <bool F>
struct put_get_helper {};

template <>
struct put_get_helper<true>
{
  template <class P, class G, class T>
  static typename property_map<P,G,int>::type get (P p, G& g, T) {
    return put_get_helper<has_property<P,G>::graph_property>::get_property(p,g);
  }

  template <class P, class G>
  static typename property_map<P,G,int>::type get_property (P p, G& g) 
  { return boost::get_property(g,p);}
};

template <>
struct put_get_helper<false>
{
  template <class P, class G, class V>
  static singleton_property_map<V> get (P, const G&, const V& v)
  { return singleton_property_map<V>(v);}

  template <class P, class G>
  static typename property_map<P,G,int>::type get_property (P p, G& g) 
  { return boost::get(p,g);}
};

} // end namespace detail

template <class P, class G, class V>
inline typename property_map<P,G,V>::type
get_or_default(P p, G& g, const V& v)
{
  return detail::put_get_helper<has_property<P,G>::any_property>::get(p,g,v);
}

template <class SRC, class DST, class PROPERTY, bool has_property>
struct copy_property_helper
{
  template <class SRCREF, class DSTREF>
  static void copy(const SRC&, const SRCREF&, DST&, const DSTREF&) {}
};

template <class SRC, class DST, class PROPERTY>
struct copy_property_helper<SRC,DST,PROPERTY,true>
{
  template <class SRCREF, class DSTREF>
  static void copy(const SRC& s, const SRCREF& sr, DST& d, const DSTREF& dr) {
    boost::put(PROPERTY(),d,dr, boost::get(PROPERTY(),s,sr));
  }
};


/*
template <class Container, class Graph, class Property>
struct container_property_map
{
  typedef typename Container::value_type value_type;
  typedef typename Container::reference reference;
  typedef typename Container::const_reference const_reference;
  typedef typename Container::iterator iterator;
  typedef typename Container::const_iterator const_iterator;
  typedef typename boost::property_map<Graph,Property>::const_type 
                     index_map_type;
  typedef boost::iterator_property_map<iterator, index_map_type,
              value_type, reference> type;
  typedef boost::iterator_property_map<const_iterator, index_map_type,
              value_type, const_reference> const_type;
};

template <class Container, class Graph, class Property>
typename container_property_map<Container,Graph,Property>::type
make_container_property_map(Container& c, const Graph& g, Property p)
{
  return typename container_property_map<Container,Graph,Property>::type
    (c.begin(),boost::get(p,g));
}

template <class Container, class Graph, class Property>
typename container_property_map<Container,Graph,Property>::const_type
make_const_container_property_map(Container& c, const Graph& g, Property p)
{
  return typename container_property_map<Container,Graph,Property>::type
    (c.begin(),boost::get(p,g));
}

*/

} // end namespace alps

#endif // ALPS_LATTICE_GRAPH_PROPERTIES_H
