/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_LATTICE_PROPERTYMAP_H
#define ALPS_LATTICE_PROPERTYMAP_H

#include <alps/config.h>
#include <boost/pending/property.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/any.hpp>
#include <boost/version.hpp>
#include <boost/mpl/if.hpp>


//
// Forward declaration
//
namespace boost {
template <class S1, class S2, class S3, class VP, class EP, class GP, class S4>
class adjacency_list;
}


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

template <class T>
struct get_const_type_member {
    typedef typename T::const_type type;
};

template <bool Found, class PropertyList, class Tag, class Default>
struct found_property_type_or_default_impl {
  BOOST_STATIC_CONSTANT(bool, found=false);
  typedef Default type;
};

template <class PropertyList, class Tag, class Default>
struct found_property_type_or_default_impl<true,PropertyList,Tag,Default> {
  BOOST_STATIC_CONSTANT(bool, found=true);
  typedef typename boost::property_value<PropertyList,Tag>::type type;
};
// helper functions to probe for graph properties


template <class PropertyList, class Tag, class Default=int>
struct found_property_type_or_default
: found_property_type_or_default_impl<
      boost::lookup_one_property<PropertyList,Tag>::found
    , PropertyList
    , Tag
    , Default
> {};



} // end namespace detail



template <class Property, class Graph, class Default=int>
struct has_property {
  BOOST_STATIC_CONSTANT(bool, vertex_property=false);
  BOOST_STATIC_CONSTANT(bool, edge_property=false);
  BOOST_STATIC_CONSTANT(bool, graph_property=false);
  BOOST_STATIC_CONSTANT(bool, any_property=false);
  BOOST_STATIC_CONSTANT(bool, site_property = vertex_property);
  BOOST_STATIC_CONSTANT(bool, bond_property = edge_property);
  typedef Default vertex_property_type;
  typedef Default edge_property_type;
  typedef Default graph_property_type;
  typedef Default property_type;
  typedef property_type type;
  typedef vertex_property_type site_property_type;
  typedef edge_property_type bond_property_type;
};

template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
struct has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>, D>
{
  typedef boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4> Graph;
  BOOST_STATIC_CONSTANT(bool, vertex_property = (detail::found_property_type_or_default<VP,P,D>::found));
  BOOST_STATIC_CONSTANT(bool, edge_property   = (detail::found_property_type_or_default<EP,P,D>::found));
  BOOST_STATIC_CONSTANT(bool, graph_property  = (detail::found_property_type_or_default<GP,P,D>::found));
  BOOST_STATIC_CONSTANT(bool, any_property = (edge_property || vertex_property || graph_property));
  BOOST_STATIC_CONSTANT(bool, site_property = vertex_property);
  BOOST_STATIC_CONSTANT(bool, bond_property = edge_property);
  typedef typename detail::found_property_type_or_default<VP,P,D>::type vertex_property_type;
  typedef typename detail::found_property_type_or_default<EP,P,D>::type edge_property_type;
  typedef typename detail::found_property_type_or_default<GP,P,D>::type graph_property_type;
  typedef typename boost::mpl::if_c<edge_property,edge_property_type,
    typename boost::mpl::if_c<vertex_property,vertex_property_type,
    graph_property_type>::type>::type property_type;
  typedef property_type type;
  typedef vertex_property_type site_property_type;
  typedef edge_property_type bond_property_type;
};

template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::vertex_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::edge_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::graph_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::any_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::site_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::bond_property;

template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
struct has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>, D>
{
  typedef boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4> Graph;
  BOOST_STATIC_CONSTANT(bool, vertex_property = (detail::found_property_type_or_default<VP,P,D>::found));
  BOOST_STATIC_CONSTANT(bool, edge_property   = (detail::found_property_type_or_default<EP,P,D>::found));
  BOOST_STATIC_CONSTANT(bool, graph_property  = (detail::found_property_type_or_default<GP,P,D>::found));
  BOOST_STATIC_CONSTANT(bool, any_property =
    (edge_property || vertex_property || graph_property));
  BOOST_STATIC_CONSTANT(bool, site_property = vertex_property);
  BOOST_STATIC_CONSTANT(bool, bond_property = edge_property);
  typedef typename detail::found_property_type_or_default<VP,P,D>::type vertex_property_type;
  typedef typename detail::found_property_type_or_default<EP,P,D>::type edge_property_type;
  typedef typename detail::found_property_type_or_default<GP,P,D>::type graph_property_type;
  typedef typename boost::mpl::if_c<edge_property,edge_property_type,
    typename boost::mpl::if_c<vertex_property,vertex_property_type,
    graph_property_type>::type>::type property_type;
  typedef property_type type;
  typedef vertex_property_type site_property_type;
  typedef edge_property_type bond_property_type;
};

template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::vertex_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::edge_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::graph_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::any_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::site_property;
template <class s1, class s2, class s3, class VP, class EP, class GP, class s4, class P, class D>
const bool has_property<P, const boost::adjacency_list<s1,s2,s3,VP,EP,GP,s4>,D>::bond_property;

template <class P, class G, class Default>
struct property_map
{
  typedef typename boost::mpl::eval_if_c<
      has_property<P,G>::graph_property
    , boost::add_reference<typename boost::graph_property<G,P>::type>
    , boost::mpl::eval_if_c<
          has_property<P,G>::any_property
        , boost::property_map<G,P>
        , boost::mpl::identity<singleton_property_map<Default> >
    >
  >::type type;

  typedef typename boost::mpl::eval_if_c<
      has_property<P,G>::graph_property
    , boost::add_reference<typename boost::graph_property<G,P>::type const>
    , boost::mpl::eval_if_c<
          has_property<P,G>::any_property
        , detail::get_const_type_member< boost::property_map<G,P> >
        , boost::mpl::identity<singleton_property_map<Default> >
    >
  >::type const_type;
};

template <class P, class G, class Default>
struct property_map<P, const G, Default>
{
  typedef typename boost::mpl::eval_if_c<
      has_property<P,G>::graph_property
    , boost::add_reference<typename boost::graph_property<G,P>::type const>
    , boost::mpl::eval_if_c<
          has_property<P,G>::any_property
        , detail::get_const_type_member< boost::property_map<G,P> >
        , boost::mpl::identity<singleton_property_map<Default> >
    >
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
  static typename property_map<P,G,T>::type get (P p, G& g, T) {
    return put_get_helper<has_property<P,G>::graph_property>::get_property(p,g);
  }

  template <class P, class G>
  static typename property_map<P,G,int>::type get_property (P p, G& g) 
  { 
    BOOST_STATIC_ASSERT((has_property<P,G>::graph_property));
    return boost::get_property(g,p);
  }
};

template <>
struct put_get_helper<false>
{
  template <class P, class G, class V>
  static singleton_property_map<V> get (P, const G&, const V& v)
  { 
    BOOST_STATIC_ASSERT((!has_property<P,G>::any_property));
    return singleton_property_map<V>(v);
  }


  template <class P, class G>
  static typename property_map<P,G,int>::type get_property (P p, G& g) 
  { 
    BOOST_STATIC_ASSERT((has_property<P,G>::any_property));
    BOOST_STATIC_ASSERT((!has_property<P,G>::graph_property));
    return boost::get(p,g);
  }
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

} // end namespace alps

#endif // ALPS_LATTICE_GRAPH_PROPERTIES_H
