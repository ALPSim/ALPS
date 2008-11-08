/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
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

// Weighted Union-Find Algorithm
// Reference:
//   D. Knuth,
//   `The Art of Computer Programming, Vol. 1, Fundamental Algorithms'
//   3rd edition (Addison Wesley, Reading, 1997) Sec 2.3.3.

#ifndef UNION_FIND_H
#define UNION_FIND_H

#include <algorithm> // for std::swap
#include <vector>

namespace union_find {

struct node {
  node() : parent(-1), id() {}
  bool is_root() const { return parent < 0; }
  int weight() const { return -parent; }
  int add_weight(node const& n) { return parent += n.parent; }
  int parent;  // negative for root fragment
  int id;
};

template<class T>
inline int add(std::vector<T>& v)
{
  v.push_back(T());
  return v.size() - 1; // return index of new node
}

template<class T>
inline int root_index(std::vector<T> const& v, int g)
{ while (!v[g].is_root()) g = v[g].parent; return g; }

template<class T>
inline T const& root(std::vector<T> const& v, int g)
{ return v[root_index(v, g)]; }

template<class T>
inline T const& root(std::vector<T> const& v, T const& n)
{ return n.is_root() ? n : root(v, n.parent); }

template<class T>
inline int cluster_id(std::vector<T> const& v, int g)
{ return root(v, g).id; }

template<class T>
inline int cluster_id(std::vector<T> const& v, T const& n)
{ return root(v, n).id; }

template<class T>
inline void update_link(std::vector<T>& v, int g, int r)
{ while (g != r) { int p = v[g].parent; v[g].parent = r; g = p; } }

template<class T>
inline int unify(std::vector<T>& v, int g0, int g1)
{
  using std::swap;
  int r0 = root_index(v, g0);
  int r1 = root_index(v, g1);
  if (r0 != r1) {
    if (v[r0].weight() < v[r1].weight()) swap(r0, r1);
    v[r0].add_weight(v[r1]);
    v[r1].parent = r0;
  }
  update_link(v, g0, r0);
  update_link(v, g1, r0);
  return r0; // return (new) root node
}

} // end namespace union_find

#endif
