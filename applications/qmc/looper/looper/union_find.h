/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
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

// Weighted Union-Find Algorithm
// Reference:
//   D. Knuth,
//   `The Art of Computer Programming, Vol. 1, Fundamental Algorithms'
//   3rd edition (Addison Wesley, Reading, 1997) Sec 2.3.3.

#ifndef LOOPER_UNION_FIND_H
#define LOOPER_UNION_FIND_H

#include <algorithm> // for std::swap
#include <vector>
#include <iostream>

namespace looper {
namespace union_find {

class node {
public:
  node() : parent_(-1) {} // root node with weight = 1
  bool is_root() const { return parent_ < 0; }
  void set_parent(int parent) { parent_ = parent; }
  int parent() const { return parent_; }
  int add_weight(node const& n) { return parent_ += n.parent_; }
  int weight() const { return -parent_; }
  void set_id(int id) { id_ = id; }
  int id() const { return id_; }
private:
  int parent_; // negative for root fragment
  int id_;
};

class node_noweight {
public:
  node_noweight() : parent_(id_mask) {} // root note with id = 0
  bool is_root() const { return parent_ < 0; }
  void set_parent(int parent) { parent_ = parent; }
  int parent() const { return parent_; }
  void set_id(int id) { parent_ = id ^ id_mask; }
  int id() const { return parent_ ^ id_mask; }
private:
  static const int id_mask = -1;
  int parent_; // negative for root fragment
};

template<class T>
inline int add(std::vector<T>& v) {
  v.push_back(T());
  return v.size() - 1; // return index of new node
}

template<class T>
inline int root_index(std::vector<T> const& v, int g) {
  while (!v[g].is_root()) g = v[g].parent(); return g;
}

template<class T>
inline T const& root(std::vector<T> const& v, int g) { return v[root_index(v, g)]; }

template<class T>
inline T const& root(std::vector<T> const& v, T const& n) {
  return n.is_root() ? n : root(v, n.parent());
}

template<class T>
inline int cluster_id(std::vector<T> const& v, int g) { return root(v, g).id(); }

template<class T>
inline int cluster_id(std::vector<T> const& v, T const& n) { return root(v, n).id(); }

template<class T>
inline void set_root(std::vector<T>& v, int g) {
  int r = root_index(v, g);
  if (r != g) {
    v[g] = v[r];
    v[r].set_parent(g);
  }
}

template<class T>
inline void update_link(std::vector<T>& v, int g, int r) {
  while (g != r) {
    int p = v[g].parent();
    v[g].set_parent(r);
    g = p;
  }
}

template<class T>
inline int unify_weight(std::vector<T>& v, int g0, int g1) {
  using std::swap;
  int r0 = root_index(v, g0);
  int r1 = root_index(v, g1);
  if (r0 != r1) {
    if (v[r0].weight() < v[r1].weight()) swap(r0, r1);
    v[r0].add_weight(v[r1]);
    v[r1].set_parent(r0);
  }
  update_link(v, g0, r0);
  update_link(v, g1, r0);
  return r0; // return (new) root node
}

template<class T>
inline int unify_noweight(std::vector<T>& v, int g0, int g1) {
  int r0 = root_index(v, g0);
  int r1 = root_index(v, g1);
  if (r0 != r1) v[r1].set_parent(r0);
  update_link(v, g0, r0);
  update_link(v, g1, r0);
  return r0; // return (new) root node
}

inline int unify(std::vector<node>& v, int g0, int g1) {
  return unify_weight(v, g0, g1);
}

inline int unify(std::vector<node_noweight>& v, int g0, int g1) {
  return unify_noweight(v, g0, g1);
}

template<class T>
inline void output(std::vector<T> const& v, std::ostream& os = std::cout) {
  for (int i = 0; i < v.size(); ++i) {
    os << "node " << i << ": ";
    int g = i;
    if (v[g].is_root()) {
      os << "root (id = " << v[g].id() << ")" << std::endl;
    } else {
      while (true) {
        if (v[g].is_root()) {
          os << g << std::endl;
          break;
        } else {
          os << g << " -> ";
          g = v[g].parent();
        }
      }
    }
  }
}

} // end namespace union_find
} // end namespace looper

#endif
