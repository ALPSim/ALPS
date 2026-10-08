/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2004 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <cmath>
#include <boost/random.hpp>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include <unordered_set>

// Track lifetime outside object storage: reading init_ before construction in
// the old test double was undefined and made sanitizer results unreliable.
template<class T>
struct non_pod {
  inline static std::unordered_set<const non_pod*> live;
  T data_;
  non_pod() : data_(0) { EXPECT_TRUE(live.insert(this).second); }
  non_pod(const T& value) : data_(value) { EXPECT_TRUE(live.insert(this).second); }
  non_pod(const non_pod& other) : data_(other.data_) {
    EXPECT_EQ(live.count(&other), 1u);
    EXPECT_TRUE(live.insert(this).second);
  }
  ~non_pod() { EXPECT_EQ(live.erase(this), 1u); }
  non_pod& operator=(const non_pod& other) {
    EXPECT_EQ(live.count(this), 1u);
    EXPECT_EQ(live.count(&other), 1u);
    data_ = other.data_;
    return *this;
  }
  non_pod& operator=(const T& value) {
    EXPECT_EQ(live.count(this), 1u);
    data_ = value;
    return *this;
  }
  bool operator==(const non_pod& other) const { return data_ == other.data_; }
  bool operator!=(const non_pod& other) const { return !(*this == other); }
  friend bool operator==(T value, const non_pod& other) { return value == other.data_; }
  friend bool operator!=(T value, const non_pod& other) { return !(value == other); }
  friend std::ostream& operator<<(std::ostream& os, const non_pod& value) {
    return os << value.data_;
  }
};

template<class Vec, class RNG>
void make_array(Vec& vec, RNG& rng, std::size_t n) {
  vec.clear();
  for (std::size_t i = 0; i < n; ++i) vec.push_back(rng());
}

template<class S, class T, class U>
void test_main(std::size_t m, std::size_t n) {
  S s;
  T t;
  U u;

  std::vector<double> v;
  boost::lagged_fibonacci607 rng(331u); // Preserve the historical default seed.

  for (std::size_t i = 0; i < m; ++i) {
#ifdef VERBOSE
    std::cout << i << ' ';
#endif
    double r = rng();
    if (r < 0.1) {
      // clear
#ifdef VERBOSE
      std::cout << "clear\n";
#endif
      s.clear();
      t.clear();
      u.clear();
    } else if (r < 0.2) {
      // push_back
      double d = rng();
      if (s.size() < n) {
#ifdef VERBOSE
        std::cout << "push_back " << d << std::endl;
#endif
        s.push_back(d);
        t.push_back(d);
        u.push_back(d);
      } else {
#ifdef VERBOSE
        std::cout << std::endl;
#endif
      }
#ifdef DEQUE
    } else if (r < 0.25) {
      // push_front
      double d = rng();
      if (s.size() < n) {
#ifdef VERBOSE
        std::cout << "push_front " << d << std::endl;
#endif
        s.push_front(d);
        t.push_front(d);
        u.push_front(d);
      } else {
#ifdef VERBOSE
        std::cout << std::endl;
#endif
      }
#endif // DEQUE
    } else if (r < 0.3) {
      // pop_back
      if (!s.empty()) {
#ifdef VERBOSE
        std::cout << "pop_back\n";
#endif
        s.pop_back();
        t.pop_back();
        u.pop_back();
      } else {
#ifdef VERBOSE
        std::cout << std::endl;
#endif
      }
#ifdef DEQUE
    } else if (r < 0.35) {
      // pop_front
      if (!s.empty()) {
#ifdef VERBOSE
        std::cout << "pop_front\n";
#endif
        s.pop_front();
        t.pop_front();
        u.pop_front();
      } else {
#ifdef VERBOSE
        std::cout << std::endl;
#endif
      }
#endif // DEQUE
    } else if (r < 0.4) {
      // resize
      double d = rng();
      std::size_t p = int(n * rng());
#ifdef VERBOSE
      std::cout << "resize " << p << ' ' << d << std::endl;
#endif
      s.resize(p, d);
      t.resize(p, d);
      u.resize(p, d);
    } else if (r < 0.7) {
      // insert
      std::size_t p = int(s.size() * rng());
      std::size_t x = int((n - s.size()) * rng());
      if (rng() < 0.5) {
        double d = rng();
#ifdef VERBOSE
        std::cout << "insert " << p << ' ' << x << ' ' << d << std::endl;
#endif
        if (x != 1) {
          s.insert(s.begin() + p, x, d);
          t.insert(t.begin() + p, x, d);
          u.insert(u.begin() + p, x, d);
        } else {
          s.insert(s.begin() + p, d);
          t.insert(t.begin() + p, d);
          u.insert(u.begin() + p, d);
        }
      } else {
#if !(__GNUC__ == 3 && __GNUC_MINOR__ == 1)
        // gcc-3.1 has a bug in std::uninitialized_copy,
        // so just skip this test.
#ifdef VERBOSE
        std::cout << "insert sequence\n";
#endif
        make_array(v, rng, x);
        s.insert(s.begin() + p, v.begin(), v.end());
        t.insert(t.begin() + p, v.begin(), v.end());
        u.insert(u.begin() + p, v.begin(), v.end());
#endif // !(__GNUC__ == 3 && __GNUC_MINOR__ == 1)
      }
    } else {
      // erase
      std::size_t p = int(s.size() * rng());
      std::size_t x = p + int((s.size() - p) * rng());
#ifdef VERBOSE
      std::cout << "erase " << p << ' ' << x << std::endl;
#endif
      if (x != 1) {
        s.erase(s.begin() + p, s.begin() + x);
        t.erase(t.begin() + p, t.begin() + x);
        u.erase(u.begin() + p, u.begin() + x);
      } else {
        s.erase(s.begin() + p);
        t.erase(t.begin() + p);
        u.erase(u.begin() + p);
      }
    }

    SCOPED_TRACE(::testing::Message() << "seed 331, operation " << i << ", draw " << r);
    ASSERT_EQ(s.size(), u.size());
    ASSERT_EQ(t.size(), u.size());
    for (std::size_t index = 0; index < u.size(); ++index) {
      EXPECT_EQ(s[index], u[index]) << "index " << index;
      EXPECT_EQ(t[index].data_, u[index]) << "non-POD index " << index;
    }
  }

  return;
}
