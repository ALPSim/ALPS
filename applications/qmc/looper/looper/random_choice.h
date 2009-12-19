/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_RANDOM_CHOICE_H
#define LOOPER_RANDOM_CHOICE_H

#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>
#include <boost/throw_exception.hpp>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace looper {

//
// double-based Walker algorithm
//

template<class IntType = unsigned int, class RealType = double>
class random_choice_walker_d {
public:
  typedef RealType input_type;
  typedef IntType result_type;

  random_choice_walker_d() {}
  template<class CONT>
  random_choice_walker_d(const CONT& weights) { init(weights); }

  template<class Engine>
  result_type operator()(Engine& eng) const {
    result_type x = result_type(RealType(size()) * eng());
    return (eng() < cutoff(x)) ? x : alias(x);
  }

  // Initialization routine with complexity O(N).
  template<class CONT>
  void init(const CONT& weights) {
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    BOOST_STATIC_ASSERT(std::numeric_limits<IntType>::is_integer);
    BOOST_STATIC_ASSERT(!std::numeric_limits<RealType>::is_integer);
#endif
    if (weights.size() == 0)
      boost::throw_exception(std::invalid_argument("random_choice_walker_d::init"));
    IntType n = weights.size();
    RealType norm = RealType(0);
    BOOST_FOREACH(RealType w, weights) {
      if (w < RealType(0))
        boost::throw_exception(std::invalid_argument("random_choice_walker_d::init"));
      norm += w;
    }
    if (norm <= RealType(0))
      boost::throw_exception(std::invalid_argument("random_choice_walker_d::init"));
    norm = n / norm;

    // Initialize arrays.  We will reorder the elements in `array', so
    // that all the negative elements precede the positive ones.
    resize(n);
    std::vector<std::pair<RealType, result_type> > array(n);
    typename std::vector<std::pair<RealType, result_type> >::iterator
      neg_p = array.begin();
    typename std::vector<std::pair<RealType, result_type> >::iterator
      pos_p = array.end();
    for (result_type i = 0; i < n; ++i) {
      RealType b = norm * weights[i] - RealType(1.);
      if (b < RealType(0.)) {
        *neg_p = std::make_pair(b, i);
        ++neg_p;
      } else {
        --pos_p;
        *pos_p = std::make_pair(b, i);
      }
    }

    // Note: now `pos_p' is pointing the first non-negative element in
    // the array.

    // Assign alias and cutoff values
    for (neg_p = array.begin(); neg_p != array.end(); ++neg_p) {
      if (pos_p != array.end()) {
        cutoff(neg_p->second) = RealType(1.) + neg_p->first;
        alias(neg_p->second) = pos_p->second;
        pos_p->first += neg_p->first;
        if (pos_p->first <= RealType(0.)) ++pos_p;
      } else {
        cutoff(neg_p->second) = RealType(1.);
        alias(neg_p->second) = neg_p->second; // never referred
      }
    }
  }

  // Original O(N^2) initialization routine given in A. W. Walker, ACM
  // Trans. Math. Software, 3, 253 (1077).
  template<class CONT>
  void init_walker1977(const CONT& weights, RealType tol = 1.0e-10) {
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    BOOST_STATIC_ASSERT(std::numeric_limits<IntType>::is_integer);
    BOOST_STATIC_ASSERT(!std::numeric_limits<RealType>::is_integer);
#endif
    if (weights.size() == 0)
      boost::throw_exception(std::invalid_argument("random_choice_walker_d::init_walker1977"));
    IntType n = weights.size();
    RealType norm = RealType(0);
    BOOST_FOREACH(RealType w, weights) {
      if (w < RealType(0))
        boost::throw_exception(std::invalid_argument("random_choice_walker_d::init_walker1977"));
      norm += w;
    }
    if (norm <= RealType(0))
      boost::throw_exception(std::invalid_argument("random_choice_walker_d::init_walker1977"));
    norm = n / norm;

    // Initialize arrays.
    resize(n);
    std::vector<RealType> b(n);
    for (result_type i = 0; i < n; ++i) {
      cutoff(i) = 1;
      alias(i) = i;
      b[i] = norm * weights[i] - 1;
    }

    for (result_type i = 0; i < n; ++i) {
      // Find the largest positive and negative differences and their
      // positions.
      RealType sum = 0;
      RealType minval = 0;
      RealType maxval = 0;
      result_type minpos, maxpos;
      for (result_type j = 0; j < n; ++j) {
        sum += std::abs(b[j]);
        if (b[j] <= minval) {
          minval = b[j];
          minpos = j;
        }
        if (b[j] >= maxval) {
          maxval = b[j];
          maxpos = j;
        }
      }

      if (sum < tol) break;

      // Assign alias and cutoff values
      cutoff(minpos) = 1 + minval;
      alias(minpos) = maxpos;
      b[maxpos] += minval;
      b[minpos] = 0;
    }
  }

  template<class CONT>
  bool check(const CONT& weights, RealType tol = 1.0e-10) const {
    bool r = true;
    tol *= size();

    RealType norm = RealType(0);
    BOOST_FOREACH(RealType w, weights) norm += w;
    norm = weights.size() / norm;

    for (result_type i = 0; i < size(); ++i) {
      RealType p = cutoff(i);
      for (result_type j = 0; j < size(); ++j)
        if (alias(j) == i) p += 1 - cutoff(j);
      if (std::abs(p - norm * weights[i]) > tol) r = false;
    }
    return r;
  }

protected:
  IntType size() const { return table_.size(); }
  void resize(IntType n) {
    table_.clear();
    table_.resize(n);
  }

  RealType& cutoff(result_type i) { return table_[i].first; }
  RealType cutoff(result_type i) const { return table_[i].first; }

  result_type& alias(result_type i) { return table_[i].second; }
  result_type alias(result_type i) const { return table_[i].second; }

private:
  std::vector<std::pair<RealType, result_type> > table_; // first element:  cutoff value
                                                         // second element: alias
};


//
// optimized integer-based version of Walker algorithm
//

template<class IntType = unsigned int, class RealType = double>
class random_choice_walker_i {
public:
  typedef IntType input_type;
  typedef IntType result_type;

  random_choice_walker_i() {}
  template<class CONT>
  random_choice_walker_i(const CONT& weights) { init(weights); }

  template<class Engine>
  result_type operator()(Engine& eng) const {
    result_type x = eng() >> bits_;
    return (eng() < cutoff(x)) ? x : alias(x);
  }

  template<class CONT>
  void init(const CONT& weights) {
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    BOOST_STATIC_ASSERT(std::numeric_limits<IntType>::is_integer);
    BOOST_STATIC_ASSERT(std::numeric_limits<IntType>::digits == 32);
    BOOST_STATIC_ASSERT(!std::numeric_limits<IntType>::is_signed);
    BOOST_STATIC_ASSERT(!std::numeric_limits<RealType>::is_integer);
#else
    if (std::numeric_limits<IntType>::min() != 0 ||
        std::numeric_limits<IntType>::max() != 4294967295ul)
      boost::throw_exception(std::range_error("random_choice_walker_i::init"));
#endif
    if (weights.size() == 0)
      boost::throw_exception(std::range_error("random_choice_walker_i::init"));
    IntType n = 2;
    bits_ = 31;
    while (n < weights.size()) {
      n <<= 1;
      bits_ -= 1;
    }

    double norm = 0;
    for (result_type i = 0; i < weights.size(); ++i) norm += weights[i];
    norm = n / norm;

    std::vector<double> w(n, 0);
    std::copy(weights.begin(), weights.end(), w.begin());

    // Initialize arrays.  We will reorder the elements in `array', so
    // that all the negative elements precede the positive ones.
    resize(n);
    std::vector<std::pair<RealType, result_type> > array(n);
    typename std::vector<std::pair<RealType, result_type> >::iterator
      neg_p = array.begin();
    typename std::vector<std::pair<RealType, result_type> >::iterator
      pos_p = array.end();
    for (result_type i = 0; i < n; ++i) {
      RealType b = norm * w[i] - RealType(1.);
      if (b < RealType(0.)) {
        *neg_p = std::make_pair(b, i);
        ++neg_p;
      } else {
        --pos_p;
        *pos_p = std::make_pair(b, i);
      }
    }

    // Note: now `pos_p' is pointing the first non-negative element in
    // the array.

    // Assign alias and cutoff values
    const double nm = std::numeric_limits<IntType>::max();
    for (neg_p = array.begin(); neg_p != array.end(); ++neg_p) {
      if (pos_p != array.end()) {
        cutoff(neg_p->second) = IntType(nm * (RealType(1.) + neg_p->first));
        alias(neg_p->second) = pos_p->second;
        pos_p->first += neg_p->first;
        if (pos_p->first <= RealType(0.)) ++pos_p;
      } else {
        cutoff(neg_p->second) = IntType(nm);
        alias(neg_p->second) = neg_p->second; // never referred
      }
    }
  }

  template<class CONT>
  bool check(const CONT& weights, RealType tol = 1.0e-10) const {
    bool r = true;
    tol *= table_.size();

    RealType norm = RealType(0);
    BOOST_FOREACH(RealType w, weights) norm += w;
    norm = weights.size() / norm;

    const double nm = RealType(1) / std::numeric_limits<IntType>::max();
    for (result_type i = 0; i < table_.size(); ++i) {
      RealType p = nm * cutoff(i);
      for (result_type j = 0; j < table_.size(); ++j)
        if (alias(j) == i) p += 1 - nm * cutoff(j);
      if (std::abs(p - norm * weights[i]) > tol) r = false;
    }
    return r;
  }

protected:
  void resize(IntType n) { table_.resize(n); }

  IntType& cutoff(IntType i) { return table_[i].first; }
  IntType cutoff(IntType i) const { return table_[i].first; }

  IntType& alias(IntType i) { return table_[i].second; }
  IntType alias(IntType i) const { return table_[i].second; }

private:
  IntType bits_; // number of bits to be disposed
  std::vector<std::pair<IntType, IntType> > table_;
};


//
// random_choice_bsearch (O(log N) algorithm using binary search algorithm)
//

template<class IntType = unsigned int, class RealType = double>
class random_choice_bsearch {
public:
  typedef RealType input_type;
  typedef IntType result_type;

  random_choice_bsearch() : accum_(0) {}
  template<class CONT>
  random_choice_bsearch(CONT const& weights) { init(weights); }

  template<class CONT>
  void init(CONT const& weights) {
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    BOOST_STATIC_ASSERT(std::numeric_limits<IntType>::is_integer);
    BOOST_STATIC_ASSERT(!std::numeric_limits<RealType>::is_integer);
#endif
    if (weights.size() == 0)
      boost::throw_exception(std::invalid_argument("random_choice_bsearch::init"));
    RealType norm = 0;
    BOOST_FOREACH(RealType w, weights) {
      if (w < RealType(0))
        boost::throw_exception(std::invalid_argument("random_choice_bsearch::init"));
      norm += w;
    }
    if (norm <= RealType(0))
      boost::throw_exception(std::invalid_argument("random_choice_bsearch::init"));
    accum_.resize(0);
    double a = 0;
    BOOST_FOREACH(RealType w, weights) {
      a += w / norm;
      accum_.push_back(a);
    }
  }

  template<class Engine>
  result_type operator()(Engine& eng) const {
    double p = eng();
    int first = 0;
    int last = accum_.size();  // pointing to the next of the last element
    int current = first + ((last - first) >> 1);
    if (last - first == 0) {
      return first;
    } else if (last - first == 1) {
      if (p < accum_[first])
        return first;
      else
        return first + 1;
    }
    while (true) {
      if (last - first > 3) {
        if (p < accum_[current]) {
          last = current;
          current = first + ((current - first) >> 1);
        } else {
          first = current;
          current += ((last - current) >> 1);
        }
      } else if (last - first == 3) {
        if (p < accum_[first])
          return first;
        else if  (p < accum_[first + 1])
          return first + 1;
        else if  (p < accum_[first + 2])
          return first + 2;
        else
          return first + 3;
      } else {
        if (last - first == 2) {
          if (p < accum_[first])
            return first;
          else if  (p < accum_[first + 1])
            return first + 1;
          else
            return first + 2;
        }
      }
    }
  }

private:
  std::vector<RealType> accum_;
};


//
// random_choice_lsearch (O(N) algorithm with naive linear search)
//

template<class IntType = unsigned int, class RealType = double>
class random_choice_lsearch {
public:
  typedef RealType input_type;
  typedef IntType result_type;

  random_choice_lsearch() : accum_(0) {}
  template<class CONT>
  random_choice_lsearch(CONT const& weights) { init(weights); }

  template<class CONT>
  void init(CONT const& weights) {
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    BOOST_STATIC_ASSERT(std::numeric_limits<IntType>::is_integer);
    BOOST_STATIC_ASSERT(!std::numeric_limits<RealType>::is_integer);
#endif
    if (weights.size() == 0)
      boost::throw_exception(std::invalid_argument("random_choice_lsearch::init"));
    double norm = 0;
    BOOST_FOREACH(RealType w, weights) {
      if (w < RealType(0))
        boost::throw_exception(std::invalid_argument("random_choice_lsearch::init"));
      norm += w;
    }
    if (norm <= RealType(0))
      boost::throw_exception(std::invalid_argument("random_choice_lsearch::init"));
    accum_.resize(0);
    double a = 0;
    BOOST_FOREACH(RealType w, weights) {
      a += w / norm;
      accum_.push_back(a);
    }
  }

  template<class Engine>
  result_type operator()(Engine& eng) const {
    RealType x = eng();
    for (result_type r = 0; r < accum_.size(); ++r) if (accum_[r] > x) return r;
    return result_type(0); // never reached
  }

private:
  std::vector<RealType> accum_;
};

typedef random_choice_walker_i<unsigned int, double> random_choice;

} // end namespace looper

#endif // LOOPER_RANDOM_CHOICE_H
