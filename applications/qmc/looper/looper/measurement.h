/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

// boostinspect:nounnamed

#ifndef LOOPER_MEASUREMENT_H
#define LOOPER_MEASUREMENT_H

#include "lattice.h"
#include "power.h"
#include <alps/alea.h>
#include <boost/call_traits.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <string>
#include <vector>

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

struct has_improved_estimator_tag {};
struct has_normal_estimator_tag {};
struct has_pre_evaluator_tag {};
struct has_evaluator_tag {};

//
// helper functions
//

inline
void add_scalar_obs(alps::ObservableSet& m, std::string const& name, bool is_signed = false) {
  if (!m.has(name)) m << make_observable(alps::RealObservable(name), is_signed);
}

inline
void add_vector_obs(alps::ObservableSet& m, std::string const& name, bool is_signed = false) {
  if (!m.has(name)) m << make_observable(alps::RealVectorObservable(name), is_signed);
}

inline
void add_vector_obs(alps::ObservableSet& m, std::string const& name,
  alps::RealVectorObservable::label_type const& label, bool is_signed = false) {
  if (!m.has(name)) m << make_observable(alps::RealVectorObservable(name, label), is_signed);
}

// for path integral
template<typename OP>
inline void proceed(boost::mpl::true_, double& t, OP const& op) { t = op.time(); }
template<typename OP>
inline void proceed(boost::mpl::false_, double, OP const&) {}

// for sse
inline void proceed(boost::mpl::true_, double& t) { t += 1; }
inline void proceed(boost::mpl::false_, double) {}


//
// measurement_set
//

namespace { struct null_measurement; }

template<typename M1,
         typename M2 = null_measurement,
         typename M3 = null_measurement,
         typename M4 = null_measurement,
         typename M5 = null_measurement,
         typename M6 = null_measurement,
         typename M7 = null_measurement,
         typename M8 = null_measurement>
struct measurement_set {};


//
// composite_measurement (forward declaration)
//

template<typename M1, typename M2>
struct composite_measurement;


//
// measurement traits
//

template<typename MEASUREMENT_SET>
struct measurement;

template<typename M1>
struct measurement<measurement_set<M1, null_measurement, null_measurement, null_measurement,
  null_measurement, null_measurement, null_measurement, null_measurement> > {
  typedef M1 type;
};

template<typename M1, typename M2>
struct measurement<measurement_set<M1, M2, null_measurement, null_measurement, null_measurement,
  null_measurement, null_measurement, null_measurement> > {
  typedef composite_measurement<M1, M2> type;
};

template<typename M1, typename M2, typename M3>
struct measurement<measurement_set<M1, M2, M3, null_measurement, null_measurement,
  null_measurement, null_measurement, null_measurement> > {
  typedef composite_measurement<
          composite_measurement<M1, M2>, M3> type;
};

template<typename M1, typename M2, typename M3, typename M4>
struct measurement<measurement_set<M1, M2, M3, M4, null_measurement, null_measurement,
  null_measurement, null_measurement> > {
  typedef composite_measurement<
          composite_measurement<
          composite_measurement<M1, M2>, M3>, M4> type;
};

template<typename M1, typename M2, typename M3, typename M4, typename M5>
struct measurement<measurement_set<M1, M2, M3, M4, M5, null_measurement, null_measurement,
  null_measurement> > {
  typedef composite_measurement<
          composite_measurement<
          composite_measurement<
          composite_measurement<M1, M2>, M3>, M4>, M5> type;
};

template<typename M1, typename M2, typename M3, typename M4, typename M5, typename M6>
struct measurement<measurement_set<M1, M2, M3, M4, M5, M6, null_measurement, null_measurement> > {
  typedef composite_measurement<
          composite_measurement<
          composite_measurement<
          composite_measurement<
          composite_measurement<M1, M2>, M3>, M4>, M5>, M6> type;
};

template<typename M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename M7>
struct measurement<measurement_set<M1, M2, M3, M4, M5, M6, M7, null_measurement> > {
  typedef composite_measurement<
          composite_measurement<
          composite_measurement<
          composite_measurement<
          composite_measurement<
          composite_measurement<M1, M2>, M3>, M4>, M5>, M6>, M7> type;
};

template<typename M1, typename M2, typename M3, typename M4, typename M5, typename M6, typename M7,
  typename M8>
struct measurement<measurement_set<M1, M2, M3, M4, M5, M6, M7, M8> > {
  typedef composite_measurement<
          composite_measurement<
          composite_measurement<
          composite_measurement<
          composite_measurement<
          composite_measurement<
          composite_measurement<M1, M2>, M3>, M4>, M5>, M6>, M7>, M8> type;
};


//
// other traits classes
//

template<typename MEASUREMENT_SET, typename MC, typename LAT, typename TIME>
struct estimator {
  typedef typename measurement<MEASUREMENT_SET>::type measurement_t;
  typedef typename measurement_t::template estimator<MC, LAT, TIME> type;
};

template<typename ESTIMATOR>
struct estimate {
  typedef typename ESTIMATOR::estimate type;
};

template<typename ESTIMATOR>
struct collector {
  template<typename BASE_ESTIMATOR>
  class basic_collector : public  BASE_ESTIMATOR::collector {
  public:
    basic_collector()
      : BASE_ESTIMATOR::collector(), nop_(0), nc_(0), range_(std::make_pair(1, 0)), noc_(0) {}
    basic_collector& operator+=(basic_collector const& coll) {
      BASE_ESTIMATOR::collector::operator+=(coll);
      nop_ += coll.nop_;
      nc_ += coll.nc_;
      range_ = std::make_pair(std::min(range_.first, coll.range_.first),
                              std::max(range_.second, coll.range_.second));
      return *this;
    }
    basic_collector& operator+=(typename BASE_ESTIMATOR::estimate const& est) {
      BASE_ESTIMATOR::collector::operator+=(est);
      return *this;
    }
    void set_num_clusters(unsigned int n) { nc_ = n; }
    void inc_num_clusters(unsigned int n) { nc_ += n; }
    double num_clusters() const { return nc_; }
    void set_num_operators(unsigned int n) { nop_ = n; }
    double num_operators() const { return nop_; }

    // for parallel QMC
    void set_num_open_clusters(unsigned int n) { noc_ = n; }
    unsigned int num_open_clusters() const { return noc_; }
    void clear_range() { range_ = std::make_pair(1, 0); }
    void set_range(int pos) { range_ = std::make_pair(pos, pos); }
    std::pair<int, int> const& range() const { return range_; }
    bool empty() const { return range_.first > range_.second; }

  private:
    double nop_; // total number of operators
    double nc_; // total number of (closed) clusters
    std::pair<int, int> range_; // configuration range (for parallel QMC)
    unsigned int noc_; // number of open clusters (for parallel QMC)
  };
  typedef basic_collector<ESTIMATOR> type;
};

template<typename MEASUREMENT>
struct pre_evaluator_selector {
private:
  template<bool, typename M>
  struct impl {
    static void pre_evaluate(alps::ObservableSet&, alps::Parameters const&,
    alps::ObservableSet const&) {}
  };
  template<typename M>
  struct impl<true, M> {
    static void pre_evaluate(alps::ObservableSet& m, alps::Parameters const& params,
      alps::ObservableSet const& m_in) {
      M::pre_evaluator::pre_evaluate(m, params, m_in);
    }
  };
public:
  static void pre_evaluate(alps::ObservableSet& m, alps::Parameters const& params,
    alps::ObservableSet const& m_in) {
    impl<boost::is_base_of<has_pre_evaluator_tag, MEASUREMENT>::value, MEASUREMENT>::
      pre_evaluate(m, params, m_in);
  }
};

template<typename MEASUREMENT>
struct evaluator_selector {
private:
  template<bool, typename M>
  struct impl {
    static void evaluate(alps::ObservableSet&, alps::Parameters const&,
    alps::ObservableSet const&) {}
  };
  template<typename M>
  struct impl<true, M> {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const& params,
      alps::ObservableSet const& m_in) {
      M::evaluator::evaluate(m, params, m_in);
    }
  };
public:
  static void evaluate(alps::ObservableSet& m, alps::Parameters const& params,
    alps::ObservableSet const& m_in) {
    impl<boost::is_base_of<has_evaluator_tag, MEASUREMENT>::value, MEASUREMENT>::
      evaluate(m, params, m_in);
  }
};


//
// free functions
//

template<typename ESTIMATOR>
typename estimate<ESTIMATOR>::type get_estimate(ESTIMATOR const& estimator) {
  typename estimate<ESTIMATOR>::type estimate;
  estimator.init_estimate(estimate);
  return estimate;
}

template<typename ESTIMATOR>
typename collector<ESTIMATOR>::type get_collector(ESTIMATOR const& estimator) {
  typename collector<ESTIMATOR>::type collector;
  estimator.init_collector(collector);
  return collector;
}


//
// accumulator for improved estimator
//

template<typename ESTIMATOR, typename FRAGMENT, typename IMPROVE>
struct accumulator {
  typedef ESTIMATOR estimator_t;
  typedef typename estimate<estimator_t>::type estimate_t;
  typedef typename estimator_t::lattice_t lattice_t;
  typedef typename boost::call_traits<typename estimator_t::time_t>::param_type time_pt;
  typedef FRAGMENT fragment_t;
  accumulator(std::vector<estimate_t> const&, int /* nc */, lattice_t const& /* lat */,
    estimator_t const& /* emt */, std::vector<fragment_t> const& /* fr */) {}
  void begin_s(int, time_pt, int, int) const {}
  void begin_b(int, int, time_pt, int, int, int, int, int) const {}
  void end_s(int, time_pt, int, int) const {}
  void end_b(int, int, time_pt, int, int, int, int, int) const {}
  void start_bottom(int, time_pt, int, int) const {}
  void start(int, time_pt, int, int) const {}
  void stop_top(int, time_pt, int, int) const {}
  void stop(int, time_pt, int, int) const {}
};

template<typename ESTIMATOR, typename FRAGMENT>
struct accumulator<ESTIMATOR, FRAGMENT, boost::mpl::true_> {
  typedef ESTIMATOR estimator_t;
  typedef typename estimate<estimator_t>::type estimate_t;
  typedef typename estimator_t::lattice_t lattice_t;
  typedef typename boost::call_traits<typename estimator_t::time_t>::param_type time_pt;
  typedef FRAGMENT fragment_t;
  accumulator(std::vector<estimate_t>& es, int nc, lattice_t const& lt, estimator_t const& emt,
    std::vector<fragment_t> const& fr) :
    estimates(es), lat(lt), estimator(emt), fragments(fr) {
    estimates.resize(nc);
    for (int i = 0; i < nc; ++i) estimator.init_estimate(estimates[i]);
  }
  void begin_s(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].begin_s(estimator, lat, t, s, c);
  }
  void begin_b(int p0, int p1, time_pt t, int b, int s0, int s1, int c0, int c1) {
    estimates[fragments[p0].id()].begin_bs(estimator, lat, t, b, s0, c0);
    estimates[fragments[p1].id()].begin_bt(estimator, lat, t, b, s1, c1);
  }
  void end_s(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].end_s(estimator, lat, t, s, c);
  }
  void end_b(int p0, int p1, time_pt t, int b, int s0, int s1, int c0, int c1) {
    estimates[fragments[p0].id()].end_bs(estimator, lat, t, b, s0, c0);
    estimates[fragments[p1].id()].end_bt(estimator, lat, t, b, s1, c1);
  }
  void start_bottom(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].start_bottom(estimator, lat, t, s, c);
  }
  void start(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].start(estimator, lat, t, s, c);
  }
  void stop(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].stop(estimator, lat, t, s, c);
  }
  void stop_top(int p, time_pt t, int s, int c) {
    estimates[fragments[p].id()].stop_top(estimator, lat, t, s, c);
  }

  std::vector<estimate_t>& estimates;
  lattice_t const& lat;
  estimator_t const& estimator;
  std::vector<fragment_t> const& fragments;
};


//
// energy measurement
//

struct energy_estimator {
  template<typename M>
  static void init_observables(M& m, bool is_signed) {
    add_scalar_obs(m, "Energy", is_signed);
    add_scalar_obs(m, "Energy Density", is_signed);
    add_scalar_obs(m, "Energy^2", is_signed);
  }

  // normal estimator

  template<typename RG>
  static void measurement(alps::ObservableSet& m, lattice_helper<RG> const& lat, double beta,
    double nop, double sign, double ene) {
    m["Energy"] << sign * ene;
    m["Energy Density"] << sign * ene / lat.volume();
    m["Energy^2"] << sign * (power2(ene) - nop / power2(beta));
  }
};

struct energy_evaluator {
  static void evaluate(alps::ObservableSet& m, alps::ObservableSet const& m_in) {
    if (m_in.has("Inverse Temperature") && m_in.has("Volume") &&
        m_in.has("Energy") && m_in.has("Energy^2")) {
      alps::RealObsevaluator beta = m_in["Inverse Temperature"];
      alps::RealObsevaluator vol = m_in["Volume"];
      alps::RealObsevaluator ene = m_in["Energy"];
      alps::RealObsevaluator ene2 = m_in["Energy^2"];
      if (beta.count() && vol.count() && ene.count() && ene2.count()) {
        alps::RealObsevaluator c("Specific Heat");
        c = beta.mean() * beta.mean() * (ene2 - ene * ene) / vol.mean();
        m.addObservable(c);
      }
    }
  }
};


//
// dumb measurement
//

template<class DUMMY>
struct dumb_measurement {
  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;
    typedef typename DUMMY::template estimator<MC, LAT, TIME> estimator_t;

    void initialize(alps::Parameters const& /* params */, lattice_t const& /* lat */,
      bool /* is_signed */, bool /* use_improved_estimator */) {}
    template<typename M>
    void init_observables(M& /* m */, bool /* is_signed */) {}

    // improved estimator

    struct estimate {
      estimate& operator+=(estimate const&) { return *this; }
      void begin_s(estimator_t const&, lattice_t const&, double, int, int) const {}
      void begin_bs(estimator_t const&, lattice_t const&, double, int, int, int) const {}
      void begin_bt(estimator_t const&, lattice_t const&, double, int, int, int) const {}
      void end_s(estimator_t const&, lattice_t const&, double, int, int) const {}
      void end_bs(estimator_t const&, lattice_t const&, double, int, int, int) const {}
      void end_bt(estimator_t const&, lattice_t const&, double, int, int, int) const {}
      void start_bottom(estimator_t const&, lattice_t const&, double, int, int) const {}
      void start(estimator_t const&, lattice_t const&, double, int, int) const {}
      void stop(estimator_t const&, lattice_t const&, double, int, int) const {}
      void stop_top(estimator_t const&, lattice_t const&, double, int, int) const {}
    };
    void init_estimate(estimate& /* est */) const {}

    struct collector {
      template<typename EST>
      collector& operator+=(EST const& /* est */) { return *this; }
      template<typename M>
      void commit(M& /* m */, lattice_t const& /* lat */, double /* beta */,
        int /* nop */, double /* sign */) const {}
    };
    void init_collector(collector& /* coll */) const {}

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& /* m */, lattice_t const& /* lat */, double /* beta */,
      double /* sign */, std::vector<int> const& /* spins */,
      std::vector<OP> const& /* operators */, std::vector<int> const& /* spins_c */,
      std::vector<FRAGMENT> const& /* fragments */, collector const& /* coll */) {}

    // normal estimator

    template<typename M, typename OP>
    void normal_measurement(M& /* m */, lattice_t const& /* lat */, double /* beta */,
      double /* sign */, std::vector<int> const& /* spins */,
      std::vector<OP> const& /* operators */, std::vector<int> const& /* spins_c */) {}
  };
};


//
// composite_measurement
//

template<typename MEASUREMENT1, typename MEASUREMENT2>
struct composite_measurement :
  public has_pre_evaluator_tag, public has_evaluator_tag {

  typedef MEASUREMENT1 measurement1;
  typedef MEASUREMENT2 measurement2;

  template<typename MC, typename LAT, typename TIME>
  struct estimator :
    public measurement1::template estimator<MC, LAT, TIME>,
    public measurement2::template estimator<MC, LAT, TIME> {

    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;

    typedef typename boost::call_traits<time_t>::param_type time_pt;

    typedef typename measurement1::template estimator<MC, LAT, TIME> estimator1;
    typedef typename measurement2::template estimator<MC, LAT, TIME> estimator2;
    typedef typename estimator1::estimate estimate1;
    typedef typename estimator2::estimate estimate2;
    typedef typename estimator1::collector collector1;
    typedef typename estimator2::collector collector2;

    void initialize(alps::Parameters const& params, lattice_t const& lat, bool is_signed,
      bool use_improved_estimator) {
      estimator1::initialize(params, lat, is_signed, use_improved_estimator);
      estimator2::initialize(params, lat, is_signed, use_improved_estimator);
    }
    template<typename M>
    void init_observables(M& m, bool is_signed) {
      estimator1::init_observables(m, is_signed);
      estimator2::init_observables(m, is_signed);
    }

    // improved estimator

    struct estimate : public estimate1, public estimate2 {
      estimate() : estimate1(), estimate2() {}
      estimate& operator+=(estimate const& rhs) {
        estimate1::operator+=(rhs);
        estimate2::operator+=(rhs);
        return *this;
      }
      void begin_s(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
        estimate1::begin_s(emt, lat, t, s, c);
        estimate2::begin_s(emt, lat, t, s, c);
      }
      void begin_bs(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s, int c) {
        estimate1::begin_bs(emt, lat, t, b, s, c);
        estimate2::begin_bs(emt, lat, t, b, s, c);
      }
      void begin_bt(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s, int c) {
        estimate1::begin_bt(emt, lat, t, b, s, c);
        estimate2::begin_bt(emt, lat, t, b, s, c);
      }
      void end_s(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
        estimate1::end_s(emt, lat, t, s, c);
        estimate2::end_s(emt, lat, t, s, c);
      }
      void end_bs(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s, int c) {
        estimate1::end_bs(emt, lat, t, b, s, c);
        estimate2::end_bs(emt, lat, t, b, s, c);
      }
      void end_bt(estimator_t const& emt, lattice_t const& lat, time_pt t, int b, int s, int c) {
        estimate1::end_bt(emt, lat, t, b, s, c);
        estimate2::end_bt(emt, lat, t, b, s, c);
      }
      void start_bottom(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
        estimate1::start_bottom(emt, lat, t, s, c);
        estimate2::start_bottom(emt, lat, t, s, c);
      }
      void start(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
        estimate1::start(emt, lat, t, s, c);
        estimate2::start(emt, lat, t, s, c);
      }
      void stop(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
        estimate1::stop(emt, lat, t, s, c);
        estimate2::stop(emt, lat, t, s, c);
      }
      void stop_top(estimator_t const& emt, lattice_t const& lat, time_pt t, int s, int c) {
        estimate1::stop_top(emt, lat, t, s, c);
        estimate2::stop_top(emt, lat, t, s, c);
      }
    };
    void init_estimate(estimate& est) const {
      estimator1::init_estimate(est);
      estimator2::init_estimate(est);
    }

    struct collector : public collector1, public collector2 {
      template<typename EST>
      collector& operator+=(EST const& est) {
        collector1::operator+=(est);
        collector2::operator+=(est);
        return *this;
      }
      template<typename M>
      void commit(M& m, lattice_t const& lat, double beta, int nop, double sign) const {
        collector1::commit(m, lat, beta, nop, sign);
        collector2::commit(m, lat, beta, nop, sign);
      }
    };
    void init_collector(collector& coll) const {
      estimator1::init_collector(coll);
      estimator2::init_collector(coll);
    }

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& spins, std::vector<OP> const& operators,
      std::vector<int> const& spins_c, std::vector<FRAGMENT> const& fragments,
      collector const& coll) {
      estimator1::improved_measurement(m, lat, beta, sign, spins, operators, spins_c, fragments,
                                       coll);
      estimator2::improved_measurement(m, lat, beta, sign, spins, operators, spins_c, fragments,
                                       coll);
    }

    // normal estimator

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& spins, std::vector<OP> const& operators, std::vector<int>& spins_c) {
      estimator1::normal_measurement(m, lat, beta, sign, spins, operators, spins_c);
      estimator2::normal_measurement(m, lat, beta, sign, spins, operators, spins_c);
    }
  };

  struct pre_evaluator {
    static void pre_evaluate(alps::ObservableSet& m, alps::Parameters const& params,
      alps::ObservableSet const& m_in) {
      pre_evaluator_selector<measurement1>::pre_evaluate(m, params, m_in);
      pre_evaluator_selector<measurement2>::pre_evaluate(m, params, m_in);
    }
  };
  struct evaluator {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const& params,
      alps::ObservableSet const& m_in) {
      evaluator_selector<measurement1>::evaluate(m, params, m_in);
      evaluator_selector<measurement2>::evaluate(m, params, m_in);
    }
  };
};

} // end namespace looper

#endif // LOOPER_MEASUREMENT_H
