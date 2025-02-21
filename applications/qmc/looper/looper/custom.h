/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_CUSTOM_H
#define LOOPER_CUSTOM_H

#include "measurement.h"
#include <alps/scheduler/measurement_operators.h>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <vector>

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

template<typename LAT>
struct custom_measurement_initializer : private alps::MeasurementOperators {
  typedef LAT lattice_t;
  typedef boost::tuple<std::string, std::vector<std::vector<double> > > s_elements_type;
  typedef boost::tuple<std::string, std::vector<std::vector<double> >,
    std::vector<std::vector<double> > > p_elements_type;
  typedef std::pair<std::string, std::string> s_expression_type;
  typedef std::pair<std::string, std::pair<std::string, std::string> > p_expression_type;

  custom_measurement_initializer(alps::Parameters const& params);

  void init(lattice_t const& lat,
    std::vector<s_elements_type>& average_elements,
    std::vector<s_elements_type>& local_elements,
    std::vector<p_elements_type>& correlation_elements,
    std::vector<p_elements_type>& strfactor_elements);

private:
  alps::Parameters const& params_;
};

struct custom_measurement {

  typedef dumb_measurement<custom_measurement> dumb;

  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef LAT lattice_t;
    typedef typename custom_measurement_initializer<lattice_t>::s_elements_type s_elements_type;
    typedef typename custom_measurement_initializer<lattice_t>::p_elements_type p_elements_type;
    typedef typename alps::property_map<coordinate_t, const typename lattice_t::real_graph_type,
              coordinate_type>::type coordinate_map_t;

    std::vector<s_elements_type> average_elements;
    std::vector<s_elements_type> local_elements;
    std::vector<p_elements_type> correlation_elements;
    std::vector<p_elements_type> strfactor_elements;
    coordinate_map_t coordinate;
    boost::optional<int> origin;
    std::valarray<double> mltplcty;
    std::vector<std::string> site_label, distance_label, momenta_label;

    // working vector
    std::valarray<double> local, corr, sfac;

    void initialize(alps::Parameters const& params, lattice_t const& lat,
      bool /* is_signed */, bool /* use_improved_estimator */) {
      custom_measurement_initializer<lattice_t> initializer(params);
      initializer.init(lat, average_elements, local_elements, correlation_elements,
                       strfactor_elements);
      if (local_elements.size()) {
        site_label = alps::site_labels(lat.rg());
        local.resize(site_label.size());
      }
      if (correlation_elements.size()) {
        if (params.defined("INITIAL_SITE"))
          origin = static_cast<int>(params["INITIAL_SITE"]);
        distance_label = distance_labels(lat, origin);
        if (!origin) {
          std::vector<unsigned int> m = distance_multiplicities(lat);
          mltplcty.resize(m.size());
          for (int i = 0; i < mltplcty.size(); ++i) mltplcty[i] = 1. / m[i];
        }
        corr.resize(distance_label.size());
      }
      if (strfactor_elements.size()) {
        coordinate = alps::get_or_default(coordinate_t(), lat.rg(), 0);
        momenta_label = momenta_labels(lat);
        sfac.resize(momenta_label.size());
      }
    }
    template<typename M>
    void init_observables(M& m, bool is_signed) {
      BOOST_FOREACH(s_elements_type const& elms, average_elements)
        add_scalar_obs(m, boost::get<0>(elms), is_signed);
      if (local_elements.size()) {
        BOOST_FOREACH(s_elements_type const& elms, local_elements)
          add_vector_obs(m, boost::get<0>(elms), site_label, is_signed);
      }
      if (correlation_elements.size()) {
        BOOST_FOREACH(p_elements_type const& elms, correlation_elements)
          add_vector_obs(m, boost::get<0>(elms), distance_label, is_signed);
      }
      if (strfactor_elements.size()) {
        BOOST_FOREACH(p_elements_type const& elms, strfactor_elements)
          add_vector_obs(m, boost::get<0>(elms), momenta_label, is_signed);
      }
    }

    typedef typename dumb::template estimator<MC, LAT, TIME>::estimate estimate;
    void init_estimate(estimate const&) const {}

    typedef typename dumb::template estimator<MC, LAT, TIME>::collector collector;
    void init_collector(collector const&) const {}

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& /* m */, lattice_t const& /* lat */, double /* beta */,
      double /* sign */, std::vector<int> const& /* spins */,
      std::vector<OP> const& /* operators */, std::vector<int> const& /* spins_c */,
      std::vector<FRAGMENT> const& /* fragments */, collector const& /* coll */) {}

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat, double /* beta */, double sign,
      std::vector<int> const& spins, std::vector<OP> const& /* operators */,
      std::vector<int>& /* spins_c */) {

      // average
      BOOST_FOREACH(s_elements_type const& elms, average_elements) {
        double v = 0;
        BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs,
          sites(lat.rg())) {
          int st = 0;
          BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs,
            sites(lat, rs))
            st += (1 - spins[vs]);
          v += boost::get<1>(elms)[get(site_type_t(), lat.rg(), rs)][st];
        }
        m[boost::get<0>(elms)] << sign * v / lat.volume();
      }

      // local
      BOOST_FOREACH(s_elements_type const& elms, local_elements) {
        local = 0;
        BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs,
          sites(lat.rg())) {
          int st = 0;
          BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs,
            sites(lat, rs))
            st += (1 - spins[vs]);
          local[rs] = sign * boost::get<1>(elms)[get(site_type_t(), lat.rg(), rs)][st];
        }
        m[boost::get<0>(elms)] << local;
      }

      // correlation
      BOOST_FOREACH(p_elements_type const& elms, correlation_elements) {
        corr = 0;
        if (origin) {
          int st0 = 0;
          BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs0,
            sites(lat, origin.get()))
            st0 += (1 - spins[vs0]);
          double v0 = boost::get<1>(elms)[get(site_type_t(), lat.rg(), origin.get())][st0];
          BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs1,
            sites(lat.rg())) {
            int st1 = 0;
            BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs1,
              sites(lat, rs1))
              st1 += (1 - spins[vs1]);
            double v1 = boost::get<2>(elms)[get(site_type_t(), lat.rg(), rs1)][st1];
            corr[rs1] = sign * v0 * v1;
          }
        } else {
          BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs0,
            sites(lat.rg())) {
            int st0 = 0;
            BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs0,
              sites(lat, rs0))
              st0 += (1 - spins[vs0]);
            double v0 = boost::get<1>(elms)[get(site_type_t(), lat.rg(), rs0)][st0];
            BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs1,
              sites(lat.rg())) {
              int st1 = 0;
              BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs1,
                sites(lat, rs1))
                st1 += (1 - spins[vs1]);
              double v1 = boost::get<2>(elms)[get(site_type_t(), lat.rg(), rs1)][st1];
              int d = distance(lat, rs0, rs1);
              corr[d] += v0 * v1;
            }
          }
          corr *= sign * mltplcty;
        }
        m[boost::get<0>(elms)] << corr;
      }

      // structure factor
      BOOST_FOREACH(p_elements_type const& elms, strfactor_elements) {
        sfac = 0;
        int k = 0;
        typename momentum_iterator<lattice_t>::type mit, mit_end;
        for (boost::tie(mit, mit_end) = momenta(lat); mit != mit_end; ++mit, ++k) {
          std::complex<double> v0, v1;
          BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& rs,
            sites(lat.rg())) {
            int t = get(site_type_t(), lat.rg(), rs);
            int st = 0;
            BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type const& vs,
              sites(lat, rs))
              st += (1 - spins[vs]);
            double phase = alps::numeric::scalar_product(momentum(*mit,lat), coordinate(rs));
            std::complex<double> cphase(std::cos(phase), std::sin(phase));
            v0 += boost::get<1>(elms)[t][st] * cphase;
            v1 += boost::get<2>(elms)[t][st] * cphase;
          }
          sfac[k] = std::real(std::conj(v0) * v1);
        }
        sfac *= (sign / lat.volume());
        m[boost::get<0>(elms)] << sfac;
      }
    }
  };
};

} // end namespace looper

#endif // LOOPER_CUSTOM_H
