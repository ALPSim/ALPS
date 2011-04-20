/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_SUSCEPTIBILITY_H
#define LOOPER_SUSCEPTIBILITY_H

#include "divide_if_positive.h"
#include "measurement.h"
#include "type.h"

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

struct susceptibility :
  public has_evaluator_tag {
  template<typename MC, typename LAT, typename TIME>
  struct estimator {
    typedef MC   mc_type;
    typedef LAT lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;
    typedef typename alps::property_map<gauge_t,
              const typename lattice_t::virtual_graph_type,
              double>::type gauge_map_t;

    bool bipartite, improved;
    gauge_map_t gauge;

    void initialize(alps::Parameters const& /* params */, lattice_t const& lat,
      bool /* is_signed */, bool use_improved_estimator) {
      bipartite = is_bipartite(lat);
      improved = use_improved_estimator;
      gauge = alps::get_or_default(gauge_t(), lat.vg(), 0);
    }
    template<typename M>
    void init_observables(M& m, bool is_signed) {
      add_scalar_obs(m, "Magnetization", is_signed);
      add_scalar_obs(m, "Magnetization Density", is_signed);
      add_scalar_obs(m, "|Magnetization|", is_signed);
      add_scalar_obs(m, "|Magnetization Density|", is_signed);
      add_scalar_obs(m, "Magnetization^2", is_signed);
      add_scalar_obs(m, "Magnetization Density^2", is_signed);
      add_scalar_obs(m, "Magnetization^4", is_signed);
      add_scalar_obs(m, "Magnetization Density^4", is_signed);
      add_scalar_obs(m, "Susceptibility", is_signed);
      if (improved) {
        add_scalar_obs(m, "Generalized Magnetization^2", is_signed);
        add_scalar_obs(m, "Generalized Magnetization Density^2", is_signed);
        add_scalar_obs(m, "Generalized Magnetization^4", is_signed);
        add_scalar_obs(m, "Generalized Magnetization Density^4", is_signed);
        add_scalar_obs(m, "Generalized Susceptibility", is_signed);
      }
      if (bipartite) {
        add_scalar_obs(m, "Staggered Magnetization", is_signed);
        add_scalar_obs(m, "Staggered Magnetization Density", is_signed);
        add_scalar_obs(m, "|Staggered Magnetization|", is_signed);
        add_scalar_obs(m, "|Staggered Magnetization Density|", is_signed);
        add_scalar_obs(m, "Staggered Magnetization^2", is_signed);
        add_scalar_obs(m, "Staggered Magnetization Density^2", is_signed);
        add_scalar_obs(m, "Staggered Magnetization^4", is_signed);
        add_scalar_obs(m, "Staggered Magnetization Density^4", is_signed);
        add_scalar_obs(m, "Staggered Susceptibility", is_signed);
        if (improved) {
          add_scalar_obs(m, "Generalized Staggered Magnetization^2", is_signed);
          add_scalar_obs(m, "Generalized Staggered Magnetization Density^2", is_signed);
          add_scalar_obs(m, "Generalized Staggered Magnetization^4", is_signed);
          add_scalar_obs(m, "Generalized Staggered Magnetization Density^4", is_signed);
          add_scalar_obs(m, "Generalized Staggered Susceptibility", is_signed);
        }
      }
    }

    // improved estimator

    struct estimate {
      double usize0, umag0, usize, umag;
      double ssize0, smag0, ssize, smag;
      estimate() : usize0(0), umag0(0), usize(0), umag(0), ssize0(0), smag0(0), ssize(0),
        smag(0) {}
      void init() {
        usize0 = umag0 = usize = umag = 0;
        ssize0 = smag0 = ssize = smag = 0;
      }
      estimate& operator+=(estimate const& rhs) {
        usize0 += rhs.usize0;
        umag0 += rhs.umag0;
        usize += rhs.usize;
        umag += rhs.umag;
        ssize0 += rhs.ssize0;
        smag0 += rhs.smag0;
        ssize += rhs.ssize;
        smag += rhs.smag;
        return *this;
      }
      void begin_s(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
        end_s(emt, lat, -t, s, c);
      }
      void begin_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
        begin_s(emt, lat, t, s, c);
      }
      void begin_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
        begin_s(emt, lat, t, s, c);
      }
      void end_s(estimator_t const& emt, lattice_t const&, double t, int s, int c) {
        usize += t * 0.5;
        umag  += t * (0.5-c);
        double gg = emt.gauge[s];
        ssize += gg * t * 0.5;
        smag  += gg * t * (0.5-c);
      }
      void end_bs(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
        end_s(emt, lat, t, s, c); }
      void end_bt(estimator_t const& emt, lattice_t const& lat, double t, int, int s, int c) {
        end_s(emt, lat, t, s, c); }
      void start_bottom(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
        begin_s(emt, lat, t, s, c);
        usize0 += 0.5;
        umag0  += (0.5-c);
        double gg = emt.gauge[s];
        ssize0 += gg * 0.5;
        smag0  += gg * (0.5-c);
      }
      void start(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
        begin_s(emt, lat, t, s, c);
      }
      void stop(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
        end_s(emt, lat, t, s, c);
      }
      void stop_top(estimator_t const& emt, lattice_t const& lat, double t, int s, int c) {
        end_s(emt, lat, t, s, c);
      }
    };
    void init_estimate(estimate& est) const { est.init(); }

    struct collector {
      double umag0, usize2, umag2, usize4, umag4, usize, umag;
      double smag0, ssize2, smag2, ssize4, smag4, ssize, smag;
      void init() {
        umag0 = usize2 = umag2 = usize4 = umag4 = usize = umag = 0;
        smag0 = ssize2 = smag2 = ssize4 = smag4 = ssize = smag = 0;
      }
      collector& operator+=(collector const& coll) {
        umag0 += coll.umag0;
        usize2 += coll.usize2;
        umag2  += coll.umag2;
        usize4 += coll.usize4;
        umag4  += coll.umag4;
        usize  += coll.usize;
        umag   += coll.umag;
        smag0 += coll.smag0;
        ssize2 += coll.ssize2;
        smag2  += coll.smag2;
        ssize4 += coll.ssize4;
        smag4  += coll.smag4;
        ssize  += coll.ssize;
        smag   += coll.smag;
        return *this;
      }
      collector& operator+=(estimate const& est) {
        umag0  += est.umag0;
        usize2 += power2(est.usize0);
        umag2  += power2(est.umag0);
        usize4 += power4(est.usize0);
        umag4  += power4(est.umag0);
        usize  += power2(est.usize);
        umag   += power2(est.umag);
        smag0  += est.smag0;
        ssize2 += power2(est.ssize0);
        smag2  += power2(est.smag0);
        ssize4 += power4(est.ssize0);
        smag4  += power4(est.smag0);
        ssize  += power2(est.ssize);
        smag   += power2(est.smag);
        return *this;
      }
      template<typename M>
      void commit(M& m, lattice_t const& lat, double beta, int nop, double sign) const {
        double vol = lat.volume();
        m["Magnetization"] << 0.0;
        m["Magnetization Density"] << 0.0;
        m["|Magnetization|"] << sign * std::abs(umag0);
        m["|Magnetization Density|"] << sign * std::abs(umag0) / vol;
        m["Magnetization^2"] << sign * umag2;
        m["Magnetization Density^2"] << sign * umag2 / power2(vol);
        m["Magnetization^4"] << sign * (3 * power2(umag2) - 2 * umag4);
        m["Magnetization Density^4"]
          << sign * (3 * power2(umag2) - 2 * umag4) / power4(vol);
        m["Susceptibility"]
          << (typename is_sse<mc_type>::type() ?
              sign * beta * (dip(umag, nop) + umag2) / (nop + 1) / vol :
              sign * beta * umag / vol);
        m["Generalized Magnetization^2"] << sign * usize2;
        m["Generalized Magnetization Density^2"]
          << sign * usize2 / power2(vol);
        m["Generalized Magnetization^4"]
          << sign * (3 * power2(usize2) - 2 * usize4);
        m["Generalized Magnetization Density^4"]
          << sign * (3 * power2(usize2) - 2 * usize4) / power4(vol);
        m["Generalized Susceptibility"]
          << (typename is_sse<mc_type>::type() ?
              sign * beta * (dip(usize, nop) + usize2) / (nop + 1) / vol :
              sign * beta * usize / vol);
        if (is_bipartite(lat)) {
          m["Staggered Magnetization"] << 0.0;
          m["Staggered Magnetization Density"] << 0.0;
          m["|Staggered Magnetization|"] << sign * std::abs(smag0);
          m["|Staggered Magnetization Density|"] << sign * std::abs(smag0) / vol;
          m["Staggered Magnetization^2"] << sign * smag2;
          m["Staggered Magnetization Density^2"] << sign * smag2 / power2(vol);
          m["Staggered Magnetization^4"]
            << sign * (3 * power2(smag2) - 2 * smag4);
          m["Staggered Magnetization Density^4"]
            << sign * (3 * power2(smag2) - 2 * smag4) / power4(vol);
          m["Staggered Susceptibility"]
            << (typename is_sse<mc_type>::type() ?
                sign * beta * (dip(smag, nop) + smag2) / (nop + 1) / vol :
                sign * beta * smag / vol);
          m["Generalized Staggered Magnetization^2"] << sign * ssize2;
          m["Generalized Staggered Magnetization Density^2"]
            << sign * ssize2 / power2(vol);
          m["Generalized Staggered Magnetization^4"]
            << sign * (3 * power2(ssize2) - 2 * ssize4);
          m["Generalized Staggered Magnetization Density^4"]
            << sign * (3 * power2(ssize2) - 2 * ssize4) / power4(vol);
          m["Generalized Staggered Susceptibility"]
            << (typename is_sse<mc_type>::type() ?
                sign * beta * (dip(ssize, nop) + ssize2) / (nop + 1) / vol :
                sign * beta * ssize / vol);
        }
      }
    };
    void init_collector(collector& coll) const { coll.init(); }

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& /* spins */, std::vector<OP> const& operators,
      std::vector<int> const& /* spins_c */, std::vector<FRAGMENT> const& /* fragments */,
      collector const& coll) {
      coll.commit(m, lat, beta, operators.size(), sign);
    }

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat, double beta, double sign,
      std::vector<int> const& spins, std::vector<OP> const& operators,
      std::vector<int>& spins_c) {
      if (improved) return;

      double vol = lat.volume();
      double nop = operators.size();
      double umag = 0;
      double smag = 0;
      typename virtual_site_iterator<lattice_t>::type si, si_end;
      for (boost::tie(si, si_end) = sites(lat.vg()); si != si_end; ++si) {
        umag += 0.5-spins[*si];
        smag += (0.5-spins[*si]) * gauge[*si];
      }
      m["Magnetization"] << sign * umag;
      m["Magnetization Density"] << sign * umag / vol;
      m["|Magnetization|"] << sign * std::abs(umag);
      m["|Magnetization Density|"] << sign * std::abs(umag) / vol;
      m["Magnetization^2"] << sign * power2(umag);
      m["Magnetization Density^2"] << sign * power2(umag / vol);
      m["Magnetization^4"] << sign * power4(umag);
      m["Magnetization Density^4"] << sign * power4(umag / vol);
      if (is_bipartite(lat)) {
        m["Staggered Magnetization"] << sign * smag;
        m["Staggered Magnetization Density"] << sign * smag / vol;
        m["|Staggered Magnetization|"] << sign * std::abs(smag);
        m["|Staggered Magnetization Density|"] << sign * std::abs(smag) / vol;
        m["Staggered Magnetization^2"] << sign * power2(smag);
        m["Staggered Magnetization Density^2"] << sign * power2(smag / vol);
        m["Staggered Magnetization^4"] << sign * power4(smag);
        m["Staggered Magnetization Density^4"] << sign * power4(smag / vol);
      }
      double umag_a = 0; /* 0 * umag; */
      double smag_a = 0; /* 0 * smag; */
      std::copy(spins.begin(), spins.end(), spins_c.begin());
      double t = 0;
      for (typename std::vector<OP>::const_iterator oi = operators.begin();
           oi != operators.end(); ++oi) {
        if (oi->is_offdiagonal()) {
          proceed(typename is_path_integral<mc_type>::type(), t, *oi);
          umag_a += t * umag;
          smag_a += t * smag;
          if (oi->is_site()) {
            unsigned int s = oi->pos();
            spins_c[s] ^= 1;
            umag += 1-2*spins_c[s];
            smag += gauge[s] * (1-2*spins_c[s]);
          } else {
            unsigned int s0 = source(oi->pos(), lat.vg());
            unsigned int s1 = target(oi->pos(), lat.vg());
            spins_c[s0] ^= 1;
            spins_c[s1] ^= 1;
            umag += 1-2*spins_c[s0] + 1-2*spins_c[s1];
            smag += gauge[s0] * (1-2*spins_c[s0])
              + gauge[s1] * (1-2*spins_c[s1]);
          }
          umag_a -= t * umag;
          smag_a -= t * smag;
        }
        proceed(typename is_sse<mc_type>::type(), t);
      }
      if (typename is_path_integral<mc_type>::type()) {
        umag_a += umag;
        m["Susceptibility"] << sign * beta * power2(umag_a) / vol;
        smag_a += smag;
        if (bipartite)
          m["Staggered Susceptibility"] << sign * beta * power2(smag_a) / vol;
      } else {
        umag_a += nop * umag;
        m["Susceptibility"]
          << sign * beta * (dip(power2(umag_a), nop) + power2(umag)) / (nop + 1) / vol;
        smag_a += nop * smag;
        if (bipartite)
          m["Staggered Susceptibility"]
            << sign * beta * (dip(power2(smag_a), nop) + power2(smag)) / (nop + 1) / vol;
      }
    }
  };

  struct evaluator {
    static void evaluate(alps::ObservableSet& m, alps::Parameters const& /* params */,
      alps::ObservableSet const& m_in) {
      if (m_in.has("Magnetization^2") && m_in.has("Magnetization^4")) {
        alps::RealObsevaluator obse_m2 = m_in["Magnetization^2"];
        alps::RealObsevaluator obse_m4 = m_in["Magnetization^4"];
        if (obse_m2.count() && obse_m4.count()) {
          alps::RealObsevaluator eval("Binder Ratio of Magnetization");
          eval = power2(obse_m2) / obse_m4;
          m.addObservable(eval);
        }
      }
      if (m_in.has("Staggered Magnetization^2") &&
          m_in.has("Staggered Magnetization^4")) {
        alps::RealObsevaluator obse_m2 = m_in["Staggered Magnetization^2"];
        alps::RealObsevaluator obse_m4 = m_in["Staggered Magnetization^4"];
        if (obse_m2.count() && obse_m4.count()) {
          alps::RealObsevaluator eval("Binder Ratio of Staggered Magnetization");
          eval = power2(obse_m2) / obse_m4;
          m.addObservable(eval);
        }
      }
      if (m_in.has("Generalized Magnetization^2") &&
          m_in.has("Generalized Magnetization^4")) {
        alps::RealObsevaluator obse_m2 = m_in["Generalized Magnetization^2"];
        alps::RealObsevaluator obse_m4 = m_in["Generalized Magnetization^4"];
        if (obse_m2.count() && obse_m4.count()) {
          alps::RealObsevaluator eval("Binder Ratio of Generalized Magnetization");
          eval = power2(obse_m2) / obse_m4;
          m.addObservable(eval);
        }
      }
      if (m_in.has("Generalized Staggered Magnetization^2") &&
          m_in.has("Generalized Staggered Magnetization^4")) {
        alps::RealObsevaluator obse_m2 =
          m_in["Generalized Staggered Magnetization^2"];
        alps::RealObsevaluator obse_m4 =
          m_in["Generalized Staggered Magnetization^4"];
        if (obse_m2.count() && obse_m4.count()) {
          alps::RealObsevaluator eval("Binder Ratio of Generalized Staggered Magnetization");
          eval = power2(obse_m2) / obse_m4;
          m.addObservable(eval);
        }
      }
    }
  };
};

} // end namespace looper

#endif // LOOPER_SUSCEPTIBILITY_H
