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

#ifndef LOOPER_STIFFNESS_H
#define LOOPER_STIFFNESS_H

#include "measurement.h"

// workaround for SuSE 11.4, which defines macro TIME in pyconfig.h
#ifdef TIME
# undef TIME
#endif

namespace looper {

//
// stiffness measurement
//

template<unsigned int MAX_DIM>
struct stiffness
{
  template<typename MC, typename LAT, typename TIME>
  struct estimator
  {
    typedef MC   mc_type;
    typedef LAT  lattice_t;
    typedef TIME time_t;
    typedef estimator<mc_type, lattice_t, time_t> estimator_t;
    typedef typename alps::property_map<real_bond_t,
              const typename lattice_t::virtual_graph_type,
              typename real_bond_descriptor<lattice_t>::type>::type
              real_bond_map_t;
    typedef typename alps::property_map<alps::bond_vector_relative_t,
              const typename lattice_t::real_graph_type,
              alps::coordinate_type>::type bond_vector_relative_map_t;

    bool improved;
    real_bond_map_t real_bond;
    bond_vector_relative_map_t bond_vector_relative;
    unsigned int dim;

    void initialize(alps::Parameters const& /* params */, lattice_t const& lat,
      bool /* is_signed */, bool use_improved_estimator) {
      improved = use_improved_estimator;
      real_bond = alps::get_or_default(real_bond_t(), lat.vg(),
                                       typename real_bond_descriptor<lattice_t>::type());
      bond_vector_relative = alps::get_or_default(bond_vector_relative_t(), lat.rg(),
                                                  coordinate_type());
      dim = get_property(lat.rg(), dimension_t());
      if (improved && dim > MAX_DIM) {
        std::cout << "WARNING: Spatial dimension (=" << dim << ") is too large.  "
                  << "Stiffness will be measured only for the first "
                  << MAX_DIM << " dimensions\n";
        dim = MAX_DIM;
      }
    }
    template<typename M>
    void init_observables(M& m, bool is_signed) {
      if (dim > 0) add_scalar_obs(m, "Stiffness", is_signed);
    }

    // improved estimator

    struct estimate {
      double winding[MAX_DIM];
      estimate()
      {
          std::fill(winding, winding+MAX_DIM, 0.0);
      }
      void init() {
        for (int i = 0; i < MAX_DIM; ++i) winding[i] = 0;
      }
      estimate& operator+=(estimate const& rhs) {
        for (int i = 0; i < MAX_DIM; ++i) winding[i] += rhs.winding[i];
        return *this;
      }
      void begin_s(estimator_t const&, lattice_t const&, double, int, int) const {}
      void begin_bs(estimator_t const& emt, lattice_t const& lat, double, int b, int, int c) {
        alps::coordinate_type const& vr =
          emt.bond_vector_relative[emt.real_bond[bond(lat.vg(), b)]];
        for (int i = 0; i < emt.dim; ++i) winding[i] -= (1-2*c) * vr[i];
      }
      void begin_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
      void end_s(estimator_t const&, lattice_t const&, double, int, int) const {}
      void end_bs(estimator_t const& emt, lattice_t const& lat, double, int b, int, int c) {
        alps::coordinate_type const& vr =
          emt.bond_vector_relative[emt.real_bond[bond(lat.vg(), b)]];
        for (int i = 0; i < emt.dim; ++i) winding[i] += (1-2*c) * vr[i];
      }
      void end_bt(estimator_t const&, lattice_t const&, double, int, int, int) {}
      void start_bottom(estimator_t const&, lattice_t const&, double, int, int) const {}
      void start(estimator_t const&, lattice_t const&, double, int, int) const {}
      void stop(estimator_t const&, lattice_t const&, double, int, int) const {}
      void stop_top(estimator_t const&, lattice_t const&, double, int, int) const {}
    };
    void init_estimate(estimate& est) const { est.init(); }

    struct collector {
      unsigned int dim;
      double w2;
      void init(unsigned int d) { dim = d; w2 = 0; }
      collector& operator+=(collector const& coll) {
        for (int i = 0; i < dim; ++i) w2 += coll.w2;
        return *this;
      }
      collector& operator+=(estimate const& est) {
        for (int i = 0; i < dim; ++i) w2 += power2(0.5 * est.winding[i]);
        return *this;
      }
      template<typename M>
      void commit(M& m, lattice_t const&, double beta, int, double sign) const {
        if (dim > 0) m["Stiffness"] << sign * w2 / (beta * dim);
      }
    };
    void init_collector(collector& coll) const { coll.init(dim); }

    template<typename M, typename OP, typename FRAGMENT>
    void improved_measurement(M& m,
                              lattice_t const& lat,
                              double beta, double sign,
                              std::vector<int> const& /* spins */,
                              std::vector<OP> const& operators,
                              std::vector<int> const& /* spins_c */,
                              std::vector<FRAGMENT> const& /* fragments */,
                              collector const& coll)
    { coll.commit(m, lat, beta, operators.size(), sign); }

    // normal estimator

    template<typename M, typename OP>
    void normal_measurement(M& m, lattice_t const& lat,
                            double beta, double sign,
                            std::vector<int> const& spins,
                            std::vector<OP> const& operators,
                            std::vector<int>& spins_c)
    {
      if (improved || dim == 0) return;

      std::valarray<double> winding(0., dim);
      std::copy(spins.begin(), spins.end(), spins_c.begin());
      for (typename std::vector<OP>::const_iterator oi = operators.begin();
           oi != operators.end(); ++oi) {
        if (oi->is_offdiagonal()) {
          if (oi->is_bond()) {
            double s = 1-2*spins_c[source(oi->pos(), lat.vg())];
            alps::coordinate_type const& vr =
              bond_vector_relative[real_bond[bond(lat.vg(), oi->pos())]];
            for (int i = 0; i < dim; ++i) winding[i] += s * vr[i];
            spins_c[source(oi->pos(), lat.vg())] ^= 1;
            spins_c[target(oi->pos(), lat.vg())] ^= 1;
          } else {
            spins_c[oi->pos()] ^= 1;
          }
        }
      }

      double w2 = 0;
      for (int i = 0; i < dim; ++i) w2 += power2(winding[i]);
      m["Stiffness"] << sign * w2 / (beta * dim);
    }
  };
};

} // end namespace looper

#endif // LOOPER_STIFFNESS_H
