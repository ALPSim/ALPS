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

#ifndef LOOPER_MODEL_IMPL_H
#define LOOPER_MODEL_IMPL_H

#include "model.h"
#include "model_parameter.h"
#include "weight_impl.h"

#include <alps/model.h>
#include <boost/foreach.hpp>

namespace looper {

template<typename RG, typename LG, typename WH>
void spinmodel_helper<RG, LG, WH>::init(alps::Parameters const& p, lattice_helper<RG>& lat,
  bool is_path_integral) {
  typedef lattice_helper<RG> lattice_t;
  alps::model_helper<short> mh(lat.graph_helper(), p);
  model_parameter mp(p, lat.graph_helper(), mh);
  quantal_ = mp.is_quantal();
  frustrated_ = mp.is_frustrated();

  // virtual graph
  lat.generate_virtual_graph(mp, mp.has_d_term());

  // weight table
  alps::Parameters params(p);
  if (frustrated_ && !params.defined("FORCE_SCATTER")) params["FORCE_SCATTER"] = 0.1;
  weight_helper_t wt(mp, lat, params);
  chooser_.init(wt, is_path_integral);

  // signs
  site_sign_.clear();
  bond_sign_.clear();
  signed_ = mp.is_signed();
  if (signed_) {
    site_sign_.resize(num_sites(lat.vg()));
    BOOST_FOREACH(typename weight_helper_t::site_weight_t const& s, wt.site_weights())
      site_sign_[s.first] = (s.second.sign < 0) ? 1 : 0;
    bond_sign_.resize(num_bonds(lat.vg()));
    BOOST_FOREACH(typename weight_helper_t::bond_weight_t const& b, wt.bond_weights())
      bond_sign_[b.first] = (b.second.sign < 0) ? 1 : 0;
  }

  // external longitudinal field
  field_.clear();
  if (mp.has_field()) {
    field_.resize(num_sites(lat.vg()));
    BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type rs, sites(lat.rg()))
      BOOST_FOREACH(typename virtual_site_descriptor<lattice_t>::type vs, sites(lat, rs))
        field_[vs] = mp.site(rs, lat.rg()).hz;
  }

  // energy offset
  site_offset_ = mp.site_energy_offset() + wt.site_energy_offset();
  bond_offset_ = mp.bond_energy_offset() + wt.bond_energy_offset();
}

} // end namespace looper

#endif // LOOPER_MODEL_IMPL_H
