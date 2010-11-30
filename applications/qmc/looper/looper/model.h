/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_MODEL_H
#define LOOPER_MODEL_H

#include "graph_impl.h"
#include "lattice.h"
#include "weight.h"

#include <alps/parameter.h>
#include <alps/random/buffered_rng.h>
#include <vector>

namespace looper {

template<typename RG, typename LG, typename WH = weight_helper<> >
class spinmodel_helper {
public:
  typedef RG real_graph_t;
  typedef LG local_graph_t;
  typedef WH weight_helper_t;

  typedef lattice_helper<real_graph_t> lattice_t;
  typedef typename local_graph_t::location_t location_t;

  spinmodel_helper(alps::Parameters const& p, lattice_t& lat, bool is_path_integral = true) {
    init(p, lat, is_path_integral);
  }

  void init(alps::Parameters const& p, lattice_t& lat, bool is_path_integral = true);

  void check_parameter(bool support_longitudinal_field, bool support_negative_sign) const {
    if (has_field()) {
      if (support_longitudinal_field) {
        std::cout << "WARNING: model has a magnetic field\n";
      } else {
        boost::throw_exception(std::invalid_argument("longitudinal magnetic field not supported"));
      }
    }
    if (is_frustrated()) {
      std::cout << "WARNING: model is classically frustrated\n";
    }
    if (is_signed()) {
      if (support_negative_sign) {
        std::cout << "WARNING: model has negative signs\n";
      } else {
        boost::throw_exception(std::invalid_argument("negative signs not supported"));
      }
    }
  }

  bool is_quantal() const { return quantal_; }
  bool is_frustrated() const { return frustrated_; }
  double energy_offset() const { return site_energy_offset() + bond_energy_offset(); }
  double site_energy_offset() const { return site_offset_; }
  double bond_energy_offset() const { return bond_offset_; }

  bool has_field() const { return field_.size(); }
  std::vector<double> const& field() const { return field_; }

  bool is_signed() const { return signed_; }
  std::vector<int> const& site_sign() const { return site_sign_; }
  int site_sign(int p) const { return site_sign_[p]; }
  std::vector<int> const& bond_sign() const { return bond_sign_; }
  int bond_sign(int p) const { return bond_sign_[p]; }

  // from graph_chooser
  template<typename RNG>
  local_graph_t const& choose_graph(RNG& rng) const { return chooser_.graph(rng); }
  template<typename RNG>
  local_graph_t choose_diagonal(RNG& rng, location_t const& loc, int c) const {
    return chooser_.diagonal(rng, loc, c);
  }
  template<typename RNG>
  local_graph_t choose_diagonal(RNG& rng, location_t const& loc, int c0, int c1) const {
    return chooser_.diagonal(rng, loc, c0, c1);
  }
  template<typename RNG>
  local_graph_t choose_offdiagonal(RNG& rng, const location_t& loc, int c0, int c1) const {
    return chooser_.offdiagonal(rng, loc, c0, c1);
  }
  double graph_weight() const { return chooser_.weight(); }

private:
  bool quantal_;
  bool frustrated_;
  double site_offset_;
  double bond_offset_;
  graph_chooser<LG> chooser_;
  std::vector<double> field_;
  bool signed_;
  std::vector<int> site_sign_;
  std::vector<int> bond_sign_;
};

} // end namespace looper

#endif // LOOPER_MODEL_H
