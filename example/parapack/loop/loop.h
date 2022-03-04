/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef PARAPACK_EXAMPLE_LOOP_LOOP_H
#define PARAPACK_EXAMPLE_LOOP_LOOP_H

#include "union_find.h"
#include <alps/parapack/worker.h>
#include <algorithm> // for std::swap
#include <vector>

enum operator_type { diagonal, offdiagonal };

struct local_operator_t {
  local_operator_t() {}
  local_operator_t(int b, double t) : type(diagonal), bond(b), time(t) {}
  void flip() { type = (type == diagonal ? offdiagonal : diagonal); }
  operator_type type;
  unsigned int bond;
  unsigned int upper_loop, lower_loop;
  double time;
};

inline alps::ODump& operator<<(alps::ODump& od, const local_operator_t& op) {
  int t = op.type;
  od << t << op.bond << op.time;
  return od;
}

inline alps::IDump& operator>>(alps::IDump& id, local_operator_t& op) {
  int t;
  id >> t >> op.bond >> op.time;
  op.type = operator_type(t);
  return id;
}

struct cluster_t {
  cluster_t(bool t = false) : to_flip(t), size(0), mag(0), length(0) {}
  bool to_flip;
  int size;
  int mag;
  double length;
};

typedef union_find::node fragment_t;

// helper function, which returns x^2
template<typename T> inline T power2(T x) { return x * x; }

class qmc_worker : public alps::parapack::lattice_mc_worker<> {
private:
  typedef alps::parapack::lattice_mc_worker<> super_type;

public:
  qmc_worker(alps::Parameters const& params) :
    super_type(params),
    beta_(params.defined("T") ? (1.0 / evaluate("T", params)) : 1.0), volume_(volume()),
    nsites_(num_sites()), nbonds_(num_bonds()),
    mcs_(params), operators_(0), spins_(nsites_, 0),
    operators_p_(), fragments_(), current_(nsites_), clusters_(),
    r_time_(engine(), boost::exponential_distribution<>(beta_ * nbonds_ / 2)) {
    if (!is_bipartite())
      boost::throw_exception(std::invalid_argument("non-bipartite lattice"));
  }

  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::SimpleRealObservable("Volume")
        << alps::RealObservable("Energy")
        << alps::RealObservable("Energy Density")
        << alps::RealObservable("Staggered Magnetization^2")
        << alps::RealObservable("Susceptibility")
        << alps::RealObservable("Staggered Susceptibility");
  }

  bool is_thermalized() const { return mcs_.is_thermalized(); }
  double progress() const { return mcs_.progress(); }

  void run(alps::ObservableSet& obs) {
    ++mcs_;
    r_time_.distribution() = boost::exponential_distribution<>(beta_ * nbonds_ / 2);

    //
    // diagonal update and cluster construction
    //

    // initialize operator information
    std::swap(operators_, operators_p_);
    operators_.resize(0);

    // initialize cluster information (setup s cluster fragments)
    fragments_.resize(nsites_);
    std::fill(fragments_.begin(), fragments_.end(), fragment_t());
    for (int i = 0; i < nsites_; ++i) current_[i] = i;

    double t = r_time_();
    for (std::vector<local_operator_t>::iterator opi = operators_p_.begin();
         t < 1 || opi != operators_p_.end();) {

      // diagonal update
      if (opi == operators_p_.end() || t < opi->time) {
        unsigned int b = (unsigned int)(nbonds_ * random_01());
        if (spins_[source(bond(b))] != spins_[target(bond(b))]) {
          operators_.push_back(local_operator_t(b, t));
          t += r_time_();
        } else {
          t += r_time_();
          continue;
        }
      } else {
        if (opi->type == diagonal) {
          ++opi;
          continue;
        } else {
          operators_.push_back(*opi);
          ++opi;
        }
      }

      // building up clusters
      std::vector<local_operator_t>::iterator oi = operators_.end() - 1;
      unsigned int s0 = source(bond(oi->bond));
      unsigned int s1 = target(bond(oi->bond));
      oi->lower_loop = unify(fragments_, current_[s0], current_[s1]);
      oi->upper_loop = current_[s0] = current_[s1] = add(fragments_);
      if (oi->type == offdiagonal) {
        spins_[s0] ^= 1;
        spins_[s1] ^= 1;
      }
    }

    // connect bottom and top cluster fragments
    for (int s = 0; s < nsites_; ++s) unify(fragments_, s, current_[s]);

    //
    // cluster flip
    //

    // assign cluster id & determine if clusters are to be flipped
    int nc = 0;
    for (std::vector<fragment_t>::iterator ci = fragments_.begin();
         ci != fragments_.end(); ++ci) if (ci->is_root()) ci->id = nc++;
    clusters_.resize(nc);
    for (std::vector<fragment_t>::iterator ci = fragments_.begin();
         ci != fragments_.end(); ++ci) ci->id = cluster_id(fragments_, *ci);
    for (std::vector<cluster_t>::iterator pi = clusters_.begin();
         pi != clusters_.end(); ++pi) *pi = cluster_t(random_01() < 0.5);

    // 'flip' operators & do improved measurements
    for (std::vector<local_operator_t>::iterator oi = operators_.begin();
         oi != operators_.end(); ++oi) {
      int id_l = fragments_[oi->lower_loop].id;
      int id_u = fragments_[oi->upper_loop].id;
      clusters_[id_l].length += 2 * oi->time;
      clusters_[id_u].length -= 2 * oi->time;
      if (clusters_[id_l].to_flip ^ clusters_[id_u].to_flip) oi->flip();
    }

    // flip spins & do improved measurements
    for (unsigned int s = 0; s < nsites_; ++s) {
      int id = fragments_[s].id;
      clusters_[id].size += 1;
      clusters_[id].mag += 1 - 2 * spins_[s];
      clusters_[id].length += 1;
      if (clusters_[id].to_flip) spins_[s] ^= 1;
    }

    //
    // measurements
    //

    // accumurate loop length and magnetization
    double s2 = 0;
    double m2 = 0;
    double l2 = 0;
    for (std::vector<cluster_t>::const_iterator pi = clusters_.begin();
         pi != clusters_.end(); ++pi) {
      s2 += power2(pi->size);
      m2 += power2(pi->mag);
      l2 += power2(pi->length);
    }

    add_constant(obs["Volume"], volume_);
    obs["Energy"] << (0.25 * nbonds_ - operators_.size() / beta_);
    obs["Energy Density"] << (0.25 * nbonds_ - operators_.size() / beta_) / volume_;
    obs["Staggered Magnetization^2"] << 0.25 * s2;
    obs["Susceptibility"] << 0.25 * beta_ * m2 / volume_;
    obs["Staggered Susceptibility"] << 0.25 * beta_ * l2 / volume_;
  }

  // for exchange Monte Carlo
  typedef alps::hybrid_weight_parameter weight_parameter_type;
  void set_beta(double beta) { beta_ = beta; }
  weight_parameter_type weight_parameter() const {
    return weight_parameter_type(0, operators_.size());
  }
  static double log_weight(weight_parameter_type const& wp, double beta) {
    return wp.log_weight(beta);
  }

  void load(alps::IDump& dp) { dp >> mcs_ >> operators_ >> spins_; }
  void save(alps::ODump& dp) const { dp << mcs_ << operators_ << spins_; }

private:
  // parameters
  double beta_;
  double volume_;
  unsigned int nsites_;
  unsigned int nbonds_;

  // configurations (need checkpointing)
  alps::mc_steps mcs_;
  std::vector<local_operator_t> operators_;
  std::vector<int> spins_;

  // working area (no checkpointing)
  std::vector<local_operator_t> operators_p_;
  std::vector<fragment_t> fragments_;
  std::vector<unsigned int> current_;
  std::vector<cluster_t> clusters_;

  // RNG
  boost::variate_generator<engine_type&, boost::exponential_distribution<> > r_time_;
};

#endif // PARAPACK_EXAMPLE_LOOP_LOOP_H
