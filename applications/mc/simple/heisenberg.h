/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS/simple-mc
*
* Copyright (C) 1997-2015 by Synge Todo <wistaria@comp-phys.org>
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

// default parameters:
//   J = 1 (positive: ferromagnetic, negative: antiferromagnetic)
//   H = 0 (external field in Z direction)
//   T = infinity

#ifndef ALPS_SIMPLE_MC_HEISENBERG_H
#define ALPS_SIMPLE_MC_HEISENBERG_H

#include <alps/parapack/worker.h>
#include <alps/parapack/util.h> // for alps::id2string
#include <alps/random/uniform_on_sphere_n.h>
#include <boost/array.hpp>

class heisenberg_spin : public boost::array<double, 3> {
private:
  typedef boost::array<double, 3> super_type;
public:
  heisenberg_spin() { clear(); }
  explicit heisenberg_spin(std::size_t i) { clear(); }
  heisenberg_spin(double x, double y, double z) {
    super_type::operator[](0) = x;
    super_type::operator[](1) = y;
    super_type::operator[](2) = z;
  }
  void clear() { for (int i = 0; i < 3; ++i) super_type::operator[](i) = 0; }
  double operator*(heisenberg_spin const& v) const {
    double res = 0;
    for (int i = 0; i < 3; ++i) res += super_type::operator[](i) * v[i];
    return res;
  }
  heisenberg_spin& operator+=(heisenberg_spin const& v) {
    for (int i = 0; i < 3; ++i) super_type::operator[](i) += v[i];
    return *this;
  }
  alps::ODump& save(alps::ODump& dump) const {
    dump << uint32_t(3);
    dump.write_array(3, &(super_type::operator[](0)));
    return dump;
  }
  alps::IDump& load(alps::IDump& dump) {
    uint32_t dummy;
    dump >> dummy;
    dump.read_array(3, &(super_type::operator[](0)));
  }
};

inline alps::ODump& operator<<(alps::ODump& dump, heisenberg_spin const& v) {
  return v.save(dump);
}

inline alps::IDump& operator>>(alps::IDump& dump, heisenberg_spin& v) {
  return v.load(dump);
}

class heisenberg_worker : public alps::parapack::lattice_mc_worker<> {
private:
  typedef alps::parapack::lattice_mc_worker<> super_type;
  typedef heisenberg_spin spin_type;
  typedef alps::uniform_on_sphere_n<3, double, spin_type> dist_type;

public:
  heisenberg_worker(alps::Parameters const& params) : super_type(params), mcs_(params),
    spins_(num_sites()) {
    beta_ = (params.defined("T") ? 1 / evaluate("T", params) : 0);
    field_ = (params.defined("H") ? evaluate("H", params) : 0);
    // table of coupling constants
    int num_types = 0;
    bond_iterator itr, itr_end;
    for (boost::tie(itr, itr_end) = bonds(); itr != itr_end; ++itr)
      num_types = std::max(num_types, int(bond_type(*itr)) + 1);
    coupling_.clear();
    coupling_.resize(num_types, (params.defined("J")) ? evaluate("J", params) : 1);
    for (int t = 0; t < num_types; ++t) {
      std::string name = std::string("J") + boost::lexical_cast<std::string>(t);
      if (params.defined(name)) coupling_[t] = evaluate(name, params);
    }

    // random initial spins
    for (int s = 0; s < num_sites(); ++s) spins_[s] = dist_(engine());
    update_energy();

    snapshot_interval_ =
      (params.defined("SNAPSHOT_INTERVAL") ? evaluate("SNAPSHOT_INTERVAL", params) : 0);
    snapshot_basename_ = params["BASE_NAME"] + ".clone" +
      alps::id2string(static_cast<int>(params["CLONE_ID"]));
  }
  virtual ~heisenberg_worker() {}

  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::SimpleRealObservable("Temperature")
        << alps::SimpleRealObservable("Inverse Temperature")
        << alps::SimpleRealObservable("Number of Sites")
        << alps::RealObservable("Energy")
        << alps::RealObservable("Energy Density")
        << alps::RealObservable("Energy^2")
        << alps::RealObservable("Magnetization Density^2")
        << alps::RealObservable("Magnetization Density^4")
        << alps::RealObservable("Magnetization Density Z")
        << alps::RealObservable("Magnetization Density Z^2")
        << alps::RealObservable("Magnetization Density Z^4");
  }

  bool is_thermalized() const { return mcs_.is_thermalized(); }
  double progress() const { return mcs_.progress(); }

  void run(alps::ObservableSet& obs) {
    ++mcs_;
    for (int s = 0; s < num_sites(); ++s) {
      spin_type spin_new = dist_(engine());
      double diff = field_ * (spins_[s][2] - spin_new[2]);
      neighbor_bond_iterator itr, itr_end;
      for (boost::tie(itr, itr_end) = neighbor_bonds(s); itr != itr_end; ++itr) {
        int t = s ^ source(*itr) ^ target(*itr); // calculate index of neighbor spin
        diff += coupling_[bond_type(*itr)] * (spins_[s] * spins_[t] - spin_new * spins_[t]);
      }
      if (uniform_01() < std::exp(- beta_ * diff)) spins_[s] = spin_new;
    }
    spin_type mag(0, 0, 0);
    for (int s = 0; s < num_sites(); ++s) mag += spins_[s];
    update_energy();
    add_constant(obs["Temperature"], 1/beta_);
    add_constant(obs["Inverse Temperature"], beta_);
    add_constant(obs["Number of Sites"], (double)num_sites());
    obs["Energy"] << energy_;
    obs["Energy Density"] << energy_ / num_sites();
    obs["Energy^2"] << energy_ * energy_;
    obs["Magnetization Density^2"] << mag * mag / num_sites() / num_sites();
    obs["Magnetization Density^4"] << std::pow(mag * mag / num_sites() / num_sites(), 2.0);
    obs["Magnetization Density Z"] << mag[2] / num_sites();
    obs["Magnetization Density Z^2"] << mag[2] * mag[2] / num_sites() / num_sites();
    obs["Magnetization Density Z^4"] << std::pow(mag[2] * mag[2] / num_sites() / num_sites(), 2.0);

    if (snapshot_interval_ > 0 && mcs_() % snapshot_interval_ == 0) {
      std::string xdrfile = snapshot_basename_ + "." + boost::lexical_cast<std::string>(mcs_()) +
        ".snap";
      alps::OXDRFileDump dp(xdrfile);
      dp << static_cast<int>(alps::scheduler::MCDump_snapshot); // magic number
      dp << static_cast<int>(dimension())  // lattice dimension
         << static_cast<int>(3)            // spin dimension
         << static_cast<int>(num_sites());
      for (int s = 0; s < num_sites(); ++s) dp << coordinate(s) << spins_[s];
    }
  }

  // for exchange Monte Carlo
  typedef double weight_parameter_type;
  void set_beta(double beta) { beta_ = beta; }
  weight_parameter_type weight_parameter() const { return energy_; }
  static double log_weight(weight_parameter_type gw, double beta) { return - beta * gw; }

  void save(alps::ODump& dp) const { dp << mcs_ << spins_ << energy_; }
  void load(alps::IDump& dp) { dp >> mcs_ >> spins_ >> energy_; }

protected:
  void update_energy() {
    bond_iterator itr, itr_end;
    energy_ = 0;
    for (boost::tie(itr, itr_end) = bonds(); itr != itr_end; ++itr)
      energy_ -= coupling_[bond_type(*itr)] * (spins_[source(*itr)] * spins_[target(*itr)]);
    for (int s = 0; s < num_sites(); ++s)
      energy_ -= field_ * spins_[s][2];
  }

private:
  // parameteters
  double beta_;
  double field_;
  std::vector<double> coupling_;
  int snapshot_interval_;
  std::string snapshot_basename_;
  // configuration (need checkpointing)
  alps::mc_steps mcs_;
  std::vector<spin_type> spins_;
  double energy_;
  // random number distribution  
  dist_type dist_;
};

#endif // ALPS_SIMPLE_MC_HEISENBERG_H
