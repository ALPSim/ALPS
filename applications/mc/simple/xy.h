/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS/simple-mc
*
* Copyright (C) 1997-2015 by Synge Todo <wistaria@comp-phys.org>
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

// default parameters:
//   J = 1 (positive: ferromagnetic, negative: antiferromagnetic)
//   H = 0 (external field in X direction)
//   T = infinity

#ifndef ALPS_SIMPLE_MC_XY_H
#define ALPS_SIMPLE_MC_XY_H

#include <alps/parapack/worker.h>
#include <alps/parapack/util.h> // for alps::id2string

class xy_worker : public alps::parapack::lattice_mc_worker<> {
private:
  typedef alps::parapack::lattice_mc_worker<> super_type;
  typedef std::complex<double> spin_type;

public:
  xy_worker(alps::Parameters const& params) : super_type(params), mcs_(params),
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
    for (int s = 0; s < num_sites(); ++s) spins_[s] = std::polar(1.0, 2 * M_PI * uniform_01());
    update_energy();

    snapshot_interval_ =
      (params.defined("SNAPSHOT_INTERVAL") ? evaluate("SNAPSHOT_INTERVAL", params) : 0);
    snapshot_basename_ = params["BASE_NAME"] + ".clone" +
      alps::id2string(static_cast<int>(params["CLONE_ID"]));
  }
  virtual ~xy_worker() {}

  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::SimpleRealObservable("Temperature")
        << alps::SimpleRealObservable("Inverse Temperature")
        << alps::SimpleRealObservable("Number of Sites")
        << alps::RealObservable("Energy")
        << alps::RealObservable("Energy Density")
        << alps::RealObservable("Energy^2")
        << alps::RealObservable("Magnetization Density^2")
        << alps::RealObservable("Magnetization Density^4")
        << alps::RealObservable("Magnetization Density X")
        << alps::RealObservable("Magnetization Density X^2")
        << alps::RealObservable("Magnetization Density X^4");
  }

  bool is_thermalized() const { return mcs_.is_thermalized(); }
  double progress() const { return mcs_.progress(); }

  void run(alps::ObservableSet& obs) {
    ++mcs_;
    for (int s = 0; s < num_sites(); ++s) {
      spin_type spin_new = std::polar(1.0, 2 * M_PI * uniform_01());
      double diff = field_ * real(spins_[s] - spin_new);
      neighbor_bond_iterator itr, itr_end;
      for (boost::tie(itr, itr_end) = neighbor_bonds(s); itr != itr_end; ++itr) {
        int t = s ^ source(*itr) ^ target(*itr); // calculate index of neighbor spin
        diff += coupling_[bond_type(*itr)] * real((spins_[s] - spin_new) * conj(spins_[t]));
      }
      if (uniform_01() < std::exp(- beta_ * diff)) spins_[s] = spin_new;
    }
    spin_type mag(0, 0);
    for (int s = 0; s < num_sites(); ++s) mag += spins_[s];
    update_energy();
    add_constant(obs["Temperature"], 1/beta_);
    add_constant(obs["Inverse Temperature"], beta_);
    add_constant(obs["Number of Sites"], (double)num_sites());
    obs["Energy"] << energy_;
    obs["Energy Density"] << energy_ / num_sites();
    obs["Energy^2"] << energy_ * energy_;
    obs["Magnetization Density^2"] << real(mag * conj(mag)) / num_sites() / num_sites();
    obs["Magnetization Density^4"] << std::pow(real(mag * conj(mag)) / num_sites() / num_sites(), 2.0);
    obs["Magnetization Density X"] << real(mag) / num_sites();
    obs["Magnetization Density X^2"] << real(mag) * real(mag) / num_sites() / num_sites();
    obs["Magnetization Density X^4"] << std::pow(real(mag) * real(mag) / num_sites() / num_sites(), 2.0);

    if (snapshot_interval_ > 0 && mcs_() % snapshot_interval_ == 0) {
      std::string xdrfile = snapshot_basename_ + "." + boost::lexical_cast<std::string>(mcs_()) +
        ".snap";
      alps::OXDRFileDump dp(xdrfile);
      dp << static_cast<int>(alps::scheduler::MCDump_snapshot); // magic number
      dp << static_cast<int>(dimension())  // lattice dimension
         << static_cast<int>(2)            // spin dimension
         << static_cast<int>(num_sites());
      std::vector<double> state(2);
      for (int s = 0; s < num_sites(); ++s) {
        state[0] = real(spins_[s]);
        state[1] = imag(spins_[s]);
        dp << coordinate(s) << state;
      }
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
      energy_ -=
        coupling_[bond_type(*itr)] * real(spins_[source(*itr)] * conj(spins_[target(*itr)]));
    for (int s = 0; s < num_sites(); ++s)
      energy_ -= field_ * real(spins_[s]);
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
};

#endif // ALPS_SIMPLE_MC_XY_H
