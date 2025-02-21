/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef PARAPACK_EXAMPLE_ISING_ISING_H
#define PARAPACK_EXAMPLE_ISING_ISING_H

#include <alps/parapack/worker.h>

class ising_worker : public alps::parapack::lattice_mc_worker<> {
private:
  typedef alps::parapack::lattice_mc_worker<> super_type;

public:
  ising_worker(alps::Parameters const& params) : super_type(params), mcs_(params),
    spins_(num_sites()) {
    beta_ = (params.defined("T") ? 1 / evaluate("T", params) : 0);
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
    // std::cerr << "coupling constants = " << alps::write_vector(coupling_) << std::endl;
    // random initial spins
    for (int s = 0; s < num_sites(); ++s) spins_[s] = (uniform_01() < 0.5 ? 1 : -1);
    update_energy();
  }
  virtual ~ising_worker() {}

  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::SimpleRealObservable("Temperature")
        << alps::SimpleRealObservable("Inverse Temperature")
        << alps::SimpleRealObservable("Number of Sites")
        << alps::RealObservable("Energy")
        << alps::RealObservable("Energy Density")
        << alps::RealObservable("Energy^2")
        << alps::RealObservable("Magnetization")
        << alps::RealObservable("Magnetization^2")
        << alps::RealObservable("Magnetization^4");
  }

  bool is_thermalized() const { return mcs_.is_thermalized(); }
  double progress() const { return mcs_.progress(); }

  void run(alps::ObservableSet& obs) {
    ++mcs_;
    for (int s = 0; s < num_sites(); ++s) {
      double diff = 0;
      neighbor_bond_iterator itr, itr_end;
      for (boost::tie(itr, itr_end) = neighbor_bonds(s); itr != itr_end; ++itr) {
        int t = s ^ source(*itr) ^ target(*itr); // calculate index of neighbor spin
        diff += 2 * coupling_[bond_type(*itr)] * spins_[s] * spins_[t];
      }
      if (uniform_01() < 0.5 * (1 + std::tanh(-0.5 * beta_ * diff))) spins_[s] = -spins_[s];
    }
    double mag = 0;
    for (int s = 0; s < num_sites(); ++s) mag += spins_[s];
    update_energy();
    add_constant(obs["Temperature"], 1/beta_);
    add_constant(obs["Inverse Temperature"], beta_);
    add_constant(obs["Number of Sites"], (double)num_sites());
    obs["Energy"] << energy_;
    obs["Energy Density"] << energy_ / num_sites();
    obs["Energy^2"] << energy_ * energy_;
    obs["Magnetization"] << mag / num_sites();
    obs["Magnetization^2"] << mag * mag / num_sites() / num_sites();
    obs["Magnetization^4"] << std::pow(mag * mag / num_sites() / num_sites(), 2.0);
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
  }

private:
  // parameteters
  double beta_;
  std::vector<double> coupling_;
  // configuration (need checkpointing)
  alps::mc_steps mcs_;
  std::vector<int> spins_;
  double energy_;
};

class ising_evaluator : public alps::parapack::simple_evaluator {
public:
  ising_evaluator(alps::Parameters const&) {}
  virtual ~ising_evaluator() {}

  void evaluate(alps::ObservableSet& obs) const {
    if (obs.has("Inverse Temperature") && obs.has("Number of Sites") &&
        obs.has("Energy") && obs.has("Energy^2")) {
      alps::RealObsevaluator beta = obs["Inverse Temperature"];
      alps::RealObsevaluator n = obs["Number of Sites"];
      alps::RealObsevaluator ene = obs["Energy"];
      alps::RealObsevaluator ene2 = obs["Energy^2"];
      if (beta.count() && n.count() && ene.count() && ene2.count()) {
        alps::RealObsevaluator c("Specific Heat");
        c = beta.mean() * beta.mean() * (ene2 - ene * ene) / n.mean();
        obs.addObservable(c);
      }
    }
    if (obs.has("Magnetization^2") && obs.has("Magnetization^4")) {
      alps::RealObsevaluator m2 = obs["Magnetization^2"];
      alps::RealObsevaluator m4 = obs["Magnetization^4"];
      if (m2.count() && m4.count()) {
        alps::RealObsevaluator binder("Binder Ratio of Magnetization");
        binder = m2 * m2 / m4;
        obs.addObservable(binder);
      }
    }
    if (obs.has("Inverse Temperature")) obs.erase("Inverse Temperature");
    if (obs.has("Temperature")) obs.erase("Temperature");
  }
};

#endif // PARAPACK_EXAMPLE_ISING_ISING_H
