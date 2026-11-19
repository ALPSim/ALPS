/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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

  void init_observables(alps::Parameters const&, alps::alea::observable_set& obs) {
    obs << alps::alea::observable("Temperature").without_bins()
        << alps::alea::observable("Inverse Temperature").without_bins()
        << alps::alea::observable("Number of Sites").without_bins()
        << alps::alea::observable("Energy")
        << alps::alea::observable("Energy Density")
        << alps::alea::observable("Energy^2")
        << alps::alea::observable("Magnetization")
        << alps::alea::observable("Magnetization^2")
        << alps::alea::observable("Magnetization^4");
  }

  bool is_thermalized() const { return mcs_.is_thermalized(); }
  double progress() const { return mcs_.progress(); }

  void run(alps::alea::observable_set& obs) {
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
    alps::add_constant(obs["Temperature"], 1/beta_);
    alps::add_constant(obs["Inverse Temperature"], beta_);
    alps::add_constant(obs["Number of Sites"], (double)num_sites());
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

  void evaluate(alps::alea::result_set& obs) const {
    if (obs.has("Inverse Temperature") && obs.has("Number of Sites") &&
        obs.has("Energy") && obs.has("Energy^2")) {
      const alps::alea::observable_result beta = obs["Inverse Temperature"];
      const alps::alea::observable_result n = obs["Number of Sites"];
      const alps::alea::observable_result ene = obs["Energy"];
      const alps::alea::observable_result ene2 = obs["Energy^2"];
      if (beta.count() && n.count() && ene.count() && ene2.count()) {
        obs.insert("Specific Heat", beta.mean() * beta.mean() * (ene2 - ene * ene) / n.mean());
      }
    }
    if (obs.has("Magnetization^2") && obs.has("Magnetization^4")) {
      const alps::alea::observable_result m2 = obs["Magnetization^2"];
      const alps::alea::observable_result m4 = obs["Magnetization^4"];
      if (m2.count() && m4.count()) {
        obs.insert("Binder Ratio of Magnetization", m2 * m2 / m4);
      }
    }
    if (obs.has("Inverse Temperature")) obs.erase("Inverse Temperature");
    if (obs.has("Temperature")) obs.erase("Temperature");
  }
};

#endif // PARAPACK_EXAMPLE_ISING_ISING_H
