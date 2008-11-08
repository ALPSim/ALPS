/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef PARAPACK_EXAMPLE_SINGLE_ISING_H
#define PARAPACK_EXAMPLE_SINGLE_ISING_H

#include <alps/parapack/serial.h>

class single_ising_worker : public alps::parapack::lattice_mc_worker<> {
private:
  typedef alps::parapack::lattice_mc_worker<> super_type;

public:
  single_ising_worker(alps::Parameters const& params) : super_type(params), mcs_(params) {
    // temperature
    if (params.defined("T")) beta_ = 1 / evaluate("T", params);
    // coupling constant
    coupling_ = (params.defined("J")) ? evaluate("J", params) : 1.0;
    // configuration
    spins_.resize(num_sites());
    BOOST_FOREACH(site_descriptor v, sites()) spins_[v] = (uniform_01() < 0.5 ? 1 : -1);
    // energy
    energy_ = 0;
    BOOST_FOREACH(bond_descriptor b, bonds())
      energy_ -= coupling_ * spins_[source(b)] * spins_[target(b)];
  }
  virtual ~single_ising_worker() {}

  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    obs << alps::SimpleRealObservable("Temperature")
        << alps::SimpleRealObservable("Inverse Temperature")
        << alps::SimpleRealObservable("Number of Sites")
        << alps::RealObservable("Energy")
        << alps::RealObservable("Energy^2")
        << alps::RealObservable("Magnetization")
        << alps::RealObservable("Magnetization^2")
        << alps::RealObservable("Magnetization^4");
  }

  bool is_thermalized() const { return mcs_.is_thermalized(); }
  double progress() const { return mcs_.progress(); }

  void run(alps::ObservableSet& obs) {
    ++mcs_;

    BOOST_FOREACH(site_descriptor s, sites()) {
      double diff = 0;
      BOOST_FOREACH(site_descriptor v, neighbors(s))
        diff += 2 * coupling_ * spins_[s] * spins_[v];
      if (uniform_01() < 0.5 * (1 + std::tanh(-0.5 * beta_ * diff))) {
        spins_[s] = -spins_[s];
        energy_ += diff;
      }
    }

    // measurements
    double mag = 0;
    BOOST_FOREACH(site_descriptor v, sites()) mag += spins_[v];
    obs["Temperature"] << 1/beta_;
    obs["Inverse Temperature"] << beta_;
    obs["Number of Sites"] << (double)num_sites();
    obs["Energy"] << energy_;
    obs["Energy^2"] << energy_ * energy_;
    obs["Magnetization"] << mag;
    obs["Magnetization^2"] << mag * mag;
    obs["Magnetization^4"] << mag * mag * mag * mag;
  }

  // for exchange Monte Carlo
  typedef double weight_parameter_type;
  void set_beta(double beta) { beta_ = beta; }
  weight_parameter_type weight_parameter() const { return -energy_; }
  static double log_weight(weight_parameter_type gw, double beta) { return beta * gw; }

  void save(alps::ODump& dp) const { dp << mcs_ << spins_ << energy_; }
  void load(alps::IDump& dp) { dp >> mcs_ >> spins_ >> energy_; }

private:
  // parameteters
  double beta_;
  double coupling_;

  // configuration (need checkpointing)
  alps::mc_steps mcs_;
  std::vector<int> spins_;
  double energy_;
};

class ising_evaluator : public alps::parapack::simple_evaluator {
public:
  static std::string version() {
    return "ALPS/parapack test program: Metropolis algorithm " PARAPACK_VERSION;
  }
  static void print_copyright(std::ostream& out) {
    out << version()
        << "\n  Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>\n";
  }

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
  }
};

#endif // PARAPACK_EXAMPLE_SINGLE_ISING_H
