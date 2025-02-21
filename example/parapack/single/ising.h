/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
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

#ifndef PARAPACK_EXAMPLE_SINGLE_ISING_H
#define PARAPACK_EXAMPLE_SINGLE_ISING_H

#include <alps/parapack/worker.h>
#include <boost/graph/sequential_vertex_coloring.hpp>

class single_ising_worker : public alps::parapack::lattice_mc_worker<> {
private:
  typedef alps::parapack::lattice_mc_worker<> super_type;

public:
  single_ising_worker(alps::Parameters const& params) : super_type(params), mcs_(params) {
    // temperature
    if (params.defined("T")) beta_ = 1 / evaluate("T", params);
    // coupling constant
    coupling_ = (params.defined("J")) ? evaluate("J", params) : 1.0;
    // lattice
    std::vector<std::size_t> color(num_sites());
    typedef boost::property_map<graph_type, alps::site_index_t>::const_type vertex_index_map;
    int nc = boost::sequential_vertex_coloring(graph(),
      boost::iterator_property_map<std::size_t*, vertex_index_map>(&color.front(),
          get(boost::vertex_index, graph())));
    sublat_.clear();
    sublat_.resize(nc);
    for (std::size_t s = 0; s < num_sites(); ++s) sublat_[color[s]].push_back(s);
    // configuration
    spins_.resize(num_sites());
    #ifdef ALPS_ENABLE_OPENMP_WORKER
    #pragma omp parallel
    #endif 
    {
      #ifdef ALPS_ENABLE_OPENMP_WORKER
      int r = alps::thread_id();
      #pragma omp for
      #else
      int r = 0;
      #endif 
      for (int s = 0; s < num_sites(); ++s) spins_[s] = (uniform_01(r) < 0.5 ? 1 : -1);
    }
    double ene = 0;
    #ifdef ALPS_ENABLE_OPENMP_WORKER
    #pragma omp parallel for reduction(+: ene)
    #endif 
    for (int b = 0; b < num_bonds(); ++b) {
      bond_descriptor bd = bond(b);
      ene -= coupling_ * spins_[source(bd)] * spins_[target(bd)];
    }
    energy_ = ene;
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

    for (std::size_t t = 0; t < sublat_.size(); ++t) {
      #ifdef ALPS_ENABLE_OPENMP_WORKER
      #pragma omp parallel
      #endif 
      {
        #ifdef ALPS_ENABLE_OPENMP_WORKER
        int r = alps::thread_id();
        #pragma omp for
        #else
        int r = 0;
        #endif
        for (int k = 0; k < sublat_[t].size(); ++k) {
          int s = sublat_[t][k];
          double diff = 0;
          neighbor_iterator itr, itr_end;
          for (boost::tie(itr, itr_end) = neighbors(s); itr != itr_end; ++itr)
            diff += 2 * coupling_ * spins_[s] * spins_[*itr];
          if (uniform_01(r) < 0.5 * (1 + std::tanh(-0.5 * beta_ * diff))) {
            spins_[s] = -spins_[s];
          }
        } // end omp for
      } // end omp parallel
    }

    // measurements
    double mag = 0;
    #ifdef ALPS_ENABLE_OPENMP_WORKER
    #pragma omp parallel for reduction(+: mag)
    #endif 
    for (int s = 0; s < num_sites(); ++s) mag += spins_[s];
    double ene = 0;
    #ifdef ALPS_ENABLE_OPENMP_WORKER
    #pragma omp parallel for reduction(+: ene)
    #endif 
    for (int b = 0; b < num_bonds(); ++b) {
      bond_descriptor bd = bond(b);
      ene -= coupling_ * spins_[source(bd)] * spins_[target(bd)];
    }
    energy_ = ene;
    add_constant(obs["Temperature"], 1/beta_);
    add_constant(obs["Inverse Temperature"], beta_);
    add_constant(obs["Number of Sites"], (double)num_sites());
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
  std::vector<std::vector<int> > sublat_;

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

#endif // PARAPACK_EXAMPLE_SINGLE_ISING_H
