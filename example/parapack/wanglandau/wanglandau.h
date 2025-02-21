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

#include <alps/parapack/worker.h>
#include <alps/parapack/wanglandau.h>

template<typename WL_TYPE>
class wanglandau_worker : public alps::parapack::lattice_mc_worker<> {
private:
  typedef alps::parapack::lattice_mc_worker<> super_type;
  typedef alps::integer_range<int> range_t;

public:
  wanglandau_worker(alps::Parameters const& params)
    : super_type(params),
      coupling(alps::numeric::double2int(evaluate("COUPLING", params))), // integer value!
      mcs(params),
      weight(range_t(params.defined("ENERGY_WALK_RANGE")
                     ? params["ENERGY_WALK_RANGE"] : params["ENERGY_RANGE"], params),
             range_t(params.defined("ENERGY_MEASURE_RANGE")
                     ? params["ENERGY_MEASURE_RANGE"] : params["ENERGY_RANGE"], params)) {
    // configuration
    spins.resize(num_sites());
    BOOST_FOREACH(site_descriptor v, sites()) { spins[v] = (uniform_01() < 0.5 ? 1 : -1); }
    energy = 0;
    BOOST_FOREACH(bond_descriptor b, bonds()) {
      energy -= coupling * spins[source(b)] * spins[target(b)];
    }
    if (alps::wanglandau::is_measurement_phase<WL_TYPE>::value) weight.load_weight(params);
  }
  virtual ~wanglandau_worker() {}

  void init_observables(alps::Parameters const&, std::vector<alps::ObservableSet>& obs) {
    if (alps::wanglandau::is_learning_phase<WL_TYPE>::value) {
      // learning phase
      obs.resize(1);
      weight.init_observables(obs[0]);
    } else {
      // measurement phase
      obs.resize(weight.size());
      weight.init_observables(obs[0]);
      BOOST_FOREACH(alps::ObservableSet& o, obs) {
        o << alps::SimpleRealObservable("Number of Sites")
          << alps::SimpleRealObservable("Energy")
          << alps::SimpleRealObservable("Magnetization")
          << alps::SimpleRealObservable("Magnetization^2")
          << alps::SimpleRealObservable("Magnetization^4");
      }
    }
  }

  bool is_thermalized() const { return mcs.is_thermalized(); }
  double progress() const { return mcs.progress(); }

  void run(std::vector<alps::ObservableSet>& obs) {
    ++mcs;
    BOOST_FOREACH(site_descriptor s, sites()) {
      int enew = energy;
      BOOST_FOREACH(site_descriptor v, neighbors(s)) enew += 2 * coupling * spins[s] * spins[v];
      if (uniform_01() < weight.relative_weight(energy, enew)) {
        spins[s] = -spins[s];
        energy = enew;
      }
      weight.visit(obs[0], energy, mcs.factor());
    }
    weight.measure(obs[0], energy);
    if (alps::wanglandau::is_learning_phase<WL_TYPE>::value) {
      // learning phase
      if (mcs.progress() >= 1) mcs.check_flatness(obs[0], weight);
    } else {
      // measurementment phase
      if (weight.measure_range().is_included(energy)) {
        int index = weight.bin2index(energy);
        double mag = 0;
        BOOST_FOREACH(site_descriptor v, sites()) mag += spins[v];
        mag /= num_sites();
        add_constant(obs[index]["Number of Sites"], (double)num_sites());
        add_constant(obs[index]["Energy"], (double)energy);
        obs[index]["Magnetization"] << mag;
        obs[index]["Magnetization^2"] << mag * mag;
        obs[index]["Magnetization^4"] << mag * mag * mag * mag;
      }
    }
    if (mcs.progress() >= 1) obs[0] << weight;
  }

  void save(alps::ODump& dp) const { dp << mcs << spins << energy << weight; }
  void load(alps::IDump& dp) { dp >> mcs >> spins >> energy >> weight; }

private:
  // parameteters
  int coupling;

  // configuration (need checkpointing)
  alps::wanglandau_steps<WL_TYPE> mcs;
  std::vector<int> spins;
  int energy;
  alps::wanglandau_weight weight;
};

class wanglandau_reweight_worker : public alps::parapack::abstract_worker {
private:
  typedef alps::integer_range<int> range_t;

public:
  wanglandau_reweight_worker(alps::Parameters const& p)
    : params(p),
      weight(range_t(params.defined("ENERGY_MEASURE_RANGE")
                     ? params["ENERGY_MEASURE_RANGE"] : params["ENERGY_RANGE"], params)),
      t_min(evaluate("T_MIN", params)), t_max(evaluate("T_MAX", params)),
      t_step(evaluate("T_STEP", params)), nsteps(alps::numeric::double2int((t_max - t_min) / t_step)),
      has_reference(false), step(0) {
    weight.load_weight(params);
    if (params.defined("REFERENCE_BIN") && params.defined("REFERENCE_LOGG")) {
      has_reference = true;
      ref_bin = alps::numeric::double2int(evaluate("REFERENCE_BIN", params));
      ref_val = evaluate("REFERENCE_LOGG", params);
    }
  }
  virtual ~wanglandau_reweight_worker() {}

  void init_observables(alps::Parameters const&, std::vector<alps::ObservableSet>& obs) {
    obs.resize(nsteps);
    BOOST_FOREACH(alps::ObservableSet& o, obs) {
      o << alps::SimpleRealObservable("Temperature")
        << alps::SimpleRealObservable("Inverse Temperature")
        << alps::SimpleRealObservable("Number of Sites")
        << alps::SimpleRealObservable("Energy")
        << alps::SimpleRealObservable("Energy Density")
        << alps::SimpleRealObservable("Energy^2")
        << alps::SimpleRealObservable("Specific Heat")
        << alps::SimpleRealObservable("Free-Energy")
        << alps::SimpleRealObservable("Free-Energy Density")
        << alps::SimpleRealObservable("Entropy")
        << alps::SimpleRealObservable("Entropy Density")
        << alps::SimpleRealObservable("Magnetization")
        << alps::SimpleRealObservable("Magnetization^2")
        << alps::SimpleRealObservable("Magnetization^4")
        << alps::SimpleRealObservable("Binder Ratio of Magnetization");
    }
  }

  bool is_thermalized() const { return true; }
  double progress() const { return 1.0 * step / nsteps; }

  void run(std::vector<alps::ObservableSet>& obs) {
    double t = t_min + step * t_step;
    double beta = 1 / t;
    alps::IXDRFileDump dp(alps::wanglandau_weight::observable_dumpfile(params));
    while (true) {
      bool has_obs(dp);
      if (!has_obs) break;
      alps::wanglandau::observable_map_type m;
      dp >> m;
      double nsites = m["Number of Sites"][0];
      alps::exp_double part_sum = 0;
      alps::exp_double ene_sum = 0;
      alps::exp_double ene2_sum = 0;
      alps::exp_double mag_sum = 0;
      alps::exp_double mag2_sum = 0;
      alps::exp_double mag4_sum = 0;
      std::vector<double> const& hist_m = m["Measurement Histogram"];
      std::vector<double> const& ene_m = m["Energy"];
      std::vector<double> const& mag_m = m["Magnetization"];
      std::vector<double> const& mag2_m = m["Magnetization^2"];
      std::vector<double> const& mag4_m = m["Magnetization^4"];
      alps::exp_double factor = 1;
      if (has_reference) {
        if (hist_m[weight.bin2index(ref_bin)] > 0) {
          factor = alps::exp_value(ref_val - log(weight[ref_bin]))
            / alps::exp_double(hist_m[weight.bin2index(ref_bin)]);
        } else {
          std::cerr << "No count in REFERENCE_BIN. Skipped.\n";
          continue;
        }
      }
      for (int i = 0; i < hist_m.size(); ++i) {
        if (hist_m[i] > 0) {
          int bin = alps::numeric::double2int(ene_m[i]);
          alps::exp_double w
            = alps::exp_double(hist_m[i]) * factor * weight[bin]
            * alps::exp_value(- beta * ene_m[i]);
          part_sum += w;
          ene_sum += ene_m[i] * w;
          ene2_sum += ene_m[i] * ene_m[i] * w;
          mag_sum += mag_m[i] * w;
          mag2_sum += mag2_m[i] * w;
          mag4_sum += mag4_m[i] * w;
        }
      }
      double ene = static_cast<double>(ene_sum / part_sum);
      double ene2 = static_cast<double>(ene2_sum / part_sum);
      double mag2 = static_cast<double>(mag2_sum / part_sum);
      double mag4 = static_cast<double>(mag4_sum / part_sum);
      obs[step]["Temperature"] << t;
      obs[step]["Inverse Temperature"] << beta;
      obs[step]["Number of Sites"] << nsites;
      obs[step]["Energy"] << ene;
      obs[step]["Energy Density"] << ene / nsites;
      obs[step]["Energy^2"] << ene2;
      obs[step]["Specific Heat"] << beta * beta * (ene2 - ene * ene) / nsites;
      obs[step]["Free-Energy"] << -t * log(part_sum);
      obs[step]["Free-Energy Density"] << -t * log(part_sum) / nsites;
      obs[step]["Entropy"] << beta * ene + log(part_sum);
      obs[step]["Entropy Density"] << (beta * ene + log(part_sum)) / nsites;
      obs[step]["Magnetization"] << static_cast<double>(mag_sum / part_sum);
      obs[step]["Magnetization^2"] << mag2;
      obs[step]["Magnetization^4"] << mag4;
      obs[step]["Binder Ratio of Magnetization"] << mag2 * mag2 / mag4;
    }
    ++step;
  }

  void save(alps::ODump& dp) const { dp << step; }
  void load(alps::IDump& dp) { dp >> step; }

private:
  // parameters
  alps::Parameters params;
  alps::wanglandau_weight weight;
  double t_min;
  double t_max;
  double t_step;
  int nsteps;
  bool has_reference;
  int ref_bin;
  double ref_val;

  // configuration
  int step;
};

class wanglandau_reweight_evaluator : public alps::parapack::simple_evaluator {
public:
  wanglandau_reweight_evaluator(alps::Parameters const&) {}
  void evaluate(std::vector<alps::ObservableSet>& obs) const {
    BOOST_FOREACH(alps::ObservableSet& o, obs) {
      if (o.has("Inverse Temperature") && o.has("Number of Sites") &&
          o.has("Energy") && o.has("Energy^2")) {
        alps::RealObsevaluator beta = o["Inverse Temperature"];
        alps::RealObsevaluator n = o["Number of Sites"];
        alps::RealObsevaluator ene = o["Energy"];
        alps::RealObsevaluator ene2 = o["Energy^2"];
        if (beta.count() && n.count() && ene.count() >= 8 && ene2.count()) {
          alps::RealObsevaluator c("Specific Heat");
          c = beta.mean() * beta.mean() * (ene2 - ene * ene) / n.mean();
          o.addObservable(c);
        }
      }
      if (o.has("Magnetization^2") && o.has("Magnetization^4")) {
        alps::RealObsevaluator m2 = o["Magnetization^2"];
        alps::RealObsevaluator m4 = o["Magnetization^4"];
        if (m2.count() >= 8 && m4.count()) {
          alps::RealObsevaluator binder("Binder Ratio of Magnetization");
          binder = m2 * m2 / m4;
          o.addObservable(binder);
        }
      }
    }
  }
};
