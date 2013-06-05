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

#ifndef PARAPACK_EXAMPLE_MULTIPLE_ISING_H
#define PARAPACK_EXAMPLE_MULTIPLE_ISING_H

#include <alps/parapack/worker.h>
#include <functional>

namespace mpi = boost::mpi;

// a vector with 'tabs' on both sides: elements at i = -1 as well as
// i = size() are also accessible
template<typename T>
class tabbed_vector {
public:
  typedef T value_type;
  typedef typename std::vector<value_type>::size_type size_type;

  explicit tabbed_vector(size_type n = 0) : vector_(n + 2) {}
  explicit tabbed_vector(size_type n, value_type x) : vector_(n + 2, x) {}
  void resize(size_type n, value_type x = value_type()) { vector_.resize(n + 2, x); }

  size_type size() const { return vector_.size() - 2; }
  value_type const& operator[](int i) const { return vector_[i+1]; }
  value_type& operator[](int i) { return vector_[i+1]; }

  void save(alps::ODump& dp) const { dp << vector_; }
  void load(alps::IDump& dp) { dp >> vector_; }

private:
  std::vector<value_type> vector_;
};

template<typename T>
alps::ODump& operator<<(alps::ODump& dp, tabbed_vector<T> const& v) { v.save(dp); return dp; }

template<typename T>
alps::IDump& operator>>(alps::IDump& dp, tabbed_vector<T>& v) { v.load(dp); return dp; }


class parallel_ising_worker : public alps::parapack::mc_worker {
private:
  typedef alps::parapack::mc_worker super_type;

public:
  parallel_ising_worker(mpi::communicator const& comm, alps::Parameters const& params)
    : super_type(params), comm_(comm), mcs_(params) {
    // temperature
    if (params.defined("T")) beta_ = 1 / evaluate("T", params);
    coupling_ = (params.defined("J")) ? evaluate("J", params) : 1.0;
    // system size and local system size
    length_ = static_cast<int>(evaluate("L", params));
    if (comm_.rank() == 0)
      loclen_ = length_ - (comm_.size() - 1) * (length_ / comm_.size());
    else
      loclen_ = length_ / comm_.size();
    if (loclen_ < 2)
      boost::throw_exception(std::runtime_error("too small system size"));

    // configuration
    spins_.resize(loclen_);
    for (int i = 0; i < loclen_; ++i) spins_[i] = (uniform_01() < 0.5 ? 1 : 0);
    copy2right();
    copy2left();
  }
  virtual ~parallel_ising_worker() {}

  void init_observables(alps::Parameters const&, alps::ObservableSet& obs) {
    if (comm_.rank() == 0)
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

    for (int i = 0; i < loclen_; ++i) {
      double diff = coupling_ * (4 * (spins_[i-1] ^ spins_[i] + spins_[i] ^ spins_[i+1]) - 4);
      if (uniform_01() < 0.5 * (1 + std::tanh(-0.5 * beta_ * diff))) spins_[i] ^= 1;
      if (i == 0) copy2left();
      if (i == loclen_ - 1) copy2right();
    }

    // measurements
    energy_ = 0;
    double mag = 0;
    for (int i = 0; i < loclen_; ++i) {
      energy_ -= coupling_ * (2 * (spins_[i] ^ spins_[i+1]) - 1);
      mag += (2 * spins_[i] - 1);
    }
    if (comm_.rank() == 0) {
      double energy_out, mag_out;
      reduce(comm_, energy_, energy_out, std::plus<double>(), 0);
      reduce(comm_, mag, mag_out, std::plus<double>(), 0);
      energy_ = energy_out;
      mag = mag_out;
    } else {
      reduce(comm_, energy_, std::plus<double>(), 0);
      reduce(comm_, mag, std::plus<double>(), 0);
    }

    if (comm_.rank() == 0) {
      add_constant(obs["Temperature"], 1/beta_);
      add_constant(obs["Inverse Temperature"], beta_);
      add_constant(obs["Number of Sites"], (double)length_);
      obs["Energy"] << energy_;
      obs["Energy^2"] << energy_ * energy_;
      obs["Magnetization"] << mag;
      obs["Magnetization^2"] << mag * mag;
      obs["Magnetization^4"] << mag * mag * mag * mag;
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
  void copy2right() {
    if (comm_.size() == 1) {
      spins_[-1] = spins_[loclen_-1];
    } else {
      comm_.send((comm_.rank() + 1) % comm_.size(), 0, spins_[loclen_-1]);
      comm_.recv((comm_.rank() + comm_.size() - 1) % comm_.size(), 0, spins_[-1]);
    }
  }

  void copy2left() {
    if (comm_.size() == 1) {
      spins_[loclen_] = spins_[0];
    } else {
      comm_.send((comm_.rank() + comm_.size() - 1) % comm_.size(), 0, spins_[0]);
      comm_.recv((comm_.rank() + 1) % comm_.size(), 0, spins_[loclen_]);
    }
  }

private:
  mpi::communicator comm_;

  // parameteters
  double beta_;
  double coupling_;
  int length_;
  int loclen_;

  // configuration (need checkpointing)
  alps::mc_steps mcs_;
  tabbed_vector<int> spins_;
  double energy_;
};

#endif // PARAPACK_EXAMPLE_MULTIPLE_ISING_H
