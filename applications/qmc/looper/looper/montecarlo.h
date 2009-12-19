/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_MONTECARLO_H
#define LOOPER_MONTECARLO_H

#include <alps/expression.h>
#include <alps/osiris.h>
#include <alps/parapack/integer_range.h>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace looper {

class mc_steps {
public:
  typedef alps::integer_range<unsigned int> range_type;

  mc_steps() : mcs_(0), sweep_("[0:]"), therm_(0) {}
  mc_steps(alps::Parameters const& p) : mcs_(0),
    sweep_(p.value_or_default("SWEEPS", "[65536:]"), p),
    therm_(p.defined("THERMALIZATION") ? static_cast<int>(alps::evaluate("THERMALIZATION", p)) :
      (sweep_.min BOOST_PREVENT_MACRO_SUBSTITUTION () >> 3)) {
  }

  void set_thermalization(int c) { therm_ = c; }
  void set_sweeps(unsigned int c) { sweep_ = range_type(c); }
  void set_sweeps(range_type const& c) { sweep_ = c; }

  mc_steps& operator++() { ++mcs_; return *this; }
  mc_steps operator++(int) { mc_steps tmp = *this; this->operator++(); return tmp; }

  unsigned int operator()() const { return mcs_; }
  bool can_work() const {
    return mcs_ < therm_ || mcs_ - therm_ < sweep_.max BOOST_PREVENT_MACRO_SUBSTITUTION ();
  }
  bool is_thermalized() const { return mcs_ >= therm_; }
  double progress() const {
    return static_cast<double>(mcs_) / (therm_ + sweep_.min BOOST_PREVENT_MACRO_SUBSTITUTION ());
  }

  int thermalization() const { return therm_; }
  range_type sweeps() const { return sweep_; }

  void save(alps::ODump& dp) const { dp << mcs_ << sweep_ << therm_; }
  void load(alps::IDump& dp) { dp >> mcs_ >> sweep_ >> therm_; }

private:
  unsigned int mcs_;
  range_type sweep_;
  unsigned int therm_;
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dp, looper::mc_steps const& mcs) {
  mcs.save(dp);
  return dp;
}

inline alps::IDump& operator>>(alps::IDump& dp, looper::mc_steps& mcs) {
  mcs.load(dp);
  return dp;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

namespace looper {

class wl_steps {
public:
  wl_steps() : mcs_(0), stage_(0) {}
  template<class T>
  wl_steps(alps::Parameters const& p, alps::integer_range<T> const& exp_range)
    : mcs_(0), stage_(0),
      multicanonical_(!p.defined("DISABLE_MULTICANONICAL_MEASUREMENT")),
      zhou_bhatt_(!p.defined("DISABLE_ZHOU_BHATT")),
      block_(zhou_bhatt_ ? exp_range.size() : static_cast<int>(
        p.value_or_default("BLOCK_SWEEPS", 65536))),
      sweep_(p.value_or_default("SWEEPS", 65536)),
      iteration_(p.value_or_default("WANG_LANDAU_ITERATIONS", 16)) {
  }

  wl_steps& operator++() { ++mcs_; return *this; }
  wl_steps operator++(int) { wl_steps tmp = *this; this->operator++(); return tmp; }
  void reset_stage() { mcs_ = 0; }
  void next_stage() {
    ++stage_;
    mcs_ = 0;
    if (use_zhou_bhatt()) block_ = static_cast<int>(1.41 * block_);
  }

  int operator()() const { return mcs_; }
  int stage() const { return stage_; }
  bool can_work() const {
    if (perform_multicanonical_measurement())
      return !doing_multicanonical() || mcs_ < sweep_;
    else
      return stage_ < iteration_;
  }
  bool doing_multicanonical() const { return stage_ == iteration_; }
  double progress() const {
    if (perform_multicanonical_measurement())
      return (doing_multicanonical() && mcs_ >= sweep_) ? 1 : 0;
    else
      return (stage_ + 1 >= iteration_ && mcs_ >= block_) ? 1 : 0;
  }

  bool perform_multicanonical_measurement() const { return multicanonical_; }
  bool use_zhou_bhatt() const { return zhou_bhatt_; }
  int block() const { return block_; }
  int sweeps() const { return sweep_; }
  int iterations() const { return iteration_; }

  void save(alps::ODump& dp) const {
    dp << mcs_ << stage_ << multicanonical_ << zhou_bhatt_ << block_ << sweep_ << iteration_;
  }
  void load(alps::IDump& dp) {
    dp >> mcs_ >> stage_ >> multicanonical_ >> zhou_bhatt_ >> block_ >> sweep_ >> iteration_;
  }

private:
  int mcs_;
  int stage_;
  bool multicanonical_;
  bool zhou_bhatt_;
  int block_;
  int sweep_;
  int iteration_;
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dp, looper::wl_steps const& mcs) {
  mcs.save(dp);
  return dp;
}

inline alps::IDump& operator>>(alps::IDump& dp, looper::wl_steps& mcs) {
  mcs.load(dp);
  return dp;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_MONTECARLO_H
