/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef LOOPER_TEMPERATURE_H
#define LOOPER_TEMPERATURE_H

#include <alps/input/looper.hpp>
#include <boost/next_prior.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace looper {

class temperature {
public:
  temperature(alps::Parameters const& p) { init(p); }

  void init(alps::Parameters const& p) { init(alps::input::read_looper_temperature_parameters(p)); }

  void init(alps::input::looper_temperature_parameters const& p) {
    seq_.clear();
    final_ = -1;

    if (!p.temperature) {
      seq_.push_back(std::make_pair(1, final_));
      return;
    }
    final_ = *p.temperature;

    int s_prev = 1;
    for (std::size_t i = 0; i < p.annealing.size(); ++i) {
      seq_.push_back(std::make_pair(s_prev, p.annealing[i].second));
      s_prev += p.annealing[i].first;
    }
    seq_.push_back(std::make_pair(s_prev, final_));
  }

  void set_beta(double beta) {
    if (seq_.size() == 1) {
      final_ = 1 / beta;
    } else {
      boost::throw_exception(std::invalid_argument("can not reset inverse temperature"));
    }
  }

  double initial() const { return seq_.front().second; }
  double final() const { return final_; }
  int annealing_steps() const { return seq_.back().first - 1; }
  double operator()(int step = 0) const {
    if (final_ < 0)
      boost::throw_exception(std::logic_error("temperature is not initialized"));
    if (step <= 0) {
      if (seq_.size() == 1) {
        return final_;
      } else {
        boost::throw_exception(std::invalid_argument("invalid MCS"));
      }
    }
    std::vector<std::pair<int, double> >::const_iterator itr =
      std::lower_bound(seq_.begin(), seq_.end(), std::make_pair(step, 0.));
    if (itr == seq_.begin()) {
      return itr->second;
    } else if (itr != seq_.end()) {
      return boost::prior(itr)->second + (step - boost::prior(itr)->first) *
        (itr->second - boost::prior(itr)->second) / (itr->first - boost::prior(itr)->first);
    } else {
      return final_;
    }
  }

protected:
  temperature() {}

private:
  std::vector<std::pair<int, double> > seq_;
  double final_;
};

} // end namespace looper

#endif // LOOPER_TEMPERATURE_H
