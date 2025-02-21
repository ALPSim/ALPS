/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2007 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_TEMPERATURE_H
#define LOOPER_TEMPERATURE_H

#include <alps/expression.h>
#include <boost/next_prior.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace looper {

class temperature {
public:
  temperature(alps::Parameters const& p) { init(p); }

  void init(alps::Parameters const& p) {
    seq_.clear();
    final_ = -1;

    if (p.defined("T")) {
      final_ = alps::evaluate("T", p);
    } else if (p.defined("TEMPERATURE")) {
      final_ = alps::evaluate("TEMPERATURE", p);
    } else if (p.defined("BETA")) {
      final_ = 1 / alps::evaluate("BETA", p);
    } else if (p.defined("INVERSE_TEMPERATURE")) {
      final_ = 1 / alps::evaluate("INVERSE_TEMPERATURE", p);
    } else {
      seq_.push_back(std::make_pair(1, final_));
      return;
    }

    int s_prev = 1;
    for (int i = 0;; ++i) {
      int td = -1;
      double ts = -1;
      std::string ns = boost::lexical_cast<std::string>(i);
      if (p.defined("T_DURATION_" + ns) && p.defined("T_START_" + ns)) {
        td = static_cast<int>(alps::evaluate("T_DURATION_" + ns, p));
        ts = alps::evaluate("T_START_" + ns, p);
      } else if (p.defined("TEMPERATURE_DURATION_" + ns) && p.defined("TEMPERATURE_START_" + ns)) {
        td = static_cast<int>(alps::evaluate("TEMPERATURE_DURATION_" + ns, p));
        ts = alps::evaluate("TEMPERATURE_START_" + ns, p);
      } else if (p.defined("BETA_DURATION_" + ns) && p.defined("BETA_START_" + ns)) {
        td = static_cast<int>(alps::evaluate("BETA_DURATION_" + ns, p));
        ts = 1 / alps::evaluate("BETA_START_" + ns, p);
      } else if (p.defined("INVERSE_TEMPERATURE_DURATION_" + ns) &&
                 p.defined("INVERSE_TEMPERATURE_START_" + ns)) {
        td = static_cast<int>(alps::evaluate("INVERSE_TEMPERATURE_DURATION_" + ns, p));
        ts = 1 / alps::evaluate("INVERSE_TEMPERATURE_START_" + ns, p);
      }
      if (td <= 0 || ts <= 0) break;
      seq_.push_back(std::make_pair(s_prev, ts));
      s_prev += td;
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
