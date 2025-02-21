/*****************************************************************************
*
* ALPS Project Applications
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

#ifndef LOOPER_CUSTOM_IMPL_H
#define LOOPER_CUSTOM_IMPL_H

#include "custom.h"
#include <alps/model.h>
#include <alps/scheduler/measurement_operators.h>
#include <alps/numeric/is_nonzero.hpp>
#include <boost/foreach.hpp>
#include <boost/next_prior.hpp>

namespace looper {

template<class I>
bool build_diagonal_matrix(std::string const& op, alps::Parameters const& p,
  alps::model_helper<I> const& mh, unsigned int t, std::vector<double>& mat) {
  alps::SiteOperator term(op);
  bool valid = true;
  boost::multi_array<alps::Expression, 2> m =
    get_matrix(alps::Expression(), term, mh.model().basis().site_basis(t), p);
  std::size_t dim = m.shape()[0];
  mat.resize(dim);
  for (std::size_t i = 0; i < dim; ++i) {
    for (std::size_t j = 0; j < dim; ++j) {
      if (m[i][j].can_evaluate()) {
        if (i == j) {
          mat[i] = alps::evaluate<double>(m[i][j]);
        } else {
          if (alps::numeric::is_nonzero(alps::evaluate<double>(m[i][j]))) valid = false;
        }
      } else {
        valid = false;
      }
    }
  }
  return valid;
}

template<typename LAT>
custom_measurement_initializer<LAT>::
custom_measurement_initializer(alps::Parameters const& params) :
  alps::MeasurementOperators(params), params_(params) {}

template<typename LAT>
void custom_measurement_initializer<LAT>::
init(lattice_t const& lat,
  std::vector<s_elements_type>& average_elements,
  std::vector<s_elements_type>& local_elements,
  std::vector<p_elements_type>& correlation_elements,
  std::vector<p_elements_type>& strfactor_elements) {

  using std::max;

  if (average_expressions.size() || local_expressions.size() ||
      correlation_expressions.size() || structurefactor_expressions.size()) {

    alps::model_helper<short> mh(lat.graph_helper(), params_);

    unsigned int max_type = 0;
    std::set<unsigned int> types;
    BOOST_FOREACH(typename real_site_descriptor<lattice_t>::type const& s, sites(lat.rg())) {
      unsigned int t = get(site_type_t(), lat.rg(), s);
      types.insert(t);
      max_type = max(max_type, t);
    }

    BOOST_FOREACH(s_expression_type const& ex, average_expressions) {
      average_elements.push_back(s_elements_type());
      s_elements_type& elms = average_elements.back();
      elms.get<0>() = ex.first;
      elms.get<1>().resize(max_type + 1);
      bool valid = true;
      BOOST_FOREACH(unsigned int t, types)
        valid &= build_diagonal_matrix(ex.second, params_, mh, t, elms.get<1>()[t]);
      if (!valid) {
        std::cout << "WARNING: \"" << ex.first
                  << "\" will not be measured since it contains an off-diagonal or "
                  << "multi-site operator\n";
        average_elements.pop_back();
      }
    }
    BOOST_FOREACH(s_expression_type const& ex, local_expressions) {
      local_elements.push_back(s_elements_type());
      s_elements_type& elms = local_elements.back();
      elms.get<0>() = ex.first;
      elms.get<1>().resize(max_type + 1);
      bool valid = true;
      BOOST_FOREACH(unsigned int t, types)
        valid &= build_diagonal_matrix(ex.second, params_, mh, t, elms.get<1>()[t]);
      if (!valid) {
        std::cout << "WARNING: \"" << ex.first
                  << "\" will not be measured since it contains an off-diagonal or "
                  << "multi-site operator\n";
        local_elements.pop_back();
      }
    }
    BOOST_FOREACH(p_expression_type const& ex, correlation_expressions) {
      correlation_elements.push_back(p_elements_type());
      p_elements_type& elms = correlation_elements.back();
      elms.get<0>() = ex.first;
      elms.get<1>().resize(max_type + 1);
      elms.get<2>().resize(max_type + 1);
      bool valid = true;
      BOOST_FOREACH(unsigned int t, types) {
        valid &= build_diagonal_matrix(ex.second.first, params_, mh, t, elms.get<1>()[t]);
        valid &= build_diagonal_matrix(ex.second.second, params_, mh, t, elms.get<2>()[t]);
      }
      if (!valid) {
        std::cout << "WARNING: \"" << ex.first
                  << "\" will not be measured since it contains an off-diagonal or "
                  << "multi-site operator\n";
        correlation_elements.pop_back();
      }
    }
    BOOST_FOREACH(p_expression_type const& ex, structurefactor_expressions) {
      strfactor_elements.push_back(p_elements_type());
      p_elements_type& elms = strfactor_elements.back();
      elms.get<0>() = ex.first;
      elms.get<1>().resize(max_type + 1);
      elms.get<2>().resize(max_type + 1);
      bool valid = true;
      BOOST_FOREACH(unsigned int t, types) {
        valid &= build_diagonal_matrix(ex.second.first, params_, mh, t, elms.get<1>()[t]);
        valid &= build_diagonal_matrix(ex.second.second, params_, mh, t, elms.get<2>()[t]);
      }
      if (!valid) {
        std::cout << "WARNING: \"" << ex.first
                  << "\" will not be measured since it contains an off-diagonal or "
                  << "multi-site operator\n";
        strfactor_elements.pop_back();
      }
    }
  }
}

} // end namespace looper

#endif // LOOPER_CUSTOM_IMPL_H
