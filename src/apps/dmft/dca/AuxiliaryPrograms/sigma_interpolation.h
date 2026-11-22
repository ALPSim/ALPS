// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

#pragma once

/// @file sigma_interpolation.h
/// @brief interpolation of a self-energy component to other Matsubara frequencies,
/// for sigma_interpolate_normal and sigma_interpolate_sc

#include <gsl/gsl_spline.h>
#include <Eigen/Dense>

#include <stdexcept>
#include <vector>

namespace sigma_interpolation {

/// The coefficients (a, b, c) of the parabola a x^2 + b x + c through the first three points.
inline Eigen::Vector3d parabola(const std::vector<double> &x, const std::vector<double> &y){
  Eigen::Matrix3d M;
  M << x[0]*x[0], x[0], 1,
       x[1]*x[1], x[1], 1,
       x[2]*x[2], x[2], 1;
  const Eigen::FullPivLU<Eigen::Matrix3d> lu(M);
  if(!lu.isInvertible()) throw std::runtime_error("sigma interpolation: the three x values are not distinct.");
  return lu.solve(Eigen::Vector3d(y[0], y[1], y[2]));
}

/// The real or imaginary part y of a self-energy component at the Matsubara
/// frequencies x (increasing), evaluated at the frequencies new_x:
/// - inside [x.front(), x.back()], with a cubic spline;
/// - below x.front(), with the parabola through the first three points;
/// - above x.back(), with the high-frequency behavior of a self-energy,
///   Sigma = Sigma_0 + Sigma_1/(i w) + ...: the real part keeps its last value,
///   the imaginary part decays as 1/w, y.back() * x.back() / w (#467).
inline std::vector<double> interpolate(const std::vector<double> &x, const std::vector<double> &y,
                                       const std::vector<double> &new_x, bool imaginary_part){
  if(x.size() < 3 || x.size() != y.size())
    throw std::invalid_argument("sigma interpolation needs a self-energy on at least 3 Matsubara frequencies.");
  std::vector<double> new_y(new_x.size());
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, x.size());
  gsl_spline_init(spline, x.data(), y.data(), x.size());
  const Eigen::Vector3d p = parabola(x, y);
  for(std::size_t i=0; i<new_x.size(); ++i){
    const double w = new_x[i];
    if(w < x.front())
      new_y[i] = p[0]*w*w + p[1]*w + p[2];
    else if(w > x.back())
      new_y[i] = imaginary_part ? y.back()*x.back()/w : y.back();
    else
      new_y[i] = gsl_spline_eval(spline, w, acc);
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  return new_y;
}

}  // namespace sigma_interpolation
