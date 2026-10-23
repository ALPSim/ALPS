/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2018 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Adapted from CQMP/Maxent src/maxent_grid.hpp, revision
// e727650a1f0a427e95efebc409d907b4c84011ec, on 2026-10-14.
// Uses ALPS's HDF5 archive and retains ALPS input/output conventions.
#pragma once
#include<vector>
#include <alps/params.hpp>
#include "maxent_matrix_def.hpp"

class grid{
public:
  grid(const alps::params &p);
  const std::vector<double> &t_array() const{return t_array_;}
  double operator()(int i)const{return t_array_[i];}
private:
  ///the number of (real) frequency point for this grid.
  int nfreq_;
  std::vector<double> t_array_;

  void initialize_linear_grid();
  void initialize_logarithmic_grid(double t_min);
  void initialize_quadratic_grid(double spread);
  void initialize_half_lorentzian_grid(double cut);
  void initialize_lorentzian_grid(double cut);
};
