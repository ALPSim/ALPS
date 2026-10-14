/*
 * Copyright (C) 1998-2018 ALPS Collaboration.
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

// Adapted from CQMP/Maxent src/maxent_grid.hpp, revision
// a3992534228b0de732e033ae23154969724814be, on 2026-10-14.
// Uses ALPS parameter APIs and retains ALPS input/output conventions;
// see ../README.md for provenance and licensing.
#pragma once
#include<vector>
#include<alps/ngs/params.hpp>
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
