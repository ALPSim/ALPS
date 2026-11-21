// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT

#ifndef ALPS_DMFT_RETARDED_INTERACTION_H
#define ALPS_DMFT_RETARDED_INTERACTION_H

/// @file retarded_interaction.h
/// @brief reads the retarded interaction K(tau) and K'(tau) (RET_INT_K)
///
/// One reader for the hybridization and ct-hyb-segment solvers and for dmft,
/// which needs K'(0) for the half-filling shift of MU. The solvers each read
/// the file themselves before, without checking its length: a shorter HDF5
/// array was read past its end, and a longer file, or one with more lines,
/// was cut to N+1 values and taken to span [0, beta].

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/// K(tau) and K'(tau) on the n_tau+1 points tau_i = i beta / n_tau.
struct retarded_interaction {
  std::vector<double> K, Kp;
};

namespace retarded_interaction_detail {
inline std::runtime_error error(std::string const& file, std::string const& what) {
  return std::runtime_error("RET_INT_K " + file + ": " + what);
}
}

/// Reads K and K' from the text file `file` (n_tau+1 lines "tau-or-index K K'",
/// the first column is not used) or, with in_hdf5, from its arrays
/// /Ret_int_K and /Ret_int_Kp, and checks them: n_tau+1 finite values of
/// each, K >= 0, and K(0) = K(beta) = 0, as K(beta - tau) = K(tau).
inline retarded_interaction read_retarded_interaction(std::string const& file, bool in_hdf5, std::size_t n_tau)
{
  using retarded_interaction_detail::error;
  const std::size_t n = n_tau + 1;
  retarded_interaction r;
  if (in_hdf5) {
    std::ifstream exists(file.c_str());
    if (!exists.good())
      throw error(file, "cannot open the file.");
    alps::hdf5::archive ar(file, "r");
    for (char const* path : {"/Ret_int_K", "/Ret_int_Kp"})
      if (!ar.is_data(path))
        throw error(file, std::string("has no array ") + path + " (K_IN_HDF5 = 1).");
    ar >> alps::make_pvp("/Ret_int_K", r.K);
    ar >> alps::make_pvp("/Ret_int_Kp", r.Kp);
    if (r.K.size() != n || r.Kp.size() != n) {
      std::ostringstream msg;
      msg << "/Ret_int_K has " << r.K.size() << " and /Ret_int_Kp " << r.Kp.size()
          << " values; with N = " << n_tau << " both need N+1 = " << n << ".";
      throw error(file, msg.str());
    }
  } else {
    std::ifstream in(file.c_str());
    if (!in.good())
      throw error(file, "cannot open the file.");
    std::vector<double> entries;
    std::string word;
    while (in >> word) {
      std::istringstream number(word);
      double x;
      if (!(number >> x) || !number.eof())
        throw error(file, "entry " + std::to_string(entries.size() + 1) + " (\"" + word + "\") is not a number.");
      entries.push_back(x);
    }
    if (entries.size() != 3 * n) {
      std::ostringstream msg;
      msg << "has " << entries.size() << " numbers; with N = " << n_tau << " it needs N+1 = " << n
          << " lines \"tau K K'\", " << 3 * n << " numbers.";
      throw error(file, msg.str());
    }
    for (std::size_t i = 0; i < n; ++i) {
      r.K.push_back(entries[3 * i + 1]);
      r.Kp.push_back(entries[3 * i + 2]);
    }
  }
  double scale = 1.;
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(r.K[i]) || !std::isfinite(r.Kp[i]))
      throw error(file, "K or K' is not finite at index " + std::to_string(i) + ".");
    scale = std::max(scale, std::abs(r.K[i]));
  }
  // rounding in a K computed from a formula
  const double tolerance = 1e-10 * scale;
  for (std::size_t i = 0; i < n; ++i)
    if (r.K[i] < -tolerance) {
      std::ostringstream msg;
      msg << "K(tau) < 0 at index " << i << " (" << r.K[i] << "); K is never negative.";
      throw error(file, msg.str());
    }
  if (std::abs(r.K[0]) > tolerance || std::abs(r.K[n - 1]) > tolerance) {
    std::ostringstream msg;
    msg << "K(0) = " << r.K[0] << " and K(beta) = " << r.K[n - 1]
        << "; both are 0, by the convention K(0) = 0 and the symmetry K(beta - tau) = K(tau).";
    throw error(file, msg.str());
  }
  return r;
}

#endif
