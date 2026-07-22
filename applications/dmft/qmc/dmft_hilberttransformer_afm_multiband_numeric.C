// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Numeric lock for GeneralFSHilbertTransformer::operator() on the
// multiband antiferromagnetic path.
//
// The AFM self-consistency loop iterated `for(f=0; f<nflavor()/2; f+=2)`
// with 2*f / 2*f+1 indexing. The nflavor()/2 bound combined with the
// stride-2 step means that for nflavor==4 only f==0 runs (flavour pair
// 0,1); f==2 fails 2<2 and pair 2,3 is silently skipped, leaving those
// flavours at their incoming G0. The fix iterates `for(f=0; f<nflavor();
// f+=2)` with f / f+1 indexing (the SemicircleFS idiom) and is identical
// at nflavor==2.
//
// Multiband AFM is reachable via a DOSFILE (DOSBandstructure) — the
// lattice bandstructures throw on nflavor!=2, but DOSBandstructure allows
// FLAVORS/2 bands. This test builds two identical bands and feeds an input
// that is identical across the two flavour pairs; the transform must then
// produce identical output for both pairs. The pre-fix code leaves pair
// (2,3) untransformed, so its output differs from pair (0,1); the fixed
// code is symmetric.

#include "hilberttransformer.h"

#include <alps/parameter.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <boost/filesystem.hpp>
#include <fstream>
#include <string>

int main() {
  // Temporary 2-band DOS file with identical bands. 11 points (odd, required
  // by the Simpson integration). Columns: eps0 dos0 eps1 dos1.
  const std::string dosfile =
      (boost::filesystem::temp_directory_path() /
       "alps_dmft_afm_multiband_dos.dat").string();
  {
    std::ofstream f(dosfile);
    const double es[11]  = {-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0};
    const double dos[11] = { 0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.5,0.4,0.3,0.2,0.1};
    for (int i = 0; i < 11; ++i)
      f << es[i] << " " << dos[i] << " " << es[i] << " " << dos[i] << "\n";
  }

  alps::Parameters parms;
  parms["ANTIFERROMAGNET"] = 1;
  parms["FLAVORS"]         = 4;   // two bands
  parms["DOSFILE"]         = dosfile;

  GeneralFSHilbertTransformer transform(parms);

  const unsigned nfreq = 60;
  const double   mu = 0.5, h = 0.1, beta = 10.0;

  matsubara_green_function_t G_omega(nfreq, 1, 4);
  matsubara_green_function_t G0_omega(nfreq, 1, 4);
  for (unsigned w = 0; w < nfreq; ++w) {
    const std::complex<double> iw(0.0, (2 * w + 1) * M_PI / beta);
    const std::complex<double> z_up = 1.0 / (iw + 0.3);
    const std::complex<double> z_dn = 1.0 / (iw - 0.3);
    const std::complex<double> g_up = 1.0 / (iw + 0.5);
    const std::complex<double> g_dn = 1.0 / (iw - 0.5);
    // flavour pair (0,1) is set identical to pair (2,3)
    G_omega(w, 0) = z_up;  G_omega(w, 2) = z_up;
    G_omega(w, 1) = z_dn;  G_omega(w, 3) = z_dn;
    G0_omega(w, 0) = g_up; G0_omega(w, 2) = g_up;
    G0_omega(w, 1) = g_dn; G0_omega(w, 3) = g_dn;
  }

  const matsubara_green_function_t out = transform(G_omega, G0_omega, mu, h, beta);

  boost::filesystem::remove(dosfile);

  double max_asym   = 0.0;  // |pair(2,3) - pair(0,1)|
  double max_change = 0.0;  // how far pair(0,1) moved from its input
  for (unsigned w = 0; w < nfreq; ++w) {
    const std::complex<double> iw(0.0, (2 * w + 1) * M_PI / beta);
    max_asym   = std::max(max_asym, std::abs(out(w, 2) - out(w, 0)));
    max_asym   = std::max(max_asym, std::abs(out(w, 3) - out(w, 1)));
    max_change = std::max(max_change, std::abs(out(w, 0) - 1.0 / (iw + 0.5)));
  }

  bool ok = true;
  if (max_asym > 1e-10) {
    std::printf("FAIL: AFM transform asymmetric across identical band pairs "
                "(max|band(2,3)-band(0,1)|=%g); flavours >=2 were skipped.\n",
                max_asym);
    ok = false;
  }
  if (max_change < 1e-6) {
    std::printf("FAIL: transform did not modify the channel (max_change=%g); "
                "test is degenerate.\n", max_change);
    ok = false;
  }
  if (ok) {
    std::printf("OK: AFM transform symmetric across band pairs "
                "(asym=%g, change=%g)\n", max_asym, max_change);
  }
  return ok ? 0 : 1;
}
