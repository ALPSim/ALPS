// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT

#ifndef DMFT_QMC_INTERACTION_EXPANSION_CHOICE_H
#define DMFT_QMC_INTERACTION_EXPANSION_CHOICE_H

/// @file interaction_expansion_choice.h
/// @brief CT-INT (Interaction Expansion) run-type selection for the dmft
///        driver's Matsubara (OMEGA_LOOP) path.
///
/// Factored out of main.C so the (FLAVORS, SITES) -> run-type mapping is a
/// single mutually-exclusive decision that can be unit-tested. It reproduces
/// the driver's previous *effective* behaviour: FLAVORS==2 selected the
/// single-site Hubbard run; otherwise (with SITES==1) the multiband-density
/// run. The previously-built single-site half-filling factory was always
/// overwritten before use and is dropped. The remaining combination
/// (FLAVORS!=2 and SITES!=1) used to fall through to ExternalSolver with the
/// literal solver name "Interaction Expansion" (a guaranteed exec failure);
/// it is now reported explicitly via the `unsupported` choice.

enum class interaction_expansion_choice {
  single_site_hubbard,  ///< FLAVORS == 2  (HubbardInteractionExpansionRun)
  multiband_density,    ///< FLAVORS != 2 and SITES == 1 (MultiBandDensity...)
  unsupported           ///< FLAVORS != 2 and SITES != 1
};

/// Select the CT-INT run type from the flavour and site counts.
inline interaction_expansion_choice
select_interaction_expansion(int flavors, int sites)
{
  if (flavors == 2) return interaction_expansion_choice::single_site_hubbard;
  if (sites == 1)   return interaction_expansion_choice::multiband_density;
  return interaction_expansion_choice::unsupported;
}

#endif
