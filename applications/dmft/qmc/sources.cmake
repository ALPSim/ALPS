# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
# Shared by the command-line programs and standalone Python source builds.
set(ALPS_CTHYB_SOURCES
    hybridization/hybsim.cpp
    hybridization/hyblocal.cpp
    hybridization/hybint.cpp
    hybridization/hybfun.cpp
    hybridization/hybretintfun.cpp
    hybridization/hybmatrix.cpp
    hybridization/hybmatrix_ft.cpp
    hybridization/hybconfig.cpp
    hybridization/hybupdates.cpp
    hybridization/hybevaluate.cpp
    hybridization/hybmeasurements.cpp)
set(ALPS_CTINT_SOURCES
    fouriertransform.C
    interaction_expansion2/auxiliary.cpp
    interaction_expansion2/observables.cpp
    interaction_expansion2/fastupdate.cpp
    interaction_expansion2/selfenergy.cpp
    interaction_expansion2/solver.cpp
    interaction_expansion2/io.cpp
    interaction_expansion2/splines.cpp
    interaction_expansion2/interaction_expansion.cpp
    interaction_expansion2/measurements.cpp
    interaction_expansion2/model.cpp)
list(TRANSFORM ALPS_CTHYB_SOURCES PREPEND "${CMAKE_CURRENT_LIST_DIR}/")
list(TRANSFORM ALPS_CTINT_SOURCES PREPEND "${CMAKE_CURRENT_LIST_DIR}/")
