# Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
set(ALPS_MAXENT_SOURCES maxent_helper.cpp maxent_simulation.cpp maxent_parms.cpp)
list(TRANSFORM ALPS_MAXENT_SOURCES PREPEND "${CMAKE_CURRENT_LIST_DIR}/")
