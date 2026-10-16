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

// Adapted from CQMP/Maxent src/maxent_backcont.hpp at
// e727650a1f0a427e95efebc409d907b4c84011ec; see ../README.md.
#pragma once
#include "maxent.hpp"

///This class contains necessities to analytically continue a 
//real function to the imaginary axis
class Backcont{
  public:
    Backcont(const MaxEntParameters *param_in);
    ///Backcontinue a given function A
    vector_type backcontinue(const vector_type &A);
    ///determine the maximum difference between the two functions
    double max_error(const vector_type &y1, const vector_type &y2);
  private:
   //pointer to Parameters, set up through MaxEntSimulation->MaxEntHelper
   const MaxEntParameters *param;
   ///holds backcontinued Green's function
   vector_type G;
   kernel_type k_type;
};
