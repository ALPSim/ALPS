/*
 * Copyright (C) 1998-2018 ALPS Collaboration.
 * All rights reserved. Use is subject to license terms. See LICENSE.TXT
 * For use in publications, see ACKNOWLEDGE.TXT
 */

// Adapted from CQMP/Maxent src/maxent_backcont.hpp at
// a3992534228b0de732e033ae23154969724814be; see ../README.md.
#pragma once
#include "maxent.hpp"
#include <alps/config.h> // needed to set up correct bindings

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
