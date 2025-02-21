/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifdef LOOP_CONFIG_HEADER
# include LOOP_CONFIG_HEADER
#else

#ifndef LOOP_CONFIG_H
#define LOOP_CONFIG_H

#include <alps/lattice.h>
#include <looper/graph.h>
#include <looper/model.h>

// measurements
#include <looper/correlation.h>
#include <looper/custom.h>
#include <looper/stiffness.h>
#include <looper/susceptibility.h>

struct loop_config {
  // lattice structure
  typedef alps::coordinate_graph_type lattice_graph_t;
  typedef looper::lattice_helper<lattice_graph_t> lattice_t;

  // imaginary time
  typedef double time_t;

  // graph for loops
  typedef looper::local_graph<> loop_graph_t;

  // model
  typedef looper::spinmodel_helper<lattice_graph_t, loop_graph_t> model_t;

  // whether longitudinal external field is supported or not
  static const bool support_longitudinal_field = true;

  // whether systems with negative signs are supported or not
  static const bool support_negative_sign = true;

  // measurements
  typedef looper::measurement_set<
    looper::correlation,
    looper::custom_measurement,
    looper::stiffness<3>,
    looper::susceptibility
  > measurement_set;
};

#endif // LOOP_CONFIG_H

#endif // LOOP_CONFIG_HEADER
