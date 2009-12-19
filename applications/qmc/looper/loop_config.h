/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
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
