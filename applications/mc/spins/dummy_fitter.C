/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2005 by Matthias Troyer <troyer@comp-phys.org>,
*                       Andreas Streich <astreich@student.ethz.ch>
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

/* $Id$ */

#include "fitter.h"

/**
 * Constructer, reads the experimental data.
 */
DummyFitter::DummyFitter(const std::string& filename)
  : AbstractFitter(filename)
{
  alps::ParameterList init_params;
  {
    std::ifstream param_stream(filenamebase.c_str());
    param_stream >> init_params;
  }
  FitterParamList initial;
  initial.parms = adjust_parameters(init_params, all_exp_res);
  initial.is_genuine = true;
  next_parameters.push_back(initial);
}

/**
 * Updates the parameters and writes a new parameter file.
 *
 * @params basename the base file name for the parameter file.
 *             The output file will have the name 
 *             <basename>.step<step>.in.xml
 * @params step the step counter.
 */
void DummyFitter::find_new_parameters()
{ 
  if (num_variable_params == 0) {
    std::cerr << "no variable parameters !!\n";
    return;
  }
  std::string var_par_name = variable_params.front().name;
  FitterParamList res;
  res.parms = update_parameters(var_par_name, 
               (double)((current_params.parms)[0][var_par_name])+0.1, 
               current_params.parms);
  
  res.is_genuine = true;
  next_parameters.push_back(res);
}

/**
 * Checks whether the fitting is done.
 */
bool DummyFitter::is_done(double* seconds) {
  *seconds = -1.;
  return (step >= max_steps);
}

/**
 * Return the name of the fitter.
 */
std::string DummyFitter::get_fitter_name() const
{
  return std::string("Dummy Fitter");
}
