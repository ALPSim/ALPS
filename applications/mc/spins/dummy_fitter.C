/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2005 by Matthias Troyer <troyer@comp-phys.org>,
*                       Andreas Streich <astreich@student.ethz.ch>
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
