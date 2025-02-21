/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* Copyright (C) 2011-2012 by Lukas Gamper <gamperl@gmail.com>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Maximilian Poprawe <poprawem@ethz.ch>
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

#include <alps/alea.h>
#include <alps/alea/mcanalyze.hpp>
#include <alps/utility/encode.hpp>
#include <alps/utility/size.hpp>

#include <alps/hdf5.hpp>

#include <iostream>
#include <string>


// This is an example of how to easily calculate and fit the autocorrelation as well as estimate the integrated autocorrelation time of data stored in a hdf5 file.

int main() {
  // this using-statement makes life easier when entering the named parameters
  using namespace alps::alea;
  using alps::size;
 
  const std::string filename = "testfile.h5";

  // create mcdata object with the correct template parameter.
  mcdata<double> obs;

  // load the variable m saved in the file testfile.h5 into the mcdata object.
  obs.load(filename, "simulation/results/" + alps::hdf5_name_encode("m"));

  // calculate the full autocorrelation (with size - 1 lags)
  mctimeseries<double> auto_corr = alps::alea::autocorrelation(obs, _distance = (size(obs) - 1));
  
  // fit the autocorrelation exponentially between the values where it is at 80% and at 20% of the value at t = 1
  std::pair<double, double> fit = exponential_autocorrelation_time(auto_corr, _min=.2, _max=.8);

  // calculate the integrated autocorrelation time by summing up the autocorrelation until a specific point
  //(here when the autocorrelation reaches 20% of its value at t = 1) and integrating the fit for the tail
  double int_autocorr_time = integrated_autocorrelation_time(cut_tail(auto_corr, _limit=0.2), fit);

  // write to std::cout
  std::cout << "The autocorrelation of m is: " << auto_corr << "\n\n";
  std::cout << "The exponential fit is: " << fit.first << " * e^( " << fit.second << " * t)\n\n";
  std::cout << "The estimated integrated autocorrelation time is: " << int_autocorr_time << "\n";

  return 0;
}

