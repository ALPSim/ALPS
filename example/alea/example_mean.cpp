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

#include <alps/hdf5.hpp>

#include <iostream>
#include <string>


// This is an example of how to easily calculate the mean of data stored in a hdf5 file.

int main() {

  const std::string filename = "testfile.h5";

  // create mcdata object with the correct template parameter.
  alps::alea::mcdata<double> obs;

  // load the variable E saved in the file testfile.h5 into the mcdata object.
  obs.load(filename, "simulation/results/" + alps::hdf5_name_encode("E"));

  // calculate the mean
  double mean = alps::alea::mean(obs);

  // write to std::cout
  std::cout << "The mean of E is: " << mean << std::endl;

  // write the result back to the file
  alps::hdf5::archive ar(filename, "a");
  ar << alps::make_pvp("simulation/results/" + alps::hdf5_name_encode("E") + "/mean/value", mean);

  return 0;
}
