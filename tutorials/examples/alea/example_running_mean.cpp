/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* Copyright (C) 2011-2012 by Lukas Gamper <gamperl@gmail.com>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Maximilian Poprawe <poprawem@ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/alea.h>
#include <alps/alea/mcanalyze.hpp>
#include <alps/utility/encode.hpp>

#include <alps/hdf5.hpp>

#include <iostream>
#include <string>


// This is an example of how to calculate the running mean and the reverse running mean of data stored in a hdf5 file.

int main() {

  const std::string filename = "testfile.h5";

  // create mcdata object with the correct template parameter.
  alps::alea::mcdata<double> obs;

  // load the variable E saved in the file testfile.h5 into the mcdata object.
  obs.load(filename, "simulation/results/" + alps::hdf5_name_encode("E"));

  // calculate the running mean
  alps::alea::mctimeseries<double> running_mean = alps::alea::running_mean(obs);

  // calculate the reverse running mean
  alps::alea::mctimeseries<double> reverse_running_mean = alps::alea::reverse_running_mean(obs);

  // write to std::cout
  std::cout << "The running mean of E is: " << running_mean << "\n\n";
  std::cout << "The reverse running mean of E is: " << reverse_running_mean << "\n\n";

  return 0;
}
