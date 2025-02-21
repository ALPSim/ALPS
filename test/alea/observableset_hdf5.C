/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010-2012 by Lukas Gamper <gamperl -at- gmail.com>,
*                            Synge Todo <wistaria@comp-phys.org>
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

#include <alps/utility/encode.hpp>
#include <alps/alea.h>

#include <alps/hdf5.hpp>

#include <boost/filesystem.hpp>
#include <boost/random.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

int main() {
  try {
    std::string const filename = "observableset.h5";
    if (boost::filesystem::exists(boost::filesystem::path(filename)))
        boost::filesystem::remove(boost::filesystem::path(filename));
    {
        alps::ObservableSet measurement;
        measurement << alps::make_observable(alps::RealObservable("Test"), true)
                    << alps::RealObservable("Sign")
                    << alps::RealObservable("No Measurements")
                    << alps::IntHistogramObservable("Histogram", 0, 10)
                    << alps::RealObservable("Test 2")
                    << alps::RealObservable("Test 3")
                    << alps::SimpleRealObservable("Test 4");
        alps::hdf5::archive oar(filename, "a");
        oar["/test/0/result"] << measurement;
    }
    {
        boost::minstd_rand0 engine;
        boost::uniform_01<boost::minstd_rand0> random(engine);
        alps::ObservableSet measurement;
        alps::hdf5::archive iar(filename, "r");
        iar["/test/0/result"] >> measurement;
        for (int i = 0; i < 10000; ++i) {
          measurement["Test"] << random();
          measurement["Sign"] << 1.0;
          measurement["Histogram"] << static_cast<int>(10*random());
          measurement["Test 2"] << random();
          measurement["Test 3"] << random();
          measurement["Test 4"] << random();
        }
        alps::RealObsevaluator e2 = measurement["Test 2"];
        alps::RealObsevaluator e4 = measurement["Test 3"];
        alps::RealObsevaluator ratio("Ratio");
        ratio = e2 / e4;
        measurement.addObservable(ratio);
        
        alps::hdf5::archive oar(filename, "a");
        oar["/test/0/result"] << measurement;
        
        alps::IntHistogramObsevaluator eval = measurement["Histogram"];
    }
    {
        alps::ObservableSet measurement;
        measurement << alps::make_observable(alps::RealObservable("Test"), true)
                    << alps::RealObservable("Sign")
                    << alps::RealObservable("No Measurements")
                    << alps::IntHistogramObservable("Histogram", 0, 10);
        alps::hdf5::archive iar(filename, "r");
        iar["/test/0/result"] >> measurement;
        std::cout << measurement;
        alps::IntHistogramObsevaluator eval = measurement["Histogram"];
    }
    {
        alps::ObservableSet measurement;
        alps::hdf5::archive iar(filename, "r");
        iar["/test/0/result"] >> measurement;
        std::cout << measurement;
        alps::IntHistogramObsevaluator eval = measurement["Histogram"];
    }
    boost::filesystem::remove(boost::filesystem::path(filename));
  }
  catch (std::exception& e) {
    std::cerr << "Fatal error: " << e.what() << "\n"; 
  }
  return 0;
}
