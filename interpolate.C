/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2009 by Sebastian Fuchs <fuchs@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
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

#include "interpolator.h"
#include "cluster.h"
#include <alps/alea.h>



  
int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif
    int i=1;
    alps::Parameters new_parms;
    while (i<argc-1 && argv[i][0]=='-') {
      new_parms[argv[i]+1]=argv[i+1];
      i+=2;
    }
    if (argc-i!=1) {
      std::cerr << "Usage: \n" << argv[0] 
		<<  "inputfile \n";
      exit(1);
    }
    alps::Parameters parms;
    std::ifstream is(argv[i]);
    if(!is.is_open()) {
      std::cerr << "parameter file argv[" << i<< " ] " <<argv[i] << "is not open! exiting!" << std::endl; 
      abort(); 
    }
    is>>parms;
    is.close();
    boost::shared_ptr<Interpolator> interptr;
    if (parms.value_or_default("ANTIFERROMAGNET", false))
      interptr.reset(new AntiferroInterpolator(parms));
    else
      interptr.reset(new ParaInterpolator(parms));
    // std::cerr << "interpolating local Green function from self-energies" << std::endl;
    // interptr->local_green_function();
    std::cerr << "interpolating self-energies and calculating dispersion" << std::endl;
    interptr->interpolate();
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& e)
    {
      std::cerr << "Caught exception: " << e.what() << "\n";
      std::exit(-5);
    }
#endif
}



