/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
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

/* $Id$ */

#include "selfconsistency.h"
#include "hirschfye.h"
#include "externalsolver.h"
#include "hilberttransformer.h"
#include <alps/parameters.h>
#include <alps/copyright.h>

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

  std::cout << "ALPS DMFT solver for the single site impurity problem.\n\n";
  alps::print_copyright(std::cout);

  alps::Parameters parms;
  std::cin >> parms;
  
  boost::shared_ptr<ImpuritySolver> solver_ptr;

  if (parms["SOLVER"]=="Hirsch-Fye")
    solver_ptr.reset(new HirschFyeSolver());
  else 
    solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"],boost::filesystem::native)));
  
  SemicircleHilbertTransformer transform(2.*sqrt(2.)*double(parms["t"]));

  selfconsistency_loop(parms,*solver_ptr,transform);
#ifndef BOOST_NO_EXCEPTIONS
  }
  catch (std::exception& exc) {
      return -1;
    }
  catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
#endif  
  return 0;  
}
