#include <alps/parameters.h>
#include "selfconsistency.h"

int main(int argc, char** argv)
{

  alps::Parameters parms;
  std::cin >> parms;
  
  boost::shared_ptr<ImpuritySolver> solver_ptr;
  if (parms["SOLVER"]=="Hirsch-Fye")
    solver_ptr.reset(new HirschFyeSolver(parms));
  else 
    solver_ptr.reset(new ExternalSolver(boost::filesystem::path(parms["SOLVER"],boost::filesystem::native)));
  
  SemicircleHilbertTransformer transform(2.*sqrt(2.)*double(parms["t"]));

  selfconsistency_loop(parms,*solver_ptr,transform);
  
}
