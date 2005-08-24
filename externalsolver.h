#include "solver.h"
#include <boost/filesystem/path.hpp>
#include <cstdio>

class ExternalSolver : public ImpuritySolver {
public:
  ExternalSolver(const boost::filesystem::path& executable) 
    : exe_(executable) {}

  multiple_vector_type solve(const multiple_vector_type& G0, const alps::Parameters& parms);
    
private:
  boost::filesystem::path exe_;
}
