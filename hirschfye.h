#include "solver.h"
#include <boost/filesystem/path.hpp>
#include <cstdio>

class HirschFyeSolver : public ImpuritySolver {
public:
  HirschFyeSolver() {}
  multiple_vector_type solve(const multiple_vector_type& G0, const alps::Parameters& parms);
}
