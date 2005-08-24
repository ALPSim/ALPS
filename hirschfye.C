#include "hirschfye.h"
#include "hirschfyesim.h" // your old hirschfye.h
#include <cstdio>


HirschFyeSolver::multiple_vector_type HirschFyeSolver::solve(const multiple_vector_type& G0, const alps::Parameters& parms)
{
  alps::ProcessList where(1,alps::local_process());
  HirschFyeSimulation sim(where,parms,0); // your old IsingSimulation
  sim.set_GO(G0); // to implement
  
  while (sim.work_done() <1.)
    sim.run();

  multiple_vector_type g = sim.get_G(); // to implement
  
  // .... maybe write a checkpoint ....
  return g;
}


