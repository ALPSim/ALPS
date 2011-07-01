#include <alps/parapack/parapack.h>
#include "alps/fortran/fortran_wrapper.h"

PARAPACK_SET_VERSION("my version");
PARAPACK_SET_COPYRIGHT("my copyright");
PARAPACK_REGISTER_WORKER(alps::fortran_wrapper, "ising");
PARAPACK_REGISTER_EVALUATOR(alps::parapack::simple_evaluator, "ising");

int main(int argc, char** argv) { return alps::parapack::start(argc, argv); }
