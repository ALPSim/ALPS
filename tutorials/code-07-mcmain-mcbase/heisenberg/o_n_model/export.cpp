#define PY_ARRAY_UNIQUE_SYMBOL ndsim_PyArrayHandle

#include "ndim_spin.hpp"

#include <alps/ngs/detail/export_sim_to_python.hpp>

BOOST_PYTHON_MODULE(ndsim_c) {
    ALPS_EXPORT_SIM_TO_PYTHON(xy_sim, ndim_spin_sim<2>);
    ALPS_EXPORT_SIM_TO_PYTHON(heisenberg_sim, ndim_spin_sim<3>);
    ALPS_EXPORT_SIM_TO_PYTHON(4d_sim, ndim_spin_sim<4>);
    ALPS_EXPORT_SIM_TO_PYTHON(5d_sim, ndim_spin_sim<5>);
}
