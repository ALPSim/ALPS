// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Generate legacy_checkpoint.h5 with the native library at ALPS revision
// f28d428017773f5794dc9896a544f95e8ef40443, never with the current SDK.
// This exercises the historical C++ serializers used by the old Python bindings;
// it does not reproduce the old bindings' Python-object conversion behavior.
//
// Reproduction (from the repository root, with Boost/HDF5/BLAS/LAPACK installed):
//   mkdir -p /tmp/alps-legacy-source
//   git archive f28d428017773f5794dc9896a544f95e8ef40443 | tar -x -C /tmp/alps-legacy-source
//   cp python/pyalps/tests/fixtures/legacy_checkpoint.cpp /tmp/alps-legacy-source/
// Append these two lines to that checkout's CMakeLists.txt:
//   add_executable(legacy_checkpoint legacy_checkpoint.cpp)
//   target_link_libraries(legacy_checkpoint alps)
// Then configure and build the historical checkout:
//   cmake -S /tmp/alps-legacy-source -B /tmp/alps-legacy-build -G Ninja \
//     -DALPS_USE_SYSTEM_BOOST=ON -DCMAKE_PREFIX_PATH="$DEPENDENCY_PREFIX" \
//     -DALPS_BUILD_PYTHON=OFF -DALPS_BUILD_APPLICATIONS=OFF -DALPS_BUILD_TESTS=OFF \
//     -DALPS_BUILD_EXAMPLES=OFF -DALPS_BUILD_FORTRAN=OFF \
//     -DALPS_ENABLE_MPI=OFF -DALPS_ENABLE_OPENMP=OFF
//   cmake --build /tmp/alps-legacy-build --target legacy_checkpoint --parallel 2
//   /tmp/alps-legacy-build/legacy_checkpoint python/pyalps/tests/fixtures/legacy_checkpoint.h5
// Original generation used GCC 15.3, Boost 1.85.0 and HDF5 1.14.6 on Linux.
#include <alps/ngs/params.hpp>
#include <alps/ngs/random01.hpp>
#include <alps/ngs/mcobservables.hpp>
#include <alps/ngs/mcresult.hpp>
#include <alps/hdf5/vector.hpp>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2) return 1;
    alps::params params;
    params["L"] = 16;
    params["T"] = 1.25;
    params["SEED"] = 42;
    params["label"] = std::string("check Ω");
    params["couplings"] = std::vector<double>{1., 2., 3.};
    alps::random01 rng(42);
    for (int i = 0; i < 7; ++i) rng();
    auto expected_rng = rng;
    std::vector<double> expected;
    for (int i = 0; i < 8; ++i) expected.push_back(expected_rng());
    alps::mcobservables observables;
    observables.create_RealObservable("Energy");
    observables.create_RealVectorObservable("Correlations");
    for (int i = 0; i < 64; ++i) {
        observables["Energy"] << double(i);
        observables["Correlations"] << std::vector<double>{double(i), 2. * i};
    }
    alps::hdf5::archive archive(argv[1], "w");
    archive["/expected_rng"] << expected;
    archive.set_context("/parameters");
    params.save(archive);
    archive.set_context("/rng");
    rng.save(archive);
    archive.set_context("/observables");
    observables.save(archive);
    archive.set_context("/result");
    alps::mcresult(observables["Energy"]).save(archive);
    archive.set_context("/vector_result");
    alps::mcresult(observables["Correlations"]).save(archive);
}
