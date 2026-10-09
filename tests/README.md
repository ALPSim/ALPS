# Native tests

Configure with `ALPS_BUILD_TESTING=ON`, build, and run `ctest --test-dir <build> --output-on-failure`. Existing Boost tests and transcript fixtures live with their native modules. Cross-module tests live under `tests/`. MPI cases are labeled `mpi`.
