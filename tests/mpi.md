# Native MPI tests

MPI correctness tests are GoogleTest executables launched by CTest with real
communicators. Configure and build the MPI preset before selecting the MPI label:

```sh
cmake --preset mpi
cmake --build --preset mpi --parallel 2
ctest --preset mpi -L '^mpi$' --output-junit mpi.xml
```

Each executable runs at two and three ranks. These topologies fit the smallest
supported CI runners and exercise both paired communication and participation by
an additional rank. `comm_mpi` uses groups of two: its three-rank run keeps the
controller outside the worker communicator. The launcher needs at least three
slots; on smaller local allocations, configure `MPIEXEC_PREFLAGS=--oversubscribe`
when using OpenMPI.

## Executed contracts

| Component / executable | Assertions |
| --- | --- |
| `legacy_parameters / parameters_mpi` | Ordered keys, decimal values, and unevaluated parameter expressions survive a relay through every rank and back to the sender. |
| `alea / observableset_mpi` | The historical seed-2873 sequence of 4,096 paired samples, observable names, counts, means, errors, convergence states, and derived ratio survive the same relay. Every rank constructs its reference independently of serialization. |
| `parapack / collect_mpi` | Uneven vector lengths and offsets distribute correctly and reconstruct the original seeded vector. |
| `parapack / comm_mpi` | Control, worker, and group-head membership, communicator sizes, local ranks, and the corresponding allocated process list. |
| `parapack / process_mpi` | Group identifiers and members, free/active accounting, allocation exhaustion, recycling, release, and rejection of a second release. |
| `parapack / halt_mpi` | Shutdown waits for active groups to be released, receives every rank's acknowledgement, and remains halted on repeated calls. |
| `parapack / filelock_mpi` | Every rank acquires in turn; other ranks fail a bounded single lock attempt while it is held. Explicit release and destructor release both permit the next owner. |
| `parapack / info_test_mpi` | Reproducible per-rank worker seeds, common disorder seed, rank-specific checkpoints, master-only phase/host metadata, and progress updates. |
| `scheduler / scheduler_sum_mpi` | A deterministic completion schedule takes 128 scalar/vector Monte Carlo measurements per rank; global counts, independently summed means, result arithmetic, and HDF5 persistence agree. |

## Execution and reports

`alps_add_mpi_gtest` obtains the launcher and its flags from CMake's FindMPI
configuration (`MPIEXEC_EXECUTABLE`, `MPIEXEC_NUMPROC_FLAG`, `MPIEXEC_PREFLAGS`,
`MPIEXEC_POSTFLAGS`). The shared main initializes/finalizes MPI, aggregates rank
failures, and rejects a zero-case selection on any rank. An intentional empty
GoogleTest filter is checked to return a failing launcher status.

CTest records each complete topology as an integration test labelled `mpi` and
with a `PROCESSORS` requirement equal to the rank count. Each topology has a
120-second timeout; process shutdown assertions also have a 10-second deadline.
OpenMP and BLAS thread counts are set to one. Tests use matching collectives and
complete required communication before rank-specific fatal assertions.

Each topology owns a separate working directory. GoogleTest writes
`test-work/<component>.<target>.np<ranks>/gtest-rank-<rank>.xml`, while CTest's JUnit
report records the launcher result and timeout failures. CI retains both,
alongside CTest logs and any actual-output diagnostics. The file-lock case shares
one unique temporary directory between local ranks and synchronizes before
cleanup. It assumes the ranks share that filesystem, as they do in the hosted
CI jobs; it is not a distributed-filesystem validation.

## Inherited manual drivers

`parapack/tests/clone_mpi.C` and `worker_mpi.C` remain unbuilt manual simulation
and checkpoint/restart drivers. Neither was registered in the previous CMake
suite. They rely on historical Ising tutorial workers, factory registration,
stdin configuration, and progress-driven simulation loops; `worker_mpi.C` also
uses obsolete include paths. Their old input/output files remain with them.

The tests above do not establish complete Ising-worker restart equivalence or
clone-proxy scheduling behavior. A future integration test for those contracts
should supply fixed inputs, isolate checkpoints, bound progress, and compare
resumed measurements/state against an uninterrupted reference. These sources
are an explicit coverage gap, not skipped passing tests.
