# Validation, 21 October 2026

## Clean public-checkout installation

[clean-install-20261021.json](clean-install-20261021.json) records the fresh
public ALPS/QDK checkouts, isolated Python 3.13 environment, rebuilt native
code, dependency versions and executed commands. Follow
[INSTALL.md](../INSTALL.md) to reproduce it. Third-party native dependencies
were preinstalled and reused; their cold provisioning was not repeated.

The adapter suite passed **95 tests**, with **6 RdmContainer-only tests skipped**.
The pinned QDK provides MPSContainer/MPSSite but not RdmContainer. All **14
top-level examples** passed, as did N₂ preparation, the 8/16/32-state solve,
plotting, 8-bit Q# preparation and QPE resource counting. Upstream QDK MPS
container tests passed **8 tests** and the citation tests passed **16 tests**.
`pip check` passed; all imported solver packages and QDK resources resolved
inside the new environment with no Python/resource overlay.

The six-site Hubbard ED/DMRG energy difference was 2.66e-15. The N₂ DMRG
overlap at 16 states agrees with the earlier result to numerical precision.
Finite-precision Q# preparation gave **99.696%** ground-state probability,
compared with **99.730%** in the earlier mixed runtime. The presentation now
uses the clean-run preparation number. The new coarse QPE resource count is
57 logical qubits, 6,014,109 rotations and 8,276 CCZ gates; QPE is not executed.
These are version-specific observations, not bitwise reproducibility claims.

The macOS watchdog initially encountered the execution sandbox's process-
inspection restriction. Both optional commands subsequently passed in a normal
process-inspection context. Its `finally` cleanup terminated the first child.
Run guarded examples from an ordinary terminal.

## Original presentation measurements

[20261021.json](20261021.json) records the freshly executed example results,
controls, dependency snapshots, input hashes and final source-file hashes used
for the accompanying scientific presentation. Measurements used an existing
Apple-silicon pyalps build and the existing QDK environments; no solver or
third-party dependency was rebuilt for this port.

The adapter suite passed with both supported development configurations:
95 passed / 6 skipped with MPS-container support, and 94 passed / 7 skipped
with RdmContainer support. Skips cover capabilities absent from the respective
QDK build. All 14 top-level calculation examples ran successfully. Molecular
preparation, nitrogen bond scans, MPS serialization, Q# state preparation,
QPE resource counting, capacity/scaling pilots and both renderer input formats
were also exercised. The six-orbital Block2 reference passed in its own process.

The default QPE example uses eight phase bits and four Trotter subdivisions:
57 logical qubits, 6,014,109 arbitrary rotations and 8,242 CCZ gates. Counting
took 4.27 s and peaked at 0.814 GB RSS. Its 49.1 mEh phase bins are a coarse
demonstration, not an energy-accuracy certificate. The original 14-bit,
32-subdivision scenario exceeded physical-footprint caps of both 12 and 16 GB.
The eight-bit MPS preparation was simulated; the ten-bit simulation exceeded
180 s and then 420 s. No accuracy result is claimed for the timed-out run.

Timings are single local observations, with other machine activity uncontrolled.
They are not repeated comparative benchmarks or asymptotic scaling fits. Peak
RSS is cumulative per process; peak scratch was not measured. Larger-system
results establish capacity and final-state/RDM consistency, not a converged
singlet ground state. In particular, the four-sweep CAS(10,26) pilot has
⟨S²⟩ ≈ 1.61.

Raw local logs and generated binary files were retained outside the repository
in `/private/tmp/qdk-classic-examples-20261021`. Regenerate them with the commands
in the parent example READMEs; the source branch contains no copied historical
figures, compiled libraries or large wavefunction binaries.
