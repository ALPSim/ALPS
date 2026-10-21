# Validation snapshot, 21 October 2026

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
