# Looper

ALPS Looper provides path-integral and stochastic series expansion (SSE) quantum Monte Carlo. Start with the [susceptibility](../03-mc/02-susceptibilities/) or [measurement](../03-mc/04-measurements/) tutorials for a guided simulation.

Follow the [ALPS build instructions](../../CONTRIBUTING.md#getting-started-with-the-code) with applications enabled. To rebuild just the solver:

```sh
cmake --build --preset default --target loop
```

Prepare input with `parameter2xml` or `pyalps.writeInputFiles`, then run `loop` on the generated job file. Select `ALGORITHM="loop; path integral"` or `ALGORITHM="loop; sse"`. Installed programs find the SDK's XML resources automatically; `ALPS_XML_PATH` overrides their location. MPI support is opt-in at build time.

The [original Looper manual](index.html) retains the algorithm, annealing, observable and citation reference, including Synge Todo's attribution. It is preserved as a historical document: its standalone build instructions, `REPRESENTATION` selector and `$PREFIX/lib/xml` paths describe an older release. Use the current build and input instructions above. See [Citing ALPS](../../CITATION.md) for citation guidance.
