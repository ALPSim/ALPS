# C++ Ising simulation

Fill in the `...` placeholders in `ising-skeleton.cpp` to implement the spin
updates and measurements. The completed program is in `solution/ising.cpp`.

With an installed ALPS SDK, configure and build the solution from this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/alps/install
cmake --build build --config Release
```

After completing the exercise, build it with
`cmake --build build --config Release --target ising-exercise`. Run either executable from a
writable working directory; it writes simulation results as HDF5 files there.
