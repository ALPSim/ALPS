# ALPSCore archive fixture

`archive_7146b9e1.h5` was written by the archive of
[ALPSCore](https://github.com/ALPSCore/ALPSCore) at commit 7146b9e1, using
`writer.cpp` in this directory. `hdf5_alpscore_fixture` reads every value back
with ALPS's archive. The values are listed once, in `items.hpp`.

The file covers what the two implementations store differently: ALPSCore
writes bool and signed char without ALPS's `__alps_type__` tag, spells
attributes `path@name` rather than `path/@name`, and marks `boost::optional`
with an `alps_hdf5_optional_empty` attribute. Keep the file frozen; a copy
regenerated with ALPS's own archive would no longer test compatibility.

## Regenerating and checking both directions

`writer.cpp` is not part of the ALPS build. With ALPSCore's `utilities` and
`hdf5` modules installed under `$ALPSCORE` (its CMake needs
`-DEIGEN3_INCLUDE_DIR` for Eigen 5), build it against ALPSCore:

```sh
c++ -std=c++17 -I"$ALPSCORE/include" -I"$EIGEN" writer.cpp -o writer_core \
  -L"$ALPSCORE/lib" -lalps-hdf5 -lalps-utilities -lhdf5
./writer_core write archive_7146b9e1.h5
```

Built the same way against an ALPS SDK (`ALPS::hdf5`), `writer write` makes
an ALPS file, and `writer_core read <file>` checks that ALPSCore reads it.
On 2026-10-13 both directions read every value correctly, and `h5dump -p`
of the two files differed only by ALPS's `__alps_type__` tags.
