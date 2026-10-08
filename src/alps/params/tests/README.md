# Parameter test contracts

The seven existing targets are GoogleTest suites, discovered separately by
CTest under `params.<target>.*`. File-writing cases each own a temporary directory.

| Legacy target | Preserved and strengthened checks |
| --- | --- |
| `param_default` | Missing-value fallback, existing-value precedence, and non-inserting lookup. |
| `param_assign` | All former scalar assignments are read back; boolean vectors and strings retain values; `find` distinguishes present and missing keys. Assertions remain enabled in Release builds. |
| `param_ordering` | Insertion order before and after HDF5 storage, plus all parameter values. |
| `param_stream` | Exact string rendering previously only printed to stdout. |
| `param_not_found` | Mutable and const missing lookups throw and name the key. The old three-line `.output` contract is now explicit assertions. |
| `param_external` | Shared external source lifetime, native Boost serialization, overflow rejection, string-list conversion, heterogeneous HDF5 lists and repeated load/save. These are three independently discovered cases. |
| `param_checkpoint` | Supported native variants and values, unsupported scalar/array storage, preservation of existing values/order on failed reload, and custom-reader behavior across reload. These are seven independently discovered cases. Storage-width-dependent exclusions are visible GoogleTest skips. |

`tests/integration/params/legacy_adapters.cpp` has separate cases for legacy text
parsing/conversion (including invalid grammar) and XML parsing with default or
explicit seeds. No existing contract was replaced with a check of a wrapped
program's exit status. Numerical checkpoint equality is exact because these
operations must preserve the stored value, with no approximate computation.
