# Parameter tests

Checkpoint tests preserve values, variant types and key order. Integer storage
widths are accepted when values fit the native `int` alternative; overflow and
unsupported types must leave existing state intact. Custom readers own their
decoding policy. Serialization comparisons are exact, without numeric tolerance.

Legacy text/XML parsing is covered separately in
`tests/integration/params/legacy_adapters.cpp` and requires `ALPS::alps`.
