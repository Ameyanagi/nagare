# Benchmark datasets

NAG-004 uses no downloaded or generated files. Every table is built in memory
by the committed pure-Mojo generator in `benchmarks/linear_fixtures.mojo`, so
the repository commit identifies the exact dataset implementation and there is
no external license or checksum to track.

For moderate-magnitude tables, the generator starts at `x[0] = -2`. Uniform
steps are `0.75`; irregular step `i` is
`0.75 + ((17 * i) mod 5) * 0.125`. Ordinates are
`1.5 * x[i] - 2 + (((7 * i) mod 3) - 1) * 0.25`. The committed sizes are 8,
1,024, and 65,536 knots. The small uniform table has an independently stated
four-element prefix in the fixture tests.

Two two-knot extreme fixtures exercise numerical branches already guaranteed
by the library contract:

- `[0, 1e-308] -> [-1e-308, 1e-308]`, extrapolated at `x = 1`;
- `[-MAX_FINITE, +MAX_FINITE] -> [-MAX_FINITE, +MAX_FINITE]`, evaluated at
  `x = 1` under the documented sub-ULP central-offset limitation.

The result checksum printed for each benchmark case is a deterministic semantic
guard, not a cryptographic dataset checksum. Any future external or committed
dataset must add its canonical source, license, retrieval date, SHA-256,
generation command, and expected dimensions here.
