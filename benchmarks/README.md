# Linear interpolation benchmarks

This is a reproducibility harness, not a performance claim. Run it from the
locked development environment:

```sh
pixi install --locked
pixi run bench-linear
```

The wrapper records UTC time, CPU, OS, architecture, Mojo version, compiler
options, and the exact Pixi command. The Mojo runner then prints a versioned
schema, warmup and sample counts, one ordered line per case, elapsed nanoseconds,
and a deterministic result checksum.

## Matrix

The baseline covers 26 cases:

- construction for 8, 1,024, and 65,536 knots;
- repeated scalar evaluation for the same sizes under `ERROR`, `CLAMP`, and
  `LINEAR`;
- uniform and deterministic irregular knot layouts at every size;
- in-domain queries for `ERROR` and alternating below/above-domain queries for
  `CLAMP` and `LINEAR`;
- tiny-ordinate linear extrapolation and nearly full-range in-domain evaluation
  as explicit extreme-arithmetic cases.

Construction timing includes copying the input lists, validation, allocation,
destruction of the previous interpolator, and one `std.benchmark.keep` barrier
per construction so overwritten instances cannot be eliminated. Evaluation
timing uses one constructed interpolator, but each public `evaluate` call
intentionally includes the current `O(n)` table revalidation followed by
interval search. Query generation and checksum accumulation are also inside the
measured loop; they must remain unchanged when comparing two revisions.

The runner performs two complete warmup rounds and seven measured samples. It
reports the minimum total elapsed nanoseconds, not a derived throughput claim.
Divide `best_elapsed_ns` by `iterations` only when comparing identical cases on
equivalent machines and builds. Checksums are retained through
`std.benchmark.keep` to prevent dead-code elimination.

## Determinism and interpretation

`bench_linear.mojo --manifest` emits the output schema without taking timings.
The locked test suite checks that it has exactly 26 ordered identities and that
the construction, policy, scale, and extreme sentinels remain present. Fixture
reference, ordering, mutation, and extreme-number contracts are covered in
`tests/test_benchmark_fixtures.mojo`.

Raw timings vary with power state, thermal state, background load, compiler,
and hardware. Do not commit raw output as a marketing comparison. Record at
least three complete runs, preserve every metadata header, and explain any
methodology change before using this harness to justify optimization.
