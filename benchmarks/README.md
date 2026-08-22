# Linear interpolation benchmarks

This is a reproducibility harness, not a performance claim. Run it from the
locked development environment:

```sh
pixi install --locked
pixi run bench-linear
```

The wrapper records UTC time, Git HEAD and clean/dirty state, the `pixi.lock`
SHA-256, CPU, OS, architecture, Mojo version, compiler options, and the exact
Pixi command. Run it from a Git checkout with `git` and either `sha256sum` or
`shasum` available; those are benchmark metadata tools, not package runtime
dependencies. The Mojo runner then prints a versioned schema, warmup and sample
counts, one ordered line per case, elapsed nanoseconds, and a deterministic
result checksum.

## Matrix

The baseline covers 32 cases:

- construction for 8, 1,024, and 65,536 knots;
- repeated scalar evaluation for the same sizes under `ERROR`, `CLAMP`, and
  `LINEAR`;
- uniform and deterministic irregular knot layouts at every size;
- in-domain queries for `ERROR` and alternating below/above-domain queries for
  `CLAMP` and `LINEAR`;
- tiny-ordinate linear extrapolation and nearly full-range in-domain evaluation
  as explicit extreme-arithmetic cases. The full-range case evaluates at
  `MAX_FINITE / 4096` and normalizes each result by `MAX_FINITE`, producing the
  finite non-zero checksum `244.140625` without overflowing checksum
  accumulation;
- scalar and monotone-sorted caller-buffer batches over 1,024 and 65,536
  irregular knots, using identical queries, outputs, repetitions, and checksum
  contracts;
- scalar and monotone-sorted batches whose 128 queries cluster in the last 128
  of 65,536 irregular knots, proving sparse batches do not walk from knot zero.

Construction timing includes copying the input lists, validation, allocation,
destruction of the previous interpolator, and one `std.benchmark.keep` barrier
per construction so overwritten instances cannot be eliminated. Evaluation
timing uses one construction-validated interpolator, and each public `evaluate`
call performs interval search without rescanning the stored table. Query
generation and checksum accumulation are also inside the measured loop; they
must remain unchanged when comparing two revisions. Every timed fixed query is
passed through `std.benchmark.black_box` directly into `evaluate`, preventing
the compiler from assuming or hoisting that input. Every timed result crosses
`std.benchmark.keep` to prevent dead-code elimination.

The runner performs two complete warmup rounds and 31 measured samples. It
reports nearest-rank p50 and p95 total elapsed nanoseconds. Earlier harness
versions used only seven samples and reported their minimum; short cases could
therefore round to zero at the clock's effective resolution and the statistic
favored unusually fast runs. Version 3 raises iteration counts until the
representative Apple M4 cases take multiple milliseconds and reports both
median and tail latency. Divide an elapsed value by `iterations` only when
comparing identical cases on equivalent machines and builds. Checksums are
retained through `std.benchmark.keep` to prevent dead-code elimination.

The reusable batch cases build queries and allocate output before timing. Each
timed repetition performs the documented caller-buffer API and crosses a
`std.benchmark.keep` barrier before the next overwrite; checksum verification
occurs after the timer. This isolates interpolation/search/output writes while
preventing the compiler from collapsing repeated calls, and still checks
identical scalar and sorted numerical results.

## Determinism and interpretation

`bench_linear.mojo --manifest` emits the output schema without taking timings.
The locked test suite checks exactly 32 ordered identities together with each
numeric knot count, iteration count, checksum contract, and expected checksum.
Warmup and measured checksums must equal the untimed expectation exactly or the
runner raises. Fixture reference, ordering, mutation, and extreme-number
contracts are covered in `tests/test_benchmark_fixtures.mojo`.

Raw timings vary with power state, thermal state, background load, compiler,
and hardware. Do not commit raw output as a marketing comparison. Record at
least three complete runs, preserve every metadata header, and explain any
methodology change before using this harness to justify optimization. Prefer a
`git_state=clean` run; if local experiments require a dirty tree, preserve that
marker rather than presenting the output as a committed revision.

## Statistical CPU profiling

Latency measurements show how much time changed; sampling profiles identify
where compiled CPU time is spent. Run the dedicated steady-state driver with:

```sh
pixi run profile-linear
```

On macOS this builds an optimized binary and attaches `/usr/bin/sample` for
three seconds to dense and sparse scalar/sorted workloads. On Linux it records
call stacks with `perf record --call-graph dwarf` and renders a text report with
`perf report`. Raw reports and workload checksums are written beneath
`.pixi/profiles/nagare-linear/` and are intentionally untracked. See
[`PROFILE.md`](PROFILE.md) for the most recent engineering profile and the SIMD
acceptance decision.
