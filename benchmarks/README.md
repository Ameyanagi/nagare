# Benchmarks

Nagare now has a correctness-first scalar interval locator and linear
interpolator. No timing result is published yet: NAG-004 adds the first
benchmark harness only after the numerical contract is stable.

The first matrix must measure construction separately from repeated evaluation
for small, medium, and large validated knot tables. It must include uniform and
irregular knots, in-domain queries, every extrapolation policy, and
large-magnitude finite data that exercises both direct-delta and
scale-normalized arithmetic.

Every run records the CPU, OS, Mojo version, compiler options, dataset
provenance, warmup, measured iterations, statistic, checksum, and exact command.
A checksum of evaluated values prevents dead-code elimination. Benchmark
programs belong in `bench_*.mojo`; raw results are development evidence, not
permanent marketing claims.
