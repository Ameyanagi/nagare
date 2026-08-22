# Changelog

This project follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and uses semantic versioning after the first public release.

## [Unreleased]

### Added

- Initial experimental repository scaffold.
- A validated, right-biased binary-search interval locator.
- An owning linear interpolator with explicit error, clamp, and linear
  extrapolation policies.
- Direct-delta and scale-normalized affine evaluation paths that avoid
  premature overflow for finite in-domain data, preserve tiny-ordinate
  extrapolation when representable, and define signed-infinity extrapolation
  overflow.
- Construction-time validation with an explicit `validate()` checkpoint and
  trusted read-only operations, plus a bounded-result contract for
  ill-conditioned full-range intervals.
- Exact, independently calculated reference, and property/invariant tests for
  the first numerical vertical slice.
- An issue-sized v0.1 execution plan with numerical contracts and validation
  gates.
- A versioned 32-case linear construction/evaluation benchmark baseline with
  deterministic fixtures, host/compiler metadata, semantic checksums, reusable
  scalar/sorted batch cases, and no cross-project performance claim.
- Benchmark manifest v3 locks workload sizes, iterations, and exact checksums;
  run metadata identifies Git/lockfile state and fixed extreme cases use
  per-iteration `black_box` inputs, retained outputs, non-zero semantic
  sentinels, and nearest-rank p50/p95 timings from 31 non-zero-duration samples.
- Allocation-returning and caller-buffer sorted-query evaluation for linear
  interpolation. It preserves scalar extrapolation and floating-point
  semantics while replacing per-query binary search with one initial search and
  an `O(log n + s + m)` monotone scan over only the spanned segment range.
- A compiled dense/sparse scalar/sorted profiling driver and macOS `sample` /
  Linux `perf` workflow, with measurements and the rejected SIMD experiment
  documented.
- Allocating and caller-buffer batch first-derivative evaluation on every
  differentiable interpolator, plus scalar second derivatives for the Hermite
  family with explicit extrapolation semantics.
- Scalar and batch `__call__` evaluation sugar on every interpolator while
  retaining `evaluate` as the primary documented name.
- Read-only `knots()` and `values()` views on every interpolator, plus
  `slopes()` for user-supplied cubic Hermite slopes.

### Changed

- **Breaking:** Renamed `AkimaInterpolator` to `MakimaInterpolator` because the
  implementation matches scipy `Akima1DInterpolator(method="makima")`, not its
  classic `method="akima"` default. No compatibility alias is provided before
  the first public release.
