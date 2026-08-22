# Linear batch profile

This note records the optimization evidence for the sorted-query API. It is an
engineering result from one machine, not a portable performance guarantee.

## Environment and method

- Apple M4, arm64, macOS Darwin 25.5.0
- Mojo 1.0.0 (`ed45d567`), `--optimization-level 3`
- 65,536 deterministic irregular knots; both 65,536 dense queries and a
  128-query cluster over the final 128 knots
- caller-owned result buffer, identical scalar and sorted checksums
- 31 latency samples after two warmups; nearest-rank p50/p95
- a `keep` barrier inside every timed repetition, preventing overwritten batch
  calls from being collapsed
- `/usr/bin/sample`, one-millisecond interval, three seconds per compiled mode

The latency runner and profiler driver are deliberately separate. The former
produces independent samples; the latter repeats a steady-state workload long
enough for statistical call-stack collection. Query generation, interpolator
construction, and output allocation occur outside both timed/profiled loops.

## Result

| Query shape / API | Evaluations/sample | p50 | p95 |
| --- | ---: | ---: | ---: |
| dense / scalar binary search | 524,288 | 45.103 ms | 45.481 ms |
| dense / sorted seeded scan | 524,288 | 5.591 ms | 5.605 ms |
| sparse late / scalar binary search | 524,288 | 40.332 ms | 49.448 ms |
| sparse late / sorted seeded scan | 524,288 | 7.754 ms | 8.670 ms |

The sorted API was 8.07x/8.11x faster at dense p50/p95 and 5.20x/5.70x faster
for the sparse late cluster while producing the exact same `Float64` checksums.
Seeding the first interior interval with one binary search is what prevents a
sparse batch from walking tens of thousands of irrelevant leading knots.

In the compiled profiles, dense scalar evaluation placed 1,803 of 2,255
steady-state main-thread samples (80.0%) in per-query interval search and 327
(14.5%) in segment evaluation. Dense sorted evaluation placed 1,632 of 2,553
samples (63.9%) in segment evaluation, with the remaining 36.1% in the inlined
validation/traversal path. Sparse scalar showed the same search hotspot: 2,038
of 2,572 samples (79.2%). Sparse sorted placed 1,596 of 2,554 samples (62.5%) in
segment evaluation and only 39 (1.5%) in the one-per-batch binary-search seed;
the remaining 36.0% covered validation and monotone traversal. This is the
expected hotspot shift without the former sparse-query regression.

## SIMD acceptance decision

A guarded four-segment `Float64x4` candidate was tested against the scalar
oracle, including exact knots, duplicates, tails, extreme finite coordinates,
and non-finite rejection. It preserved results but increased the compiled
dense sorted workload from 5.85 seconds to 10.88 seconds (1.86x slower) on the
same machine. Lane assembly and eligibility/fallback checks outweighed
arithmetic savings, so the SIMD candidate was removed. The shipped
implementation keeps the faster scalar scan and all tail/extreme/non-finite
tests.
