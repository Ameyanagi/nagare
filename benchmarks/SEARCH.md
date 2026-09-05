# Repeated interval lookup

Run from the source checkout:

```sh
pixi run --locked bench-search
```

The task builds `benchmarks/bench_search.mojo` with Mojo 1.0.0 at optimization
level 3, then measures the compiled executable. Compilation, allocation,
`KnotIndex` construction, and its explicit input copy are excluded from query
timing. Construction performs one O(n) scan; the benchmark measures subsequent
queries, where the one-shot API repeats that scan and `KnotIndex` does not.

Each case uses irregular exact dyadic knots `k + (k % 7) / 16`. Queries select
interval `(q * 8191 + 17) % (n - 1)` and use its midpoint. Both paths receive
identical unsorted queries. A black-box query and retained result/checksum keep
the workload observable; each measured batch must match the independently
computed sum of interval indices. Two warmups precede 31 batches. Reported
nanoseconds per query divide the median batch time by its query count.

## Initial measurement

2026-09-05, Apple M4, macOS 26.5.1 (25F80), Mojo 1.0.0. Other ecosystem builds
were running on the same host, so these are representative observations, not a
portable latency guarantee or a CI performance threshold.

| Knots | Queries per batch | One-shot ns/query | Reusable ns/query | Speedup |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 4,096 | 40.28 | 13.18 | 3.06× |
| 1,024 | 1,024 | 3,758.79 | 52.73 | 71.28× |
| 65,536 | 256 | 254,281.25 | 128.91 | 1,972.61× |

The improvement comes from amortizing knot validation, changing Q repeated
queries from O(Q n) validation plus search to O(n + Q log n). No SIMD kernel or
changed interval semantics is involved. The existing `locate_interval` keeps
its validated one-shot contract, which remains useful when the input table is
not retained across calls.
