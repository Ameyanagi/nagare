"""Long-running compiled workload for statistical CPU profiling.

This is deliberately separate from the latency benchmark: profiler sampling
needs seconds of steady-state work, while the benchmark needs independent
short samples. Build with optimization and invoke one dense or sparse mode.
"""

from nagare import LinearInterpolator
from std.benchmark import keep
from std.collections import List
from std.sys import argv

from linear_fixtures import make_linear_fixture


comptime KNOT_COUNT = 65_536
comptime SCALAR_PASSES = 4_096
comptime SORTED_PASSES = 32_768
comptime SPARSE_SCALAR_PASSES = 524_288
comptime SPARSE_SORTED_PASSES = 4_194_304


def _queries(knots: List[Float64]) -> List[Float64]:
    var queries = List[Float64](capacity=len(knots))
    for index in range(len(knots) - 1):
        var left = knots[index]
        queries.append(left + 0.375 * (knots[index + 1] - left))
    queries.append(knots[len(knots) - 1])
    return queries^


def _sparse_queries(knots: List[Float64]) -> List[Float64]:
    comptime query_count = 128
    var start = len(knots) - query_count
    var queries = List[Float64](capacity=query_count)
    for index in range(start, len(knots) - 1):
        var left = knots[index]
        queries.append(left + 0.375 * (knots[index + 1] - left))
    queries.append(knots[len(knots) - 1])
    return queries^


def _checksum(values: List[Float64]) -> Float64:
    var result = 0.0
    for value in values:
        result += value
    return result


def main() raises:
    var args = argv()
    if len(args) != 2:
        raise Error("usage: profile_linear <scalar|sorted|sparse-scalar|sparse-sorted>")
    var mode = String(args[1])
    if (
        mode != "scalar"
        and mode != "sorted"
        and mode != "sparse-scalar"
        and mode != "sparse-sorted"
    ):
        raise Error("usage: profile_linear <scalar|sorted|sparse-scalar|sparse-sorted>")

    var fixture = make_linear_fixture(KNOT_COUNT, True)
    var interpolator = LinearInterpolator(fixture.knots.copy(), fixture.values.copy())
    var sparse = mode.startswith("sparse-")
    var queries = _sparse_queries(fixture.knots) if sparse else _queries(fixture.knots)
    var results = List[Float64](length=len(queries), fill=0.0)
    var use_scalar = mode == "scalar" or mode == "sparse-scalar"
    var passes = SCALAR_PASSES
    if mode == "sorted":
        passes = SORTED_PASSES
    elif mode == "sparse-scalar":
        passes = SPARSE_SCALAR_PASSES
    elif mode == "sparse-sorted":
        passes = SPARSE_SORTED_PASSES
    for _ in range(passes):
        if use_scalar:
            interpolator.evaluate_into(queries, results)
        else:
            interpolator.evaluate_sorted_into(queries, results)
        keep(results)
    print("mode=", mode, " passes=", passes, " checksum=", _checksum(results), sep="")
