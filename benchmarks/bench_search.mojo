"""Compare validated one-shot lookup with an owned reusable knot index."""

from nagare import KnotIndex, locate_interval
from std.benchmark import black_box, keep
from std.collections import List
from std.time import perf_counter_ns

comptime SAMPLE_COUNT = 31


def _median(mut timings: List[Int]) -> Int:
    for index in range(1, len(timings)):
        var value = timings[index]
        var position = index
        while position > 0 and timings[position - 1] > value:
            timings[position] = timings[position - 1]
            position -= 1
        timings[position] = value
    return timings[SAMPLE_COUNT // 2]


def _measure[
    reusable: Bool
](
    knots: List[Float64], index: KnotIndex, queries: List[Float64], expected: Int
) raises -> Float64:
    var timings = List[Int](capacity=SAMPLE_COUNT)
    # Two warmups precede the 31 independently timed batches.
    for sample in range(SAMPLE_COUNT + 2):
        var checksum = 0
        var started = perf_counter_ns()
        for query in queries:
            var value: Int
            comptime if reusable:
                value = index.locate(black_box(query))
            else:
                value = locate_interval(knots, black_box(query))
            keep(value)
            checksum += value
        var elapsed = perf_counter_ns() - started
        keep(checksum)
        if checksum != expected:
            raise Error(
                String("search checksum must be ", expected, "; got ", checksum)
            )
        if sample >= 2:
            timings.append(elapsed)
    var ns_per_query = Float64(_median(timings)) / Float64(len(queries))
    print(
        "mode=",
        "reusable" if reusable else "one-shot",
        " knots=",
        len(knots),
        " queries=",
        len(queries),
        " p50_ns_per_query=",
        ns_per_query,
        " checksum=",
        expected,
        sep="",
    )
    return ns_per_query


def _case(count: Int, query_count: Int) raises:
    var knots = List[Float64](capacity=count)
    for k in range(count):
        # Dyadic irregular increments keep fixture comparisons exact.
        knots.append(Float64(k) + Float64(k % 7) / 16.0)
    var queries = List[Float64](capacity=query_count)
    var expected = 0
    for q in range(query_count):
        var interval = (q * 8191 + 17) % (count - 1)
        queries.append((knots[interval] + knots[interval + 1]) / 2.0)
        expected += interval
    # Construction and the explicitly requested copy are outside query timing.
    var index = KnotIndex(knots.copy())
    var one_shot = _measure[False](knots, index, queries, expected)
    var reusable = _measure[True](knots, index, queries, expected)
    print("knots=", count, " speedup=", one_shot / reusable, sep="")


def main() raises:
    print("schema=nagare-interval-search-v1 samples=31 warmups=2 statistic=p50")
    _case(8, 4096)
    _case(1024, 1024)
    _case(65536, 256)
