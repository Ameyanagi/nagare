"""Reproducible construction, scalar, and sorted-batch benchmark baseline."""

from std.benchmark import black_box, keep
from std.sys import argv
from std.time import perf_counter_ns

from linear_fixtures import (
    LinearBenchmarkFixture,
    make_full_range_fixture,
    make_linear_fixture,
    make_tiny_ordinate_fixture,
)
from nagare import ExtrapolationPolicy, LinearInterpolator

comptime WARMUP_ROUNDS = 2
comptime SAMPLE_COUNT = 31


def _sort_timings(mut values: List[Int]):
    # The fixed 31-sample set makes this benchmark-local insertion sort tiny.
    for index in range(1, len(values)):
        var value = values[index]
        var position = index
        while position > 0 and values[position - 1] > value:
            values[position] = values[position - 1]
            position -= 1
        values[position] = value


def _print_case(
    identity: String,
    knot_count: Int,
    iterations: Int,
    expected_checksum: Float64,
    manifest_only: Bool,
    p50_elapsed_ns: Int = 0,
    p95_elapsed_ns: Int = 0,
    checksum: Float64 = 0.0,
):
    if manifest_only:
        print(
            "case=",
            identity,
            " knots=",
            knot_count,
            " iterations=",
            iterations,
            " checksum_contract=exact_float64 expected_checksum=",
            expected_checksum,
            sep="",
        )
        return
    print(
        "case=",
        identity,
        " knots=",
        knot_count,
        " iterations=",
        iterations,
        " checksum_contract=exact_float64 expected_checksum=",
        expected_checksum,
        " p50_elapsed_ns=",
        p50_elapsed_ns,
        " p95_elapsed_ns=",
        p95_elapsed_ns,
        " statistic=nearest_rank_p50_p95_of_",
        SAMPLE_COUNT,
        " checksum=",
        checksum,
        sep="",
    )


def _measure_construction(
    identity: String,
    fixture: LinearBenchmarkFixture,
    iterations: Int,
    manifest_only: Bool,
) raises:
    if manifest_only:
        _print_case(
            identity,
            len(fixture.knots),
            iterations,
            Float64(len(fixture.knots)),
            True,
        )
        return

    var last = LinearInterpolator(fixture.knots.copy(), fixture.values.copy())
    for _ in range(WARMUP_ROUNDS):
        for _ in range(iterations):
            last = LinearInterpolator(fixture.knots.copy(), fixture.values.copy())
            keep(last)

    var timings = List[Int](capacity=SAMPLE_COUNT)
    for _ in range(SAMPLE_COUNT):
        var started = perf_counter_ns()
        for _ in range(iterations):
            last = LinearInterpolator(fixture.knots.copy(), fixture.values.copy())
            keep(last)
        timings.append(perf_counter_ns() - started)
    _sort_timings(timings)
    _print_case(
        identity,
        len(fixture.knots),
        iterations,
        Float64(len(fixture.knots)),
        False,
        timings[SAMPLE_COUNT // 2],
        timings[(SAMPLE_COUNT * 95 + 99) // 100 - 1],
        Float64(last.knot_count()),
    )


def _in_domain_query(fixture: LinearBenchmarkFixture, iteration: Int) -> Float64:
    var index = iteration % (len(fixture.knots) - 1)
    var left = fixture.knots[index]
    return left + 0.375 * (fixture.knots[index + 1] - left)


def _out_of_domain_query(fixture: LinearBenchmarkFixture, iteration: Int) -> Float64:
    if iteration % 2 == 0:
        return fixture.knots[0] - 0.5 * (fixture.knots[1] - fixture.knots[0])
    var final_index = len(fixture.knots) - 1
    return fixture.knots[final_index] + 0.5 * (
        fixture.knots[final_index] - fixture.knots[final_index - 1]
    )


def _measure_evaluation(
    identity: String,
    fixture: LinearBenchmarkFixture,
    policy: ExtrapolationPolicy,
    iterations: Int,
    use_out_of_domain_queries: Bool,
    manifest_only: Bool,
) raises:
    var interpolator = LinearInterpolator(
        fixture.knots.copy(),
        fixture.values.copy(),
        extrapolation=policy,
    )
    var expected_checksum = 0.0
    for iteration in range(iterations):
        var query = _out_of_domain_query(
            fixture, iteration
        ) if use_out_of_domain_queries else _in_domain_query(fixture, iteration)
        expected_checksum += interpolator.evaluate(query)
    if manifest_only:
        _print_case(
            identity,
            len(fixture.knots),
            iterations,
            expected_checksum,
            True,
        )
        return

    var checksum: Float64
    for _ in range(WARMUP_ROUNDS):
        checksum = 0.0
        for iteration in range(iterations):
            var query = _out_of_domain_query(
                fixture, iteration
            ) if use_out_of_domain_queries else _in_domain_query(fixture, iteration)
            keep(query)
            var value = interpolator.evaluate(query)
            keep(value)
            checksum += value
        keep(checksum)
        if checksum != expected_checksum:
            raise Error("evaluation checksum changed during warmup")

    var timings = List[Int](capacity=SAMPLE_COUNT)
    var measured_checksum = 0.0
    for _ in range(SAMPLE_COUNT):
        checksum = 0.0
        var started = perf_counter_ns()
        for iteration in range(iterations):
            var query = _out_of_domain_query(
                fixture, iteration
            ) if use_out_of_domain_queries else _in_domain_query(fixture, iteration)
            keep(query)
            var value = interpolator.evaluate(query)
            keep(value)
            checksum += value
        timings.append(perf_counter_ns() - started)
        keep(checksum)
        if checksum != expected_checksum:
            raise Error("evaluation checksum changed during measurement")
        measured_checksum = checksum

    _sort_timings(timings)
    _print_case(
        identity,
        len(fixture.knots),
        iterations,
        expected_checksum,
        False,
        timings[SAMPLE_COUNT // 2],
        timings[(SAMPLE_COUNT * 95 + 99) // 100 - 1],
        measured_checksum,
    )


def _measure_fixed_query(
    identity: String,
    fixture: LinearBenchmarkFixture,
    policy: ExtrapolationPolicy,
    query: Float64,
    checksum_scale: Float64,
    iterations: Int,
    manifest_only: Bool,
) raises:
    var interpolator = LinearInterpolator(
        fixture.knots.copy(),
        fixture.values.copy(),
        extrapolation=policy,
    )
    var expected_checksum = 0.0
    for _ in range(iterations):
        expected_checksum += interpolator.evaluate(query) / checksum_scale
    if manifest_only:
        _print_case(
            identity,
            len(fixture.knots),
            iterations,
            expected_checksum,
            True,
        )
        return

    var checksum: Float64
    for _ in range(WARMUP_ROUNDS):
        checksum = 0.0
        for _ in range(iterations):
            var value = interpolator.evaluate(black_box(query))
            keep(value)
            checksum += value / checksum_scale
        keep(checksum)
        if checksum != expected_checksum:
            raise Error("extreme checksum changed during warmup")

    var timings = List[Int](capacity=SAMPLE_COUNT)
    var measured_checksum = 0.0
    for _ in range(SAMPLE_COUNT):
        checksum = 0.0
        var started = perf_counter_ns()
        for _ in range(iterations):
            var value = interpolator.evaluate(black_box(query))
            keep(value)
            checksum += value / checksum_scale
        timings.append(perf_counter_ns() - started)
        keep(checksum)
        if checksum != expected_checksum:
            raise Error("extreme checksum changed during measurement")
        measured_checksum = checksum

    _sort_timings(timings)
    _print_case(
        identity,
        len(fixture.knots),
        iterations,
        expected_checksum,
        False,
        timings[SAMPLE_COUNT // 2],
        timings[(SAMPLE_COUNT * 95 + 99) // 100 - 1],
        measured_checksum,
    )


def _sorted_queries(fixture: LinearBenchmarkFixture) -> List[Float64]:
    """Return one interior query per segment followed by the final knot."""
    var queries = List[Float64](capacity=len(fixture.knots))
    for index in range(len(fixture.knots) - 1):
        var left = fixture.knots[index]
        queries.append(left + 0.375 * (fixture.knots[index + 1] - left))
    queries.append(fixture.knots[len(fixture.knots) - 1])
    return queries^


def _sparse_late_queries(fixture: LinearBenchmarkFixture) -> List[Float64]:
    """Return 128 sorted queries clustered at the end of a large table."""
    comptime query_count = 128
    var start = len(fixture.knots) - query_count
    var queries = List[Float64](capacity=query_count)
    for index in range(start, len(fixture.knots) - 1):
        var left = fixture.knots[index]
        queries.append(left + 0.375 * (fixture.knots[index + 1] - left))
    queries.append(fixture.knots[len(fixture.knots) - 1])
    return queries^


def _result_checksum(results: List[Float64]) -> Float64:
    var checksum = 0.0
    for value in results:
        checksum += value
    return checksum


def _measure_reusable_batch(
    identity: String,
    fixture: LinearBenchmarkFixture,
    repetitions: Int,
    use_sorted_api: Bool,
    manifest_only: Bool,
    sparse_late_cluster: Bool = False,
) raises:
    var interpolator = LinearInterpolator(fixture.knots.copy(), fixture.values.copy())
    var queries = _sparse_late_queries(
        fixture
    ) if sparse_late_cluster else _sorted_queries(fixture)
    var results = List[Float64](length=len(queries), fill=0.0)
    interpolator.evaluate_into(queries, results)
    var expected_checksum = _result_checksum(results)
    var evaluations = repetitions * len(queries)
    if manifest_only:
        _print_case(
            identity,
            len(fixture.knots),
            evaluations,
            expected_checksum,
            True,
        )
        return

    for _ in range(WARMUP_ROUNDS):
        for _ in range(repetitions):
            if use_sorted_api:
                interpolator.evaluate_sorted_into(queries, results)
            else:
                interpolator.evaluate_into(queries, results)
            keep(results)
        if _result_checksum(results) != expected_checksum:
            raise Error("batch checksum changed during warmup")

    var timings = List[Int](capacity=SAMPLE_COUNT)
    for _ in range(SAMPLE_COUNT):
        var started = perf_counter_ns()
        for _ in range(repetitions):
            if use_sorted_api:
                interpolator.evaluate_sorted_into(queries, results)
            else:
                interpolator.evaluate_into(queries, results)
            keep(results)
        timings.append(perf_counter_ns() - started)
        if _result_checksum(results) != expected_checksum:
            raise Error("batch checksum changed during measurement")

    _sort_timings(timings)
    _print_case(
        identity,
        len(fixture.knots),
        evaluations,
        expected_checksum,
        False,
        timings[SAMPLE_COUNT // 2],
        timings[(SAMPLE_COUNT * 95 + 99) // 100 - 1],
        expected_checksum,
    )


def _run_size(
    size_name: String,
    count: Int,
    construction_iterations: Int,
    evaluation_iterations: Int,
    manifest_only: Bool,
) raises:
    for irregular in [False, True]:
        var shape = String("irregular" if irregular else "uniform")
        var fixture = make_linear_fixture(count, irregular)
        var prefix = String(size_name, ".", shape)
        _measure_construction(
            String("construction.", prefix),
            fixture,
            construction_iterations,
            manifest_only,
        )
        _measure_evaluation(
            String("evaluation.error.", prefix),
            fixture,
            ExtrapolationPolicy.ERROR,
            evaluation_iterations,
            False,
            manifest_only,
        )
        _measure_evaluation(
            String("evaluation.clamp.", prefix),
            fixture,
            ExtrapolationPolicy.CLAMP,
            evaluation_iterations,
            True,
            manifest_only,
        )
        _measure_evaluation(
            String("evaluation.linear.", prefix),
            fixture,
            ExtrapolationPolicy.LINEAR,
            evaluation_iterations,
            True,
            manifest_only,
        )


def main() raises:
    var raw_args = argv()
    var manifest_only = len(raw_args) == 2 and String(raw_args[1]) == "--manifest"
    if len(raw_args) > 1 and not manifest_only:
        raise Error("usage: bench_linear.mojo [--manifest]")

    print("schema=nagare-linear-benchmark-v3")
    print("warmup_rounds=", WARMUP_ROUNDS, sep="")
    print("samples=", SAMPLE_COUNT, sep="")
    print("statistic=nearest-rank p50/p95 elapsed nanoseconds across samples")

    _run_size("small", 8, 50_000, 500_000, manifest_only)
    _run_size("medium", 1_024, 4_096, 500_000, manifest_only)
    _run_size("large", 65_536, 64, 500_000, manifest_only)
    var medium_batch = make_linear_fixture(1_024, True)
    _measure_reusable_batch(
        "batch.scalar.medium.irregular", medium_batch, 256, False, manifest_only
    )
    _measure_reusable_batch(
        "batch.sorted.medium.irregular", medium_batch, 256, True, manifest_only
    )
    var large_batch = make_linear_fixture(65_536, True)
    _measure_reusable_batch(
        "batch.scalar.large.irregular", large_batch, 8, False, manifest_only
    )
    _measure_reusable_batch(
        "batch.sorted.large.irregular", large_batch, 8, True, manifest_only
    )
    _measure_reusable_batch(
        "batch.scalar.large.irregular.sparse_late",
        large_batch,
        4_096,
        False,
        manifest_only,
        True,
    )
    _measure_reusable_batch(
        "batch.sorted.large.irregular.sparse_late",
        large_batch,
        4_096,
        True,
        manifest_only,
        True,
    )
    _measure_fixed_query(
        "evaluation.linear.extreme.tiny_ordinates",
        make_tiny_ordinate_fixture(),
        ExtrapolationPolicy.LINEAR,
        1.0,
        1.0,
        1_000_000,
        manifest_only,
    )
    _measure_fixed_query(
        "evaluation.error.extreme.full_range",
        make_full_range_fixture(),
        ExtrapolationPolicy.ERROR,
        Float64.MAX_FINITE / 4_096.0,
        Float64.MAX_FINITE,
        1_000_000,
        manifest_only,
    )
