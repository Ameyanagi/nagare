"""Reproducible construction and scalar-evaluation benchmark baseline."""

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
comptime SAMPLE_COUNT = 7


def _print_case(
    identity: String,
    knot_count: Int,
    iterations: Int,
    expected_checksum: Float64,
    manifest_only: Bool,
    best_elapsed_ns: Int = 0,
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
        " best_elapsed_ns=",
        best_elapsed_ns,
        " statistic=min_of_",
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

    var best_elapsed_ns = 0
    for sample in range(SAMPLE_COUNT):
        var started = perf_counter_ns()
        for _ in range(iterations):
            last = LinearInterpolator(fixture.knots.copy(), fixture.values.copy())
            keep(last)
        var elapsed_ns = perf_counter_ns() - started
        if sample == 0 or elapsed_ns < best_elapsed_ns:
            best_elapsed_ns = elapsed_ns
    _print_case(
        identity,
        len(fixture.knots),
        iterations,
        Float64(len(fixture.knots)),
        False,
        best_elapsed_ns,
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

    var best_elapsed_ns = 0
    var best_checksum = 0.0
    for sample in range(SAMPLE_COUNT):
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
        var elapsed_ns = perf_counter_ns() - started
        keep(checksum)
        if checksum != expected_checksum:
            raise Error("evaluation checksum changed during measurement")
        if sample == 0 or elapsed_ns < best_elapsed_ns:
            best_elapsed_ns = elapsed_ns
            best_checksum = checksum

    _print_case(
        identity,
        len(fixture.knots),
        iterations,
        expected_checksum,
        False,
        best_elapsed_ns,
        best_checksum,
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

    var best_elapsed_ns = 0
    var best_checksum = 0.0
    for sample in range(SAMPLE_COUNT):
        checksum = 0.0
        var started = perf_counter_ns()
        for _ in range(iterations):
            var value = interpolator.evaluate(black_box(query))
            keep(value)
            checksum += value / checksum_scale
        var elapsed_ns = perf_counter_ns() - started
        keep(checksum)
        if checksum != expected_checksum:
            raise Error("extreme checksum changed during measurement")
        if sample == 0 or elapsed_ns < best_elapsed_ns:
            best_elapsed_ns = elapsed_ns
            best_checksum = checksum

    _print_case(
        identity,
        len(fixture.knots),
        iterations,
        expected_checksum,
        False,
        best_elapsed_ns,
        best_checksum,
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

    print("schema=nagare-linear-benchmark-v2")
    print("warmup_rounds=", WARMUP_ROUNDS, sep="")
    print("samples=", SAMPLE_COUNT, sep="")
    print("statistic=minimum elapsed nanoseconds across samples")

    _run_size("small", 8, 2_000, 10_000, manifest_only)
    _run_size("medium", 1_024, 64, 128, manifest_only)
    _run_size("large", 65_536, 2, 4, manifest_only)
    _measure_fixed_query(
        "evaluation.linear.extreme.tiny_ordinates",
        make_tiny_ordinate_fixture(),
        ExtrapolationPolicy.LINEAR,
        1.0,
        1.0,
        1_024,
        manifest_only,
    )
    _measure_fixed_query(
        "evaluation.error.extreme.full_range",
        make_full_range_fixture(),
        ExtrapolationPolicy.ERROR,
        Float64.MAX_FINITE / 4_096.0,
        Float64.MAX_FINITE,
        1_024,
        manifest_only,
    )
