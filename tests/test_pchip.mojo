from nagare import (
    CubicSplineInterpolator,
    ExtrapolationPolicy,
    PchipInterpolator,
)
from std.collections import List
from std.testing import TestSuite, assert_equal, assert_raises, assert_true


def assert_close(
    actual: Float64,
    expected: Float64,
    atol: Float64 = 1e-9,
    rtol: Float64 = 1e-12,
) raises:
    var tolerance = max(atol, rtol * abs(expected))
    assert_true(
        abs(actual - expected) <= tolerance,
        msg=String("expected ", expected, ", got ", actual),
    )


def reference_interpolator(
    policy: ExtrapolationPolicy = ExtrapolationPolicy.ERROR,
) raises -> PchipInterpolator:
    return PchipInterpolator(
        [0.0, 1.0, 2.5, 3.0, 4.5, 6.0],
        [0.0, 2.0, 1.0, 3.5, 3.0, 4.0],
        extrapolation=policy,
    )


def metadata_matches_reference(interpolator: PchipInterpolator) -> Bool:
    """Exercise every metadata accessor from a non-raising function."""
    return (
        interpolator.knot_count() == 6
        and interpolator.domain_start() == 0.0
        and interpolator.domain_end() == 6.0
    )


def test_scipy_reference_values_and_derivatives() raises:
    var interpolator = reference_interpolator()
    var queries: List[Float64] = [
        0.25,
        0.5,
        1.2,
        2.6,
        2.9,
        3.7,
        4.4,
        5.1,
        5.9,
    ]

    # Reference values from scipy 1.18.0 (numpy 2.5.2): PchipInterpolator(X, Y)(Q)
    var expected_values: List[Float64] = [
        0.74375,
        1.3833333333333333,
        1.9514074074074075,
        1.2600000000000005,
        3.239999999999999,
        3.274962962962963,
        3.00637037037037,
        3.1839999999999997,
        3.8856296296296304,
    ]
    # Reference values from scipy 1.18.0 (numpy 2.5.2):
    # PchipInterpolator(X, Y).derivative()(Q)
    var expected_derivatives: List[Float64] = [
        2.825,
        2.2333333333333334,
        -0.46222222222222215,
        4.800000000000003,
        4.800000000000004,
        -0.4977777777777778,
        -0.12444444444444414,
        0.5866666666666663,
        1.1200000000000006,
    ]

    for index in range(len(queries)):
        assert_close(interpolator.evaluate(queries[index]), expected_values[index])
        assert_close(
            interpolator.derivative(queries[index]), expected_derivatives[index]
        )


def test_monotone_data_never_overshoots_but_natural_cubic_does() raises:
    var knots: List[Float64] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    var values: List[Float64] = [0.0, 0.1, 0.2, 5.0, 9.9, 10.0]
    var pchip = PchipInterpolator(knots.copy(), values.copy())
    var cubic = CubicSplineInterpolator(knots^, values^)

    var previous = pchip.evaluate(0.0)
    for step in range(251):
        var value = pchip.evaluate(Float64(step) * 0.02)
        assert_true(value >= 0.0 and value <= 10.0)
        assert_true(value >= previous)
        previous = value

    assert_true(cubic.evaluate(1.5) < 0.0)

    var queries: List[Float64] = [2.25, 2.5, 2.75, 3.25, 3.5, 3.75]
    # Reference values from scipy 1.18.0 (numpy 2.5.2):
    # PchipInterpolator(XM, YM)(QM)
    var expected: List[Float64] = [
        0.7502314327792973,
        2.0183042289080584,
        3.5772249105827907,
        6.438396262886598,
        8.03168556701031,
        9.334132087628866,
    ]
    for index in range(len(queries)):
        assert_close(pchip.evaluate(queries[index]), expected[index])


def test_metadata_exact_knots_and_two_knot_special_case() raises:
    var interpolator = reference_interpolator()
    assert_true(metadata_matches_reference(interpolator))
    interpolator.validate()

    var knots: List[Float64] = [0.0, 1.0, 2.5, 3.0, 4.5, 6.0]
    var values: List[Float64] = [0.0, 2.0, 1.0, 3.5, 3.0, 4.0]
    for index in range(len(knots)):
        assert_equal(interpolator.evaluate(knots[index]), values[index])

    var two_knots = PchipInterpolator([-2.0, 3.0], [4.0, -6.0])
    assert_close(two_knots.derivative(-2.0), -2.0)
    assert_close(two_knots.derivative(0.0), -2.0)
    assert_close(two_knots.derivative(3.0), -2.0)


def test_constructor_uses_shared_table_validation() raises:
    with assert_raises(contains="at least two"):
        _ = PchipInterpolator([0.0], [1.0])
    with assert_raises(contains="len(knots) = 2, len(values) = 1"):
        _ = PchipInterpolator([0.0, 1.0], [1.0])
    with assert_raises(contains="knots[1] = 0.0 <= knots[0] = 0.0"):
        _ = PchipInterpolator([0.0, 0.0], [1.0, 2.0])
    with assert_raises(contains="values[1] is nan"):
        _ = PchipInterpolator([0.0, 1.0], [1.0, Float64("nan")])


def test_explicit_validate_rechecks_derived_slopes() raises:
    var valid = reference_interpolator()
    valid.validate()

    var non_finite = reference_interpolator()
    non_finite._engine._slopes[2] = Float64("nan")
    with assert_raises(contains="slopes[2] is nan"):
        non_finite.validate()


def test_error_clamp_linear_and_fill_extrapolation() raises:
    var error = reference_interpolator()
    with assert_raises(contains="query -1.0 is outside the knot domain [0.0, 6.0]"):
        _ = error.evaluate(-1.0)
    with assert_raises(contains="outside the knot domain"):
        _ = error.derivative(7.0)

    var clamp = reference_interpolator(ExtrapolationPolicy.CLAMP)
    assert_equal(clamp.evaluate(-1.0), 0.0)
    assert_equal(clamp.evaluate(7.0), 4.0)
    assert_equal(clamp.derivative(-1.0), 0.0)
    assert_equal(clamp.derivative(7.0), 0.0)

    # scipy-compatible endpoint estimates are 46/15 on the left and 7/6 on
    # the right for the shared uneven table.
    var linear = reference_interpolator(ExtrapolationPolicy.LINEAR)
    assert_close(linear.evaluate(-1.0), -46.0 / 15.0)
    assert_close(linear.evaluate(7.0), 4.0 + 7.0 / 6.0)
    assert_close(linear.derivative(-1.0), 46.0 / 15.0)
    assert_close(linear.derivative(7.0), 7.0 / 6.0)

    var default_fill = reference_interpolator(ExtrapolationPolicy.FILL)
    assert_true(default_fill.evaluate(-1.0) != default_fill.evaluate(-1.0))
    assert_true(default_fill.derivative(7.0) != default_fill.derivative(7.0))

    var custom_fill = reference_interpolator(ExtrapolationPolicy.fill(-3.25))
    assert_equal(custom_fill.evaluate(-1.0), -3.25)
    assert_equal(custom_fill.evaluate(7.0), -3.25)
    assert_equal(custom_fill.derivative(-1.0), -3.25)
    assert_equal(custom_fill.derivative(7.0), -3.25)


def test_non_finite_queries_always_raise() raises:
    var interpolator = reference_interpolator(ExtrapolationPolicy.FILL)
    with assert_raises(contains="query must be finite: received nan"):
        _ = interpolator.evaluate(Float64("nan"))
    with assert_raises(contains="query must be finite"):
        _ = interpolator.derivative(Float64("-inf"))


def test_span_wrapper_and_evaluate_into_agree() raises:
    var queries: List[Float64] = [-1.0, 0.0, 0.5, 2.6, 4.4, 6.0, 7.0]
    var interpolator = reference_interpolator(ExtrapolationPolicy.fill(-3.25))
    var results = List[Float64](length=len(queries), fill=99.0)
    interpolator.evaluate_into(queries, results)
    var allocated = interpolator.evaluate(queries)

    for index in range(len(queries)):
        assert_equal(results[index], allocated[index])
        assert_equal(results[index], interpolator.evaluate(queries[index]))

    var short_results = List[Float64](length=1, fill=0.0)
    with assert_raises(contains="len(queries) = 7, len(results) = 1"):
        interpolator.evaluate_into(queries, short_results)


def test_closed_form_integrals_and_equality_writable_surface() raises:
    var first = reference_interpolator()
    var second = reference_interpolator()
    var clamp = reference_interpolator(ExtrapolationPolicy.CLAMP)
    # scipy 1.18.0 / numpy 2.5.2: PchipInterpolator(X, Y).integrate(...)
    assert_close(first.integrate(0.0, 6.0), 14.536805555555556)
    assert_close(first.integrate(1.2, 4.8), 8.758168518518518)
    assert_true(first == second)
    assert_true(first != clamp)
    assert_true(
        String(first).startswith(
            "PchipInterpolator(6 knots on [0.0, 6.0], extrapolation=ERROR)"
        )
    )


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
