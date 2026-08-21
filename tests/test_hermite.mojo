from nagare import CubicHermiteInterpolator, ExtrapolationPolicy
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
) raises -> CubicHermiteInterpolator:
    return CubicHermiteInterpolator(
        [0.0, 1.0, 2.5, 3.0, 4.5, 6.0],
        [0.0, 2.0, 1.0, 3.5, 3.0, 4.0],
        [1.0, -0.5, 2.0, 0.0, 1.5, -1.0],
        extrapolation=policy,
    )


def metadata_matches_reference(interpolator: CubicHermiteInterpolator) -> Bool:
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

    # Reference values from scipy 1.18.0 (numpy 2.5.2):
    # CubicHermiteSpline(X, Y, S)(Q)
    var expected_values: List[Float64] = [
        0.4765625,
        1.1875,
        1.8300740740740742,
        1.3880000000000006,
        3.272,
        3.013629629629629,
        2.875703703703704,
        3.8199999999999994,
        4.083703703703703,
    ]
    # Reference values from scipy 1.18.0 (numpy 2.5.2):
    # CubicHermiteSpline(X, Y, S).derivative()(Q)
    var expected_derivatives: List[Float64] = [
        2.59375,
        2.875,
        -1.1488888888888888,
        5.440000000000002,
        4.240000000000002,
        -0.9177777777777776,
        0.9955555555555575,
        1.1000000000000003,
        -0.6777777777777787,
    ]

    for index in range(len(queries)):
        assert_close(interpolator.evaluate(queries[index]), expected_values[index])
        assert_close(
            interpolator.derivative(queries[index]), expected_derivatives[index]
        )


def test_metadata_exact_knots_and_supplied_derivatives() raises:
    var interpolator = reference_interpolator()
    assert_true(metadata_matches_reference(interpolator))
    interpolator.validate()

    var knots: List[Float64] = [0.0, 1.0, 2.5, 3.0, 4.5, 6.0]
    var values: List[Float64] = [0.0, 2.0, 1.0, 3.5, 3.0, 4.0]
    var slopes: List[Float64] = [1.0, -0.5, 2.0, 0.0, 1.5, -1.0]
    for index in range(len(knots)):
        assert_equal(interpolator.evaluate(knots[index]), values[index])
        assert_close(interpolator.derivative(knots[index]), slopes[index], atol=1e-12)


def test_constructor_rejects_invalid_slope_tables_with_details() raises:
    with assert_raises(contains="at least two"):
        _ = CubicHermiteInterpolator([0.0], [1.0], [0.0])
    with assert_raises(contains="len(knots) = 3, len(values) = 3, len(slopes) = 2"):
        _ = CubicHermiteInterpolator(
            [0.0, 1.0, 2.0],
            [1.0, 2.0, 3.0],
            [0.0, 1.0],
        )
    with assert_raises(contains="slopes[1] is nan"):
        _ = CubicHermiteInterpolator(
            [0.0, 1.0],
            [1.0, 2.0],
            [0.0, Float64("nan")],
        )
    with assert_raises(contains="slopes[0] is inf"):
        _ = CubicHermiteInterpolator(
            [0.0, 1.0],
            [1.0, 2.0],
            [Float64("inf"), 0.0],
        )


def test_explicit_validate_rechecks_slopes() raises:
    var valid = reference_interpolator()
    valid.validate()

    var mismatched = reference_interpolator()
    mismatched._slopes = [0.0]
    with assert_raises(contains="len(slopes) = 1"):
        mismatched.validate()

    var non_finite = reference_interpolator()
    non_finite._slopes[2] = Float64("nan")
    with assert_raises(contains="slopes[2] is nan"):
        non_finite.validate()


def test_error_clamp_linear_and_fill_extrapolation() raises:
    var error = CubicHermiteInterpolator(
        [0.0, 1.0, 2.0],
        [1.0, 2.0, 0.0],
        [0.5, -1.0, 2.0],
    )
    with assert_raises(contains="query -0.25 is outside the knot domain [0.0, 2.0]"):
        _ = error.evaluate(-0.25)
    with assert_raises(contains="outside the knot domain"):
        _ = error.derivative(2.25)

    var clamp = CubicHermiteInterpolator(
        [0.0, 1.0, 2.0],
        [1.0, 2.0, 0.0],
        [0.5, -1.0, 2.0],
        extrapolation=ExtrapolationPolicy.CLAMP,
    )
    assert_equal(clamp.evaluate(-2.0), 1.0)
    assert_equal(clamp.evaluate(3.0), 0.0)
    assert_equal(clamp.derivative(-2.0), 0.0)
    assert_equal(clamp.derivative(3.0), 0.0)

    var linear = CubicHermiteInterpolator(
        [0.0, 1.0, 2.0],
        [1.0, 2.0, 0.0],
        [0.5, -1.0, 2.0],
        extrapolation=ExtrapolationPolicy.LINEAR,
    )
    assert_close(linear.evaluate(-2.0), 0.0)
    assert_close(linear.evaluate(3.0), 2.0)
    assert_equal(linear.derivative(-2.0), 0.5)
    assert_equal(linear.derivative(3.0), 2.0)

    var default_fill = reference_interpolator(ExtrapolationPolicy.FILL)
    assert_true(default_fill.evaluate(-1.0) != default_fill.evaluate(-1.0))
    assert_true(default_fill.derivative(7.0) != default_fill.derivative(7.0))

    var custom_fill = reference_interpolator(ExtrapolationPolicy.fill(-7.5))
    assert_equal(custom_fill.evaluate(-1.0), -7.5)
    assert_equal(custom_fill.evaluate(7.0), -7.5)
    assert_equal(custom_fill.derivative(-1.0), -7.5)
    assert_equal(custom_fill.derivative(7.0), -7.5)


def test_non_finite_queries_always_raise() raises:
    var interpolator = reference_interpolator(ExtrapolationPolicy.FILL)
    with assert_raises(contains="query must be finite"):
        _ = interpolator.evaluate(Float64("nan"))
    with assert_raises(contains="query must be finite"):
        _ = interpolator.derivative(Float64("inf"))


def test_span_wrapper_and_evaluate_into_agree() raises:
    var queries: List[Float64] = [-1.0, 0.0, 0.25, 2.5, 5.9, 6.0, 7.0]
    var interpolator = reference_interpolator(ExtrapolationPolicy.fill(-7.5))
    var results = List[Float64](length=len(queries), fill=99.0)
    interpolator.evaluate_into(queries, results)
    var allocated = interpolator.evaluate(queries)

    assert_equal(len(allocated), len(results))
    for index in range(len(results)):
        assert_equal(results[index], allocated[index])
        assert_equal(results[index], interpolator.evaluate(queries[index]))

    var short_results = List[Float64](length=2, fill=0.0)
    with assert_raises(contains="len(queries) = 7, len(results) = 2"):
        interpolator.evaluate_into(queries, short_results)


def test_closed_form_integral_and_equality_writable_surface() raises:
    var first = reference_interpolator()
    var second = reference_interpolator()
    var clamp = reference_interpolator(ExtrapolationPolicy.CLAMP)
    var different_slopes = CubicHermiteInterpolator(
        [0.0, 1.0, 2.5, 3.0, 4.5, 6.0],
        [0.0, 2.0, 1.0, 3.5, 3.0, 4.0],
        [1.0, -0.5, 2.0, 0.0, 1.5, 0.0],
    )
    # scipy 1.18.0 / numpy 2.5.2: CubicHermiteSpline(X, Y, S).integrate(0, 6)
    assert_close(first.integrate(0.0, 6.0), 14.385416666666666)
    assert_close(first.integrate(6.0, 0.0), -14.385416666666666)
    assert_equal(first.integrate(3.0, 3.0), 0.0)
    assert_true(first == second)
    assert_true(first != clamp)
    assert_true(first != different_slopes)
    assert_true(
        String(first).startswith(
            "CubicHermiteInterpolator(6 knots on [0.0, 6.0], extrapolation=ERROR)"
        )
    )


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
