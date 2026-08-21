from nagare import AkimaInterpolator, ExtrapolationPolicy
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
) raises -> AkimaInterpolator:
    return AkimaInterpolator(
        [0.0, 1.0, 2.5, 3.0, 4.5, 6.0],
        [0.0, 2.0, 1.0, 3.5, 3.0, 4.0],
        extrapolation=policy,
    )


def metadata_matches_reference(interpolator: AkimaInterpolator) -> Bool:
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
    # Akima1DInterpolator(X, Y, method="makima")(Q)
    var expected_values: List[Float64] = [
        0.6583835341365462,
        1.2371709058456046,
        2.0538225088907818,
        1.3215039281705951,
        3.2338945005611675,
        3.268489456185627,
        2.9697179791216204,
        3.313637548891786,
        3.908509343763581,
    ]
    # Reference values from scipy 1.18.0 (numpy 2.5.2):
    # Akima1DInterpolator(X, Y, method="makima").derivative()(Q)
    var expected_derivatives: List[Float64] = [
        2.500334672021419,
        2.103971441320839,
        -0.2477757494624961,
        5.035914702581372,
        4.620426487093159,
        -0.6991281300132525,
        0.16585059635907173,
        0.6069534984789221,
        0.8941329856584095,
    ]

    for index in range(len(queries)):
        assert_close(interpolator.evaluate(queries[index]), expected_values[index])
        assert_close(
            interpolator.derivative(queries[index]), expected_derivatives[index]
        )


def test_metadata_exact_knots_and_affine_reproduction() raises:
    var interpolator = reference_interpolator()
    assert_true(metadata_matches_reference(interpolator))
    interpolator.validate()

    var knots: List[Float64] = [0.0, 1.0, 2.5, 3.0, 4.5, 6.0]
    var values: List[Float64] = [0.0, 2.0, 1.0, 3.5, 3.0, 4.0]
    for index in range(len(knots)):
        assert_equal(interpolator.evaluate(knots[index]), values[index])

    var affine = AkimaInterpolator(
        [-2.0, -0.5, 1.0, 4.0, 8.0],
        [-3.0, 0.0, 3.0, 9.0, 17.0],
    )
    var queries: List[Float64] = [-2.0, -1.25, 0.0, 2.5, 6.0, 8.0]
    for query in queries:
        assert_close(affine.evaluate(query), 2.0 * query + 1.0)
        assert_close(affine.derivative(query), 2.0)


def test_tiny_weight_guard_is_relative_to_the_global_maximum() raises:
    # Unit secants at the left have a nonzero local weighted estimate, but the
    # remote 1e12 secants make their combined weight fall below scipy's global
    # 1e-9 relative threshold. The specified guarded slope is exactly zero.
    var interpolator = AkimaInterpolator(
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 1000000000005.0, 2000000000005.0],
    )
    assert_equal(interpolator.derivative(0.0), 0.0)


def test_constructor_requires_four_valid_knots() raises:
    with assert_raises(contains="at least four knots: received 0"):
        _ = AkimaInterpolator(List[Float64](), List[Float64]())
    with assert_raises(contains="at least four knots: received 3"):
        _ = AkimaInterpolator([0.0, 1.0, 2.0], [1.0, 2.0, 3.0])
    with assert_raises(contains="len(knots) = 4, len(values) = 3"):
        _ = AkimaInterpolator([0.0, 1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
    with assert_raises(contains="knots[2] = 1.0 <= knots[1] = 1.0"):
        _ = AkimaInterpolator(
            [0.0, 1.0, 1.0, 3.0],
            [1.0, 2.0, 3.0, 4.0],
        )
    with assert_raises(contains="values[2] is nan"):
        _ = AkimaInterpolator(
            [0.0, 1.0, 2.0, 3.0],
            [1.0, 2.0, Float64("nan"), 4.0],
        )


def test_explicit_validate_rechecks_derived_slopes() raises:
    var valid = reference_interpolator()
    valid.validate()

    var non_finite = reference_interpolator()
    non_finite._engine._slopes[3] = Float64("inf")
    with assert_raises(contains="slopes[3] is inf"):
        non_finite.validate()

    var short = reference_interpolator()
    short._engine._knots = [0.0, 1.0, 2.0]
    short._engine._values = [0.0, 1.0, 2.0]
    short._engine._slopes = [1.0, 1.0, 1.0]
    with assert_raises(contains="at least four knots: received 3"):
        short.validate()


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

    var linear = reference_interpolator(ExtrapolationPolicy.LINEAR)
    var left_slope = linear.derivative(0.0)
    var right_slope = linear.derivative(6.0)
    assert_close(linear.evaluate(-1.0), -left_slope)
    assert_close(linear.evaluate(7.0), 4.0 + right_slope)
    assert_equal(linear.derivative(-1.0), left_slope)
    assert_equal(linear.derivative(7.0), right_slope)

    var default_fill = reference_interpolator(ExtrapolationPolicy.FILL)
    assert_true(default_fill.evaluate(-1.0) != default_fill.evaluate(-1.0))
    assert_true(default_fill.derivative(7.0) != default_fill.derivative(7.0))

    var custom_fill = reference_interpolator(ExtrapolationPolicy.fill(-4.5))
    assert_equal(custom_fill.evaluate(-1.0), -4.5)
    assert_equal(custom_fill.evaluate(7.0), -4.5)
    assert_equal(custom_fill.derivative(-1.0), -4.5)
    assert_equal(custom_fill.derivative(7.0), -4.5)


def test_non_finite_queries_always_raise() raises:
    var interpolator = reference_interpolator(ExtrapolationPolicy.FILL)
    with assert_raises(contains="query must be finite"):
        _ = interpolator.evaluate(Float64("nan"))
    with assert_raises(contains="query must be finite"):
        _ = interpolator.derivative(Float64("inf"))


def test_span_wrapper_and_evaluate_into_agree() raises:
    var queries: List[Float64] = [-1.0, 0.0, 0.5, 2.6, 4.4, 6.0, 7.0]
    var interpolator = reference_interpolator(ExtrapolationPolicy.fill(-4.5))
    var results = List[Float64](length=len(queries), fill=99.0)
    interpolator.evaluate_into(queries, results)
    var allocated = interpolator.evaluate(queries)

    for index in range(len(queries)):
        assert_equal(results[index], allocated[index])
        assert_equal(results[index], interpolator.evaluate(queries[index]))

    var short_results = List[Float64](length=3, fill=0.0)
    with assert_raises(contains="len(queries) = 7, len(results) = 3"):
        interpolator.evaluate_into(queries, short_results)


def test_closed_form_integral_and_equality_writable_surface() raises:
    var first = reference_interpolator()
    var second = reference_interpolator()
    var clamp = reference_interpolator(ExtrapolationPolicy.CLAMP)
    # scipy 1.18.0 / numpy 2.5.2:
    # Akima1DInterpolator(X, Y, method="makima").integrate(0, 6)
    assert_close(first.integrate(0.0, 6.0), 14.525352304928946)
    assert_true(first == second)
    assert_true(first != clamp)
    assert_true(
        String(first).startswith(
            "AkimaInterpolator(6 knots on [0.0, 6.0], extrapolation=ERROR)"
        )
    )


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
