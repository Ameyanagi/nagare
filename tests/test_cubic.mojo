from nagare import CubicSplineInterpolator, ExtrapolationPolicy, LinearInterpolator
from std.collections import List
from std.testing import TestSuite, assert_equal, assert_raises, assert_true


# Tolerances used throughout this file:
# - exact equality for knot reproduction and clamped endpoint values;
# - atol=1e-12, rtol=0 for hand-derived and affine fixtures below magnitude 32;
# - atol=1e-6, rtol=0 for two-sided limit checks at distance h=1e-8.
#   On this fixture the largest relevant derivative scale is below 16, so the
#   expected O(h) separation is below 3.2e-7; 1e-6 allows ordinary rounding
#   without concealing a discontinuity.


def assert_close(
    actual: Float64,
    expected: Float64,
    atol: Float64 = 1e-12,
    rtol: Float64 = 0.0,
) raises:
    var tolerance = max(atol, rtol * abs(expected))
    assert_true(
        abs(actual - expected) <= tolerance,
        msg=String("expected ", expected, ", got ", actual),
    )


def hand_fixture(
    policy: ExtrapolationPolicy = ExtrapolationPolicy.ERROR,
) raises -> CubicSplineInterpolator:
    return CubicSplineInterpolator(
        [0.0, 1.0, 2.0, 3.0],
        [0.0, 1.0, 0.0, 1.0],
        extrapolation=policy,
    )


def metadata_matches_fixture(interpolator: CubicSplineInterpolator) -> Bool:
    """Exercise every metadata accessor from a non-raising function."""
    return (
        interpolator.knot_count() == 4
        and interpolator.domain_start() == 0.0
        and interpolator.domain_end() == 3.0
    )


def test_metadata_and_exact_knot_reproduction() raises:
    var interpolator = hand_fixture()
    assert_true(metadata_matches_fixture(interpolator))
    assert_equal(interpolator.evaluate(0.0), 0.0)
    assert_equal(interpolator.evaluate(1.0), 1.0)
    assert_equal(interpolator.evaluate(2.0), 0.0)
    assert_equal(interpolator.evaluate(3.0), 1.0)


def test_hand_derived_natural_spline_values() raises:
    # With unit knot spacing and natural moments M0=M3=0, the interior system
    # derived from continuity of S' is:
    #
    #     4 M1 + M2 = 6((-1) - 1) = -12
    #       M1 + 4 M2 = 6(1 - (-1)) = 12
    #
    # Hence M1=-4 and M2=4. Substitution into
    #   b_i = slope_i - h_i(2 M_i + M_{i+1})/6,
    #   c_i = M_i/2, d_i = (M_{i+1}-M_i)/(6 h_i)
    # gives the independently evaluated shifted segments:
    #   S0(t) = (5/3)t - (2/3)t^3
    #   S1(t) = 1 - (1/3)t - 2t^2 + (4/3)t^3
    #   S2(t) = -(1/3)t + 2t^2 - (2/3)t^3.
    var interpolator = hand_fixture()
    assert_close(interpolator.evaluate(0.5), 3.0 / 4.0)
    assert_close(interpolator.evaluate(1.25), 13.0 / 16.0)
    assert_close(interpolator.evaluate(1.5), 1.0 / 2.0)
    assert_close(interpolator.evaluate(2.5), 1.0 / 4.0)
    assert_close(interpolator.second_derivative(1.0), -4.0)
    assert_close(interpolator.second_derivative(2.0), 4.0)


def test_affine_data_is_reproduced_over_irregular_knots() raises:
    var interpolator = CubicSplineInterpolator(
        [-4.0, -0.5, 0.0, 1.25, 8.0],
        [-14.0, -3.5, -2.0, 1.75, 22.0],
    )
    for step in range(97):
        var x = -4.0 + Float64(step) * 0.125
        assert_close(interpolator.evaluate(x), 3.0 * x - 2.0)

    var probes: List[Float64] = [-4.0, -2.25, -0.5, 0.625, 4.0, 8.0]
    for probe in probes:
        assert_close(interpolator.second_derivative(probe), 0.0)


def test_c0_c1_c2_continuity_at_every_interior_knot() raises:
    var interpolator = hand_fixture()
    var knots: List[Float64] = [1.0, 2.0]
    var h = 1e-8
    var continuity_atol = 1e-6
    for knot in knots:
        assert_close(
            interpolator.evaluate(knot - h),
            interpolator.evaluate(knot + h),
            atol=continuity_atol,
        )
        assert_close(
            interpolator.derivative(knot - h),
            interpolator.derivative(knot + h),
            atol=continuity_atol,
        )
        assert_close(
            interpolator.second_derivative(knot - h),
            interpolator.second_derivative(knot + h),
            atol=continuity_atol,
        )


def test_natural_boundary_second_derivatives_are_zero() raises:
    var interpolator = hand_fixture()
    assert_close(interpolator.second_derivative(interpolator.domain_start()), 0.0)
    assert_close(interpolator.second_derivative(interpolator.domain_end()), 0.0)


def test_error_policy_rejects_both_exterior_regions() raises:
    var interpolator = hand_fixture()
    with assert_raises(contains="outside the knot domain"):
        _ = interpolator.evaluate(-0.25)
    with assert_raises(contains="outside the knot domain"):
        _ = interpolator.evaluate(3.25)
    with assert_raises(contains="outside the knot domain"):
        _ = interpolator.derivative(-0.25)
    with assert_raises(contains="outside the knot domain"):
        _ = interpolator.derivative(3.25)
    with assert_raises(contains="outside the knot domain"):
        _ = interpolator.second_derivative(-0.25)
    with assert_raises(contains="outside the knot domain"):
        _ = interpolator.second_derivative(3.25)


def test_clamp_policy_is_constant_outside_both_ends() raises:
    var interpolator = hand_fixture(ExtrapolationPolicy.CLAMP)
    assert_equal(interpolator.evaluate(-2.0), 0.0)
    assert_equal(interpolator.evaluate(5.0), 1.0)
    assert_equal(interpolator.derivative(-2.0), 0.0)
    assert_equal(interpolator.derivative(5.0), 0.0)
    assert_equal(interpolator.second_derivative(-2.0), 0.0)
    assert_equal(interpolator.second_derivative(5.0), 0.0)


def test_linear_policy_uses_endpoint_tangent_rays() raises:
    # The hand-derived first and final segments both have endpoint slope 5/3.
    var interpolator = hand_fixture(ExtrapolationPolicy.LINEAR)
    var endpoint_slope = 5.0 / 3.0
    assert_close(interpolator.evaluate(-1.0), 0.0 + endpoint_slope * -1.0)
    assert_close(interpolator.evaluate(4.0), 1.0 + endpoint_slope * 1.0)
    assert_close(interpolator.derivative(-1.0), endpoint_slope)
    assert_close(interpolator.derivative(4.0), endpoint_slope)
    assert_equal(interpolator.second_derivative(-1.0), 0.0)
    assert_equal(interpolator.second_derivative(4.0), 0.0)


def test_fill_policy_applies_to_values_and_derivatives_outside_domain() raises:
    var default_fill = hand_fixture(ExtrapolationPolicy.FILL)
    assert_true(default_fill.evaluate(-1.0) != default_fill.evaluate(-1.0))
    assert_true(default_fill.evaluate(4.0) != default_fill.evaluate(4.0))
    assert_true(default_fill.derivative(-1.0) != default_fill.derivative(-1.0))
    assert_true(default_fill.derivative(4.0) != default_fill.derivative(4.0))
    assert_true(
        default_fill.second_derivative(-1.0) != default_fill.second_derivative(-1.0)
    )
    assert_true(
        default_fill.second_derivative(4.0) != default_fill.second_derivative(4.0)
    )

    var custom_fill = hand_fixture(ExtrapolationPolicy.fill(-1.5))
    assert_equal(custom_fill.evaluate(-1.0), -1.5)
    assert_equal(custom_fill.evaluate(4.0), -1.5)
    assert_equal(custom_fill.derivative(-1.0), -1.5)
    assert_equal(custom_fill.derivative(4.0), -1.5)
    assert_equal(custom_fill.second_derivative(-1.0), -1.5)
    assert_equal(custom_fill.second_derivative(4.0), -1.5)

    var reference = hand_fixture()
    assert_equal(custom_fill.evaluate(1.5), reference.evaluate(1.5))
    assert_equal(custom_fill.derivative(1.5), reference.derivative(1.5))
    assert_equal(
        custom_fill.second_derivative(1.5),
        reference.second_derivative(1.5),
    )

    with assert_raises(contains="query must be finite"):
        _ = custom_fill.evaluate(Float64("nan"))
    with assert_raises(contains="query must be finite"):
        _ = custom_fill.derivative(Float64("inf"))
    with assert_raises(contains="query must be finite"):
        _ = custom_fill.second_derivative(Float64("-inf"))


def test_non_finite_query_is_rejected_under_every_policy() raises:
    var error = hand_fixture(ExtrapolationPolicy.ERROR)
    var clamp = hand_fixture(ExtrapolationPolicy.CLAMP)
    var linear = hand_fixture(ExtrapolationPolicy.LINEAR)
    with assert_raises(contains="query must be finite"):
        _ = error.evaluate(Float64("nan"))
    with assert_raises(contains="query must be finite"):
        _ = clamp.evaluate(Float64("inf"))
    with assert_raises(contains="query must be finite"):
        _ = linear.evaluate(Float64("-inf"))
    with assert_raises(contains="query must be finite"):
        _ = error.derivative(Float64("nan"))
    with assert_raises(contains="query must be finite"):
        _ = linear.second_derivative(Float64("inf"))


def test_constructor_rejects_invalid_tables() raises:
    with assert_raises(contains="equal length"):
        _ = CubicSplineInterpolator([0.0, 1.0], [2.0])
    with assert_raises(contains="strictly increasing"):
        _ = CubicSplineInterpolator([0.0, 0.0], [1.0, 2.0])
    with assert_raises(contains="values must be finite"):
        _ = CubicSplineInterpolator([0.0, 1.0], [1.0, Float64("nan")])


def test_explicit_validate_rechecks_table_and_coefficients() raises:
    var valid = hand_fixture()
    valid.validate()

    var duplicate_knots = hand_fixture()
    duplicate_knots._knots[1] = duplicate_knots._knots[0]
    with assert_raises(contains="strictly increasing"):
        duplicate_knots.validate()

    var non_finite_values = hand_fixture()
    non_finite_values._values[1] = Float64("inf")
    with assert_raises(contains="values must be finite"):
        non_finite_values.validate()

    var short_coefficients = hand_fixture()
    short_coefficients._a = [0.0]
    with assert_raises(contains="buffers must match interval count"):
        short_coefficients.validate()

    var non_finite_coefficients = hand_fixture()
    non_finite_coefficients._d[1] = Float64("nan")
    with assert_raises(contains="coefficients must be finite"):
        non_finite_coefficients.validate()


def test_batch_matches_scalar_and_handles_errors_and_empty_input() raises:
    var queries: List[Float64] = [
        -1.0,
        -0.01,
        0.0,
        0.25,
        0.999,
        1.0,
        1.5,
        2.0,
        2.75,
        3.0,
        3.01,
        4.0,
    ]
    var clamp = hand_fixture(ExtrapolationPolicy.CLAMP)
    var clamped = clamp.evaluate(queries)
    for index in range(len(queries)):
        assert_equal(clamped[index], clamp.evaluate(queries[index]))

    var linear = hand_fixture(ExtrapolationPolicy.LINEAR)
    var extended = linear.evaluate(queries)
    for index in range(len(queries)):
        assert_equal(extended[index], linear.evaluate(queries[index]))

    var error = hand_fixture(ExtrapolationPolicy.ERROR)
    with assert_raises(contains="outside the knot domain"):
        _ = error.evaluate(queries)

    var empty = List[Float64]()
    var empty_results = error.evaluate(empty)
    assert_equal(len(empty_results), 0)


def test_evaluate_into_matches_allocating_wrapper() raises:
    var queries: List[Float64] = [-1.0, 0.0, 0.5, 2.0, 3.0, 4.0]
    var interpolator = hand_fixture(ExtrapolationPolicy.fill(-1.5))
    var results = List[Float64](length=len(queries), fill=99.0)
    interpolator.evaluate_into(queries, results)

    var allocated = interpolator.evaluate(queries)
    assert_equal(len(results), len(allocated))
    for index in range(len(results)):
        assert_equal(results[index], allocated[index])

    var short_results = List[Float64](length=1, fill=0.0)
    with assert_raises(contains="len(queries) = 6, len(results) = 1"):
        interpolator.evaluate_into(queries, short_results)


def test_two_knot_spline_degenerates_to_a_straight_segment() raises:
    var spline = CubicSplineInterpolator([-2.0, 3.0], [4.0, -6.0])
    var line = LinearInterpolator([-2.0, 3.0], [4.0, -6.0])
    var queries: List[Float64] = [-2.0, -1.5, 0.0, 1.25, 3.0]
    for query in queries:
        assert_close(spline.evaluate(query), line.evaluate(query))
        assert_close(spline.second_derivative(query), 0.0)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
