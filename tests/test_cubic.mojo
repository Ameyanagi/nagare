from nagare import (
    BoundaryCondition,
    CubicSplineInterpolator,
    ExtrapolationPolicy,
    LinearInterpolator,
)
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
        boundary=BoundaryCondition.NATURAL,
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
    var expected_knots: List[Float64] = [0.0, 1.0, 2.0, 3.0]
    var expected_values: List[Float64] = [0.0, 1.0, 0.0, 1.0]
    var knots = interpolator.knots()
    var values = interpolator.values()
    for index in range(len(expected_knots)):
        assert_equal(
            interpolator.evaluate(expected_knots[index]), expected_values[index]
        )
        assert_equal(knots[index], expected_knots[index])
        assert_equal(values[index], expected_values[index])


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

    with assert_raises(contains="query must be finite: received nan"):
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
    with assert_raises(
        contains="expected 3, len(a) = 1, len(b) = 3, len(c) = 3, len(d) = 3"
    ):
        short_coefficients.validate()

    var non_finite_coefficients = hand_fixture()
    non_finite_coefficients._d[1] = Float64("nan")
    with assert_raises(contains="coefficients must be finite: d[1] is nan"):
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


def test_batch_derivative_matches_scalar_and_checks_buffer_lengths() raises:
    var queries: List[Float64] = [0.0, 0.5, 2.0, 3.0]
    var interpolator = hand_fixture()
    var allocated = interpolator.derivative(queries)
    var results = List[Float64](length=len(queries), fill=99.0)
    interpolator.derivative_into(queries, results)
    for index in range(len(queries)):
        assert_equal(allocated[index], interpolator.derivative(queries[index]))
        assert_equal(results[index], allocated[index])

    var short_results = List[Float64](length=1, fill=0.0)
    with assert_raises(
        contains=(
            "query and result buffers must have equal length: "
            "len(queries) = 4, len(results) = 1"
        )
    ):
        interpolator.derivative_into(queries, short_results)


def test_two_knot_spline_degenerates_to_a_straight_segment() raises:
    var spline = CubicSplineInterpolator([-2.0, 3.0], [4.0, -6.0])
    var line = LinearInterpolator([-2.0, 3.0], [4.0, -6.0])
    var queries: List[Float64] = [-2.0, -1.5, 0.0, 1.25, 3.0]
    for query in queries:
        assert_close(spline.evaluate(query), line.evaluate(query))
        assert_close(spline.second_derivative(query), 0.0)


def shared_table(
    boundary: BoundaryCondition = BoundaryCondition.NOT_A_KNOT,
) raises -> CubicSplineInterpolator:
    return CubicSplineInterpolator(
        [0.0, 1.0, 2.5, 3.0, 4.5, 6.0],
        [0.0, 2.0, 1.0, 3.5, 3.0, 4.0],
        boundary=boundary,
    )


def test_default_not_a_knot_matches_scipy_values_derivatives_and_integral() raises:
    var interpolator = shared_table()
    assert_true(interpolator.boundary() == BoundaryCondition.NOT_A_KNOT)
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
    # scipy 1.18.0 / numpy 2.5.2: CubicSpline(X, Y)(Q)
    var expected_values: List[Float64] = [
        1.3992578125,
        2.0954166666666665,
        1.6268600000000002,
        1.4436066666666674,
        3.026676666666666,
        4.555117777777777,
        3.2432088888888884,
        1.9212549999999997,
        3.4216838888888903,
    ]
    # scipy 1.18.0 / numpy 2.5.2: CubicSpline(X, Y).derivative()(Q)
    var expected_derivatives: List[Float64] = [
        4.087239583333333,
        1.585625,
        -2.1121833333333333,
        4.853950000000001,
        5.076825000000001,
        -0.7446083333333333,
        -2.4425166666666662,
        -0.7733083333333355,
        5.25735833333334,
    ]
    for index in range(len(queries)):
        assert_close(
            interpolator.evaluate(queries[index]),
            expected_values[index],
            atol=1e-9,
            rtol=1e-12,
        )
        assert_close(
            interpolator.derivative(queries[index]),
            expected_derivatives[index],
            atol=1e-9,
            rtol=1e-12,
        )
    assert_close(interpolator.integrate(0.0, 6.0), 13.815937499999999, atol=1e-9)


def test_natural_and_clamped_scipy_fixtures() raises:
    var natural = shared_table(BoundaryCondition.NATURAL)
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
    # scipy 1.18.0 / numpy 2.5.2: CubicSpline(X, Y, bc_type="natural")(Q)
    var natural_values: List[Float64] = [
        0.7711738484398216,
        1.433878157503715,
        1.8156138902647077,
        1.4380128776622096,
        3.0317543338286277,
        4.32340278465687,
        3.1374253480821066,
        2.861372956909361,
        3.840237301194212,
    ]
    for index in range(len(queries)):
        assert_close(
            natural.evaluate(queries[index]),
            natural_values[index],
            atol=1e-9,
            rtol=1e-12,
        )
    assert_close(natural.integrate(0.0, 6.0), 14.493024599636785, atol=1e-9)

    var clamped_boundary = BoundaryCondition.clamped(0.5, -1.0)
    var clamped = shared_table(clamped_boundary)
    # scipy 1.18.0 / numpy 2.5.2:
    # CubicSpline(X, Y, bc_type=((1, 0.5), (1, -1.0)))(Q)
    var clamped_values: List[Float64] = [
        0.35939777696793007,
        1.0000607385811469,
        1.9397227081308714,
        1.4335898931000977,
        3.0337764820213793,
        4.1974300831443685,
        3.079304610733182,
        3.3793411078717197,
        4.071009826152684,
    ]
    # scipy 1.18.0 / numpy 2.5.2: the same spline's derivative()(Q)
    var clamped_derivatives: List[Float64] = [
        2.187651846452867,
        2.7501214771622937,
        -1.0029931972789115,
        4.835879494655005,
        5.080631681243927,
        -1.2340945902170402,
        -1.0269517330741804,
        1.3448104956268212,
        -0.43296728215095737,
    ]
    for index in range(len(queries)):
        assert_close(
            clamped.evaluate(queries[index]),
            clamped_values[index],
            atol=1e-9,
            rtol=1e-12,
        )
        assert_close(
            clamped.derivative(queries[index]),
            clamped_derivatives[index],
            atol=1e-9,
            rtol=1e-12,
        )
    assert_close(clamped.derivative(0.0), 0.5, atol=1e-12)
    assert_close(clamped.derivative(6.0), -1.0, atol=1e-12)
    assert_close(clamped.integrate(0.0, 6.0), 14.849186103012634, atol=1e-9)


def test_periodic_scipy_fixture_and_endpoint_continuity() raises:
    var interpolator = CubicSplineInterpolator(
        [0.0, 1.0, 2.5, 3.0, 4.5, 6.0],
        [0.0, 2.0, 1.0, 3.5, 3.0, 0.0],
        boundary=BoundaryCondition.PERIODIC,
    )
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
    # scipy 1.18.0 / numpy 2.5.2: CubicSpline(X, YP, bc_type="periodic")(Q)
    var expected_values: List[Float64] = [
        0.45625000000000004,
        1.0999999999999999,
        1.917511111111111,
        1.4200000000000004,
        3.0079999999999996,
        4.729303703703703,
        3.3104296296296285,
        1.096800000000001,
        -0.0802074074074075,
    ]
    # scipy 1.18.0 / numpy 2.5.2: the same periodic spline's derivative()(Q)
    var expected_derivatives: List[Float64] = [
        2.358333333333333,
        2.6333333333333333,
        -1.0639999999999998,
        4.693333333333334,
        5.233333333333334,
        -0.5782222222222231,
        -3.0328888888888894,
        -2.862666666666667,
        0.4795555555555586,
    ]
    for index in range(len(queries)):
        assert_close(
            interpolator.evaluate(queries[index]),
            expected_values[index],
            atol=1e-9,
            rtol=1e-12,
        )
        assert_close(
            interpolator.derivative(queries[index]),
            expected_derivatives[index],
            atol=1e-9,
            rtol=1e-12,
        )
    assert_close(interpolator.integrate(0.0, 6.0), 11.566666666666666, atol=1e-9)
    assert_close(interpolator.derivative(0.0), interpolator.derivative(6.0), atol=1e-9)
    assert_close(
        interpolator.second_derivative(0.0),
        interpolator.second_derivative(6.0),
        atol=1e-9,
    )

    with assert_raises(contains="values[0] = 0.0, values[2] = 2.0"):
        _ = CubicSplineInterpolator(
            [0.0, 1.0, 2.0],
            [0.0, 1.0, 2.0],
            boundary=BoundaryCondition.PERIODIC,
        )
    with assert_raises(contains="at least three knots"):
        _ = CubicSplineInterpolator(
            [0.0, 1.0],
            [2.0, 2.0],
            boundary=BoundaryCondition.PERIODIC,
        )


def test_not_a_knot_three_points_form_one_parabola() raises:
    var interpolator = CubicSplineInterpolator(
        [0.0, 1.0, 3.0],
        [1.0, 2.0, 10.0],
    )
    # The unique parabola is x^2 + 1.
    assert_close(interpolator.evaluate(0.5), 1.25)
    assert_close(interpolator.evaluate(2.0), 5.0)
    assert_close(interpolator.second_derivative(0.0), 2.0)
    assert_close(interpolator.second_derivative(3.0), 2.0)


def test_clamped_two_knots_and_periodic_three_knots() raises:
    var clamped = CubicSplineInterpolator(
        [-2.0, 3.0],
        [4.0, -6.0],
        boundary=BoundaryCondition.clamped(0.5, -1.0),
    )
    assert_close(clamped.derivative(-2.0), 0.5, atol=1e-12)
    assert_close(clamped.derivative(3.0), -1.0, atol=1e-12)

    var periodic = CubicSplineInterpolator(
        [0.0, 1.0, 2.0],
        [0.0, 1.0, 0.0],
        boundary=BoundaryCondition.PERIODIC,
    )
    assert_close(periodic.derivative(0.0), periodic.derivative(2.0), atol=1e-9)
    assert_close(
        periodic.second_derivative(0.0),
        periodic.second_derivative(2.0),
        atol=1e-9,
    )


def test_default_not_a_knot_reproduces_cubic_polynomial_integrals() raises:
    # f(x) = 2x^3 - 3x^2 + x - 5. The exact antiderivative is
    # F(x) = x^4/2 - x^3 + x^2/2 - 5x; values were also cross-checked against
    # scipy 1.18.0 / numpy 2.5.2 CubicSpline.integrate.
    var interpolator = CubicSplineInterpolator(
        [0.0, 1.0, 2.0, 3.0, 4.0],
        [-5.0, -5.0, 1.0, 25.0, 79.0],
    )
    assert_close(interpolator.integrate(0.0, 4.0), 52.0, atol=1e-9)
    assert_close(interpolator.integrate(0.3, 3.7), 32.87800000000001, atol=1e-9)


def test_boundary_validation_equality_and_writable_shape() raises:
    assert_true(BoundaryCondition.NOT_A_KNOT != BoundaryCondition.NATURAL)
    assert_true(BoundaryCondition.PERIODIC != BoundaryCondition.NATURAL)
    assert_true(
        BoundaryCondition.clamped(0.5, -1.0) == BoundaryCondition.clamped(0.5, -1.0)
    )
    assert_true(
        BoundaryCondition.clamped(0.5, -1.0) != BoundaryCondition.clamped(0.5, 1.0)
    )
    with assert_raises(contains="clamped boundary slopes must be finite"):
        _ = BoundaryCondition.clamped(Float64("nan"), 0.0)

    var first = shared_table()
    var second = shared_table()
    var natural = shared_table(BoundaryCondition.NATURAL)
    assert_true(first == second)
    assert_true(first != natural)
    assert_equal(String(BoundaryCondition.NOT_A_KNOT), "NOT_A_KNOT")
    assert_equal(String(BoundaryCondition.NATURAL), "NATURAL")
    assert_equal(String(BoundaryCondition.PERIODIC), "PERIODIC")
    assert_equal(
        String(BoundaryCondition.clamped(0.5, -1.0)),
        "CLAMPED(start=0.5, end=-1.0)",
    )
    assert_true(
        String(first).startswith(
            "CubicSplineInterpolator(6 knots on [0.0, 6.0], "
            "boundary=NOT_A_KNOT, extrapolation=ERROR)"
        )
    )


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
