from linear_fixtures import (
    make_full_range_fixture,
    make_linear_fixture,
    make_tiny_ordinate_fixture,
)
from nagare import ExtrapolationPolicy, LinearInterpolator
from std.testing import TestSuite, assert_equal, assert_raises, assert_true


def _assert_finite(value: Float64) raises:
    assert_true(value == value and value - value == 0.0)


def test_small_uniform_fixture_has_independent_reference_values() raises:
    var fixture = make_linear_fixture(4, False)
    assert_equal(fixture.knots, [-2.0, -1.25, -0.5, 0.25])
    assert_equal(fixture.values, [-5.25, -3.875, -2.5, -1.875])


def test_fixture_sizes_and_order_are_deterministic_properties() raises:
    for count in [8, 1_024, 65_536]:
        for irregular in [False, True]:
            var fixture = make_linear_fixture(count, irregular)
            assert_equal(len(fixture.knots), count)
            assert_equal(len(fixture.values), count)
            for index in range(count):
                _assert_finite(fixture.knots[index])
                _assert_finite(fixture.values[index])
                if index > 0:
                    assert_true(fixture.knots[index] > fixture.knots[index - 1])


def test_mutated_fixture_is_rejected_by_benchmarked_constructor() raises:
    var fixture = make_linear_fixture(8, True)
    fixture.knots[2] = fixture.knots[1]
    with assert_raises(contains="strictly increasing"):
        _ = LinearInterpolator(fixture.knots.copy(), fixture.values.copy())

    fixture = make_linear_fixture(8, True)
    fixture.values[4] = Float64("inf")
    with assert_raises(contains="values must be finite"):
        _ = LinearInterpolator(fixture.knots.copy(), fixture.values.copy())


def test_extreme_fixtures_preserve_documented_numeric_contracts() raises:
    var tiny = make_tiny_ordinate_fixture()
    var tiny_interpolator = LinearInterpolator(
        tiny.knots.copy(),
        tiny.values.copy(),
        extrapolation=ExtrapolationPolicy.LINEAR,
    )
    var extrapolated = tiny_interpolator.evaluate(1.0)
    _assert_finite(extrapolated)
    assert_true(abs(extrapolated - 2.0) <= 2e-15)

    var full_range = make_full_range_fixture()
    var full_range_interpolator = LinearInterpolator(
        full_range.knots.copy(),
        full_range.values.copy(),
    )
    var central = full_range_interpolator.evaluate(1.0)
    _assert_finite(central)
    assert_true(central >= full_range.values[0])
    assert_true(central <= full_range.values[1])
    assert_equal(central, 0.0)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
