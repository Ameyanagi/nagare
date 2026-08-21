from nagare import ExtrapolationPolicy, StepInterpolator, StepMode
from std.collections import List
from std.testing import TestSuite, assert_equal, assert_raises, assert_true


def step_fixture(
    mode: StepMode = StepMode.PREVIOUS,
    policy: ExtrapolationPolicy = ExtrapolationPolicy.ERROR,
) raises -> StepInterpolator:
    return StepInterpolator(
        [0.0, 1.0, 3.0, 4.0],
        [10.0, 20.0, 30.0, 40.0],
        mode=mode,
        extrapolation=policy,
    )


def metadata_matches_fixture(interpolator: StepInterpolator) -> Bool:
    """Exercise every metadata accessor from a non-raising function."""
    return (
        interpolator.knot_count() == 4
        and interpolator.domain_start() == 0.0
        and interpolator.domain_end() == 4.0
    )


def test_previous_mode_is_zero_order_hold_and_reproduces_knots() raises:
    var interpolator = step_fixture()
    assert_true(metadata_matches_fixture(interpolator))
    interpolator.validate()

    assert_equal(interpolator.evaluate(0.0), 10.0)
    assert_equal(interpolator.evaluate(1.0), 20.0)
    assert_equal(interpolator.evaluate(3.0), 30.0)
    assert_equal(interpolator.evaluate(4.0), 40.0)
    assert_equal(interpolator.evaluate(0.25), 10.0)
    assert_equal(interpolator.evaluate(2.0), 20.0)
    assert_equal(interpolator.evaluate(3.75), 30.0)


def test_next_mode_selects_following_value_and_reproduces_knots() raises:
    var interpolator = step_fixture(StepMode.NEXT)
    assert_equal(interpolator.evaluate(0.0), 10.0)
    assert_equal(interpolator.evaluate(1.0), 20.0)
    assert_equal(interpolator.evaluate(3.0), 30.0)
    assert_equal(interpolator.evaluate(4.0), 40.0)
    assert_equal(interpolator.evaluate(0.25), 20.0)
    assert_equal(interpolator.evaluate(2.0), 30.0)
    assert_equal(interpolator.evaluate(3.75), 40.0)


def test_nearest_mode_uses_distance_and_selects_left_on_ties() raises:
    var interpolator = step_fixture(StepMode.NEAREST)
    assert_equal(interpolator.evaluate(0.0), 10.0)
    assert_equal(interpolator.evaluate(1.0), 20.0)
    assert_equal(interpolator.evaluate(3.0), 30.0)
    assert_equal(interpolator.evaluate(4.0), 40.0)

    assert_equal(interpolator.evaluate(0.25), 10.0)
    assert_equal(interpolator.evaluate(0.75), 20.0)
    assert_equal(interpolator.evaluate(2.5), 30.0)
    assert_equal(interpolator.evaluate(3.75), 40.0)

    # Exact midpoint ties in intervals of different widths all select left.
    assert_equal(interpolator.evaluate(0.5), 10.0)
    assert_equal(interpolator.evaluate(2.0), 20.0)
    assert_equal(interpolator.evaluate(3.5), 30.0)


def test_error_clamp_and_fill_extrapolation_policies() raises:
    var error = step_fixture()
    with assert_raises(contains="pass extrapolation= to allow this"):
        _ = error.evaluate(-0.25)
    with assert_raises(contains="query 4.25 is outside the knot domain [0.0, 4.0]"):
        _ = error.evaluate(4.25)

    var clamp = step_fixture(policy=ExtrapolationPolicy.CLAMP)
    assert_equal(clamp.evaluate(-10.0), 10.0)
    assert_equal(clamp.evaluate(10.0), 40.0)

    var default_fill = step_fixture(policy=ExtrapolationPolicy.FILL)
    assert_true(default_fill.evaluate(-1.0) != default_fill.evaluate(-1.0))
    assert_true(default_fill.evaluate(5.0) != default_fill.evaluate(5.0))

    var custom_fill = step_fixture(policy=ExtrapolationPolicy.fill(-1.5))
    assert_equal(custom_fill.evaluate(-1.0), -1.5)
    assert_equal(custom_fill.evaluate(5.0), -1.5)
    assert_equal(custom_fill.evaluate(1.0), 20.0)

    with assert_raises(contains="query must be finite: received nan"):
        _ = custom_fill.evaluate(Float64("nan"))


def test_linear_extrapolation_is_rejected_at_construction_and_validation() raises:
    with assert_raises(contains="LINEAR extrapolation is undefined"):
        _ = step_fixture(policy=ExtrapolationPolicy.LINEAR)

    var interpolator = step_fixture()
    interpolator._extrapolation = ExtrapolationPolicy.LINEAR
    with assert_raises(contains="use ERROR, CLAMP, or FILL"):
        interpolator.validate()


def test_constructor_uses_shared_table_validation() raises:
    with assert_raises(contains="len(knots) = 2, len(values) = 1"):
        _ = StepInterpolator([0.0, 1.0], [10.0])
    with assert_raises(contains="knots[1] = 0.0 <= knots[0] = 0.0"):
        _ = StepInterpolator([0.0, 0.0], [10.0, 20.0])
    with assert_raises(contains="values[1] is nan"):
        _ = StepInterpolator([0.0, 1.0], [10.0, Float64("nan")])


def test_mode_constants_are_distinct() raises:
    assert_true(StepMode.PREVIOUS != StepMode.NEXT)
    assert_true(StepMode.PREVIOUS != StepMode.NEAREST)
    assert_true(StepMode.NEXT != StepMode.NEAREST)


def test_batch_and_evaluate_into_agree() raises:
    var queries: List[Float64] = [-1.0, 0.0, 0.75, 2.0, 4.0, 5.0]
    var interpolator = step_fixture(StepMode.NEAREST, ExtrapolationPolicy.fill(-1.5))
    var results = List[Float64](length=len(queries), fill=99.0)
    interpolator.evaluate_into(queries, results)

    var allocated = interpolator.evaluate(queries)
    var called = interpolator(queries)
    var expected: List[Float64] = [-1.5, 10.0, 20.0, 20.0, 40.0, -1.5]
    assert_equal(len(results), len(expected))
    for index in range(len(results)):
        assert_equal(results[index], expected[index])
        assert_equal(allocated[index], results[index])
        assert_equal(called[index], results[index])

    assert_equal(interpolator(0.75), interpolator.evaluate(0.75))

    var short_results = List[Float64](length=2, fill=0.0)
    with assert_raises(contains="len(queries) = 6, len(results) = 2"):
        interpolator.evaluate_into(queries, short_results)


def test_step_equality_and_writable_shapes() raises:
    var first = step_fixture()
    var second = step_fixture()
    var next = step_fixture(StepMode.NEXT)
    var clamp = step_fixture(policy=ExtrapolationPolicy.CLAMP)
    assert_true(first == second)
    assert_true(first != next)
    assert_true(first != clamp)
    assert_equal(String(StepMode.PREVIOUS), "PREVIOUS")
    assert_equal(String(StepMode.NEXT), "NEXT")
    assert_equal(String(StepMode.NEAREST), "NEAREST")
    assert_true(
        String(first).startswith(
            "StepInterpolator(4 knots on [0.0, 4.0], mode=PREVIOUS, "
            "extrapolation=ERROR)"
        )
    )


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
