from nagare import locate_interval
from std.collections import List
from std.testing import TestSuite, assert_equal, assert_raises, assert_true


def test_locate_interval_boundaries_and_right_bias() raises:
    var knots: List[Float64] = [0.0, 1.0, 3.0, 4.0]

    assert_equal(locate_interval(knots, 0.0), 0)
    assert_equal(locate_interval(knots, 0.5), 0)
    assert_equal(locate_interval(knots, 1.0), 1)
    assert_equal(locate_interval(knots, 2.0), 1)
    assert_equal(locate_interval(knots, 3.0), 2)
    assert_equal(locate_interval(knots, 4.0), 2)


def test_locate_interval_covers_every_interval() raises:
    var knots: List[Float64] = [-3.0, -1.5, -0.25, 2.0, 9.0, 20.0]
    for index in range(len(knots) - 1):
        var midpoint = (knots[index] + knots[index + 1]) / 2.0
        assert_equal(locate_interval(knots, midpoint), index)
        assert_true(knots[index] <= midpoint)
        assert_true(midpoint <= knots[index + 1])


def test_locate_interval_rejects_invalid_knots() raises:
    with assert_raises(contains="at least two values: received 0"):
        _ = locate_interval(List[Float64](), 0.0)
    with assert_raises(contains="at least two values: received 1"):
        _ = locate_interval([0.0], 0.0)
    with assert_raises(contains="knots[2] = 1.0 <= knots[1] = 1.0"):
        _ = locate_interval([0.0, 1.0, 1.0], 0.5)
    with assert_raises(contains="strictly increasing"):
        _ = locate_interval([0.0, 2.0, 1.0], 0.5)
    with assert_raises(contains="knots[1] is nan"):
        _ = locate_interval([0.0, Float64("nan")], 0.0)
    with assert_raises(contains="knots must be finite"):
        _ = locate_interval([Float64("-inf"), 0.0], 0.0)


def test_locate_interval_rejects_invalid_queries() raises:
    var knots: List[Float64] = [0.0, 1.0, 2.0]
    with assert_raises(contains="query -0.01 is outside the knot domain [0.0, 2.0]"):
        _ = locate_interval(knots, -0.01)
    with assert_raises(contains="outside the knot domain"):
        _ = locate_interval(knots, 2.01)
    with assert_raises(contains="query must be finite: received nan"):
        _ = locate_interval(knots, Float64("nan"))
    with assert_raises(contains="query must be finite: received inf"):
        _ = locate_interval(knots, Float64("inf"))


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
