from nagare import KnotIndex, locate_interval
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


def test_knot_index_matches_one_shot_boundaries_and_unsorted_queries() raises:
    var knots: List[Float64] = [-3.0, -1.5, -0.25, 2.0, 9.0, 20.0]
    var index = KnotIndex(knots.copy())
    var queries: List[Float64] = [20.0, -3.0, 9.0, -1.5, 2.0, -0.25, 0.0, 9.0]
    var results = List[Int](length=len(queries), fill=-1)
    index.locate_into(queries, results)
    for position in range(len(queries)):
        assert_equal(
            index.locate(queries[position]), locate_interval(knots, queries[position])
        )
        assert_equal(results[position], locate_interval(knots, queries[position]))
    for position in range(len(knots) - 1):
        var midpoint = (knots[position] + knots[position + 1]) / 2.0
        assert_equal(index.locate(midpoint), position)
    var empty_queries = List[Float64]()
    var empty_results = List[Int]()
    index.locate_into(empty_queries, empty_results)
    assert_equal(index.knot_count(), len(knots))
    assert_equal(index.knots()[0], -3.0)
    assert_equal(index.knots()[5], 20.0)
    assert_true(index == index.copy())
    assert_true(index != KnotIndex([0.0, 1.0]))
    assert_equal(String(index), "KnotIndex(knots=6, domain=[-3.0, 20.0])")
    var pair = KnotIndex([1.0, 2.0])
    assert_equal(pair.locate(1.0), 0)
    assert_equal(pair.locate(2.0), 0)


def test_knot_index_rejects_invalid_construction_and_explicit_validation() raises:
    with assert_raises(contains="at least two values: received 0"):
        _ = KnotIndex(List[Float64]())
    with assert_raises(contains="at least two values: received 1"):
        _ = KnotIndex([0.0])
    with assert_raises(contains="strictly increasing"):
        _ = KnotIndex([0.0, 1.0, 1.0])
    with assert_raises(contains="strictly increasing"):
        _ = KnotIndex([0.0, 2.0, 1.0])
    with assert_raises(contains="knots must be finite"):
        _ = KnotIndex([0.0, Float64("nan")])
    with assert_raises(contains="knots must be finite"):
        _ = KnotIndex([Float64("-inf"), 0.0])
    with assert_raises(contains="knots must be finite"):
        _ = KnotIndex([0.0, Float64("inf")])
    var index = KnotIndex([0.0, 1.0, 2.0])
    index.validate()
    index._knots[1] = 2.0
    with assert_raises(contains="strictly increasing"):
        index.validate()


def test_knot_index_rejects_invalid_queries_and_buffer_lengths() raises:
    var index = KnotIndex([0.0, 1.0, 2.0])
    var bad_queries: List[Float64] = [Float64("nan"), Float64("inf"), Float64("-inf")]
    for query in bad_queries:
        with assert_raises(contains="query must be finite"):
            _ = index.locate(query)
    with assert_raises(contains="query -0.01 is outside the knot domain [0.0, 2.0]"):
        _ = index.locate(-0.01)
    with assert_raises(contains="query 2.01 is outside the knot domain [0.0, 2.0]"):
        _ = index.locate(2.01)
    var queries: List[Float64] = [0.0, 1.5, 2.01]
    var results = List[Int](length=3, fill=-1)
    with assert_raises(contains="outside the knot domain"):
        index.locate_into(queries, results)
    queries[2] = Float64("nan")
    with assert_raises(contains="query must be finite"):
        index.locate_into(queries, results)
    var too_short = List[Int](length=2, fill=-1)
    with assert_raises(contains="len(queries) = 3, len(results) = 2"):
        index.locate_into(queries, too_short)


def test_knot_index_owns_storage_and_supports_span_windows() raises:
    var source: List[Float64] = [0.0, 1.0, 2.0]
    var index = KnotIndex(source.copy())
    source[1] = 99.0
    assert_equal(index.locate(1.5), 1)
    var queries: List[Float64] = [-10.0, 0.25, 1.0, 2.0, 10.0]
    var results = List[Int](length=5, fill=-1)
    index.locate_into(Span(queries)[1:4], Span(results)[1:4])
    assert_equal(results[0], -1)
    assert_equal(results[1], 0)
    assert_equal(results[2], 1)
    assert_equal(results[3], 1)
    assert_equal(results[4], -1)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
