"""Validated binary-search interval location."""

from std.collections import List


def _is_finite(value: Float64) -> Bool:
    # NaN fails the first comparison. Either infinity produces NaN when
    # subtracted from itself and fails the second comparison.
    return value == value and value - value == 0.0


def _validate_knots(knots: List[Float64]) raises:
    if len(knots) < 2:
        raise Error("knot sequence must contain at least two values")
    if not _is_finite(knots[0]):
        raise Error("knots must be finite")
    for index in range(1, len(knots)):
        if not _is_finite(knots[index]):
            raise Error("knots must be finite")
        if knots[index] <= knots[index - 1]:
            raise Error("knots must be strictly increasing")


def _locate_interval_in_domain(knots: List[Float64], x: Float64) -> Int:
    """Locate a known in-domain finite query without validating the table."""
    if x == knots[len(knots) - 1]:
        return len(knots) - 2

    # upper_bound(x) - 1 makes exact interior knots right-biased.
    var low = 0
    var high = len(knots)
    while low < high:
        var middle = low + (high - low) // 2
        if knots[middle] <= x:
            low = middle + 1
        else:
            high = middle
    return low - 1


def locate_interval(knots: List[Float64], x: Float64) raises -> Int:
    """Return the left index of the interval containing `x`.

    Knots must be finite and strictly increasing. The query must be finite and
    inside the inclusive knot domain. Exact interior knots select the interval
    to their right; the final knot selects the final interval.

    Complexity is O(n) validation followed by O(log n) binary search.
    """
    _validate_knots(knots)
    if not _is_finite(x):
        raise Error("query must be finite")
    if x < knots[0] or x > knots[len(knots) - 1]:
        raise Error("query is outside the knot domain")
    return _locate_interval_in_domain(knots, x)
