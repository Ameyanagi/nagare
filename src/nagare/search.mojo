"""Validated binary-search interval location."""

from std.collections import List
from std.io import Writable, Writer


def _is_finite(value: Float64) -> Bool:
    # NaN fails the first comparison. Either infinity produces NaN when
    # subtracted from itself and fails the second comparison.
    return value == value and value - value == 0.0


def _validate_knots(knots: List[Float64]) raises:
    if len(knots) < 2:
        raise Error(
            String(
                "knot sequence must contain at least two values: received ",
                len(knots),
            )
        )
    if not _is_finite(knots[0]):
        raise Error(String("knots must be finite: knots[0] is ", knots[0]))
    for index in range(1, len(knots)):
        if not _is_finite(knots[index]):
            raise Error(
                String(
                    "knots must be finite: knots[",
                    index,
                    "] is ",
                    knots[index],
                )
            )
        if knots[index] <= knots[index - 1]:
            raise Error(
                String(
                    "knots must be strictly increasing: knots[",
                    index,
                    "] = ",
                    knots[index],
                    " <= knots[",
                    index - 1,
                    "] = ",
                    knots[index - 1],
                )
            )


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


def _validate_query(knots: List[Float64], x: Float64) raises:
    if not _is_finite(x):
        raise Error(String("query must be finite: received ", x))
    if x < knots[0] or x > knots[len(knots) - 1]:
        raise Error(
            String(
                "query ",
                x,
                " is outside the knot domain [",
                knots[0],
                ", ",
                knots[len(knots) - 1],
                "]",
            )
        )


def locate_interval(knots: List[Float64], x: Float64) raises -> Int:
    """Return the left index of the interval containing `x`.

    Knots must be finite and strictly increasing. The query must be finite and
    inside the inclusive knot domain. Exact interior knots select the interval
    to their right; the final knot selects the final interval.

    Complexity is O(n) validation followed by O(log n) binary search. Use an
    owned `KnotIndex` to validate once and reuse O(log n) query lookup.
    """
    _validate_knots(knots)
    _validate_query(knots, x)
    return _locate_interval_in_domain(knots, x)


struct KnotIndex(Copyable, Equatable, Writable):
    """Owned validated knots for custom interval operations.

    Construction takes ownership and validates finite, strictly increasing
    knots in O(n). Lookup validates only the new query and binary-searches in
    O(log n), without allocation or rescanning knots. Direct mutation of
    underscore-prefixed storage is outside the contract; `validate()` provides
    an explicit checkpoint after unusual access.
    """

    var _knots: List[Float64]

    def __init__(out self, var knots: List[Float64]) raises:
        """Take ownership of at least two finite, strictly increasing knots."""
        _validate_knots(knots)
        self._knots = knots^

    def validate(self) raises:
        """Explicitly recheck stored knot finiteness and strict ordering."""
        _validate_knots(self._knots)

    def knot_count(self) -> Int:
        """Return the number of knots without revalidating storage."""
        return len(self._knots)

    def knots(self) -> Span[Float64, origin_of(self._knots)]:
        """Return a read-only view without copying or revalidating knots."""
        return self._knots

    def locate(self, x: Float64) raises -> Int:
        """Return the interval for a finite query in the inclusive domain.

        Exact interior knots select the interval to their right; the final
        knot selects the final interval, matching `locate_interval`.
        """
        _validate_query(self._knots, x)
        return _locate_interval_in_domain(self._knots, x)

    def locate_into(
        self,
        queries: Span[Float64, _],
        results: Span[mut=True, Int, _],
    ) raises:
        """Locate unsorted queries into caller-owned storage without allocation.

        Buffer lengths must match. Empty input is valid. Query validation and
        O(log n) search run once per query; no knot validation is repeated.
        Results are unspecified after an invalid query raises, as with the
        interpolators' `evaluate_into` APIs.
        """
        if len(queries) != len(results):
            raise Error(
                String(
                    "query and result buffers must have equal length: ",
                    "len(queries) = ",
                    len(queries),
                    ", len(results) = ",
                    len(results),
                )
            )
        for index in range(len(queries)):
            results[index] = self.locate(queries[index])

    def __eq__(self, other: Self) -> Bool:
        """Compare the validated knot sequences exactly."""
        if len(self._knots) != len(other._knots):
            return False
        for index in range(len(self._knots)):
            if self._knots[index] != other._knots[index]:
                return False
        return True

    def write_to[W: Writer](self, mut writer: W):
        """Write the knot count and inclusive domain."""
        writer.write(
            "KnotIndex(knots=",
            self.knot_count(),
            ", domain=[",
            self._knots[0],
            ", ",
            self._knots[len(self._knots) - 1],
            "])",
        )
