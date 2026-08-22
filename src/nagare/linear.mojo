"""Validated one-dimensional piecewise-linear interpolation."""

from std.collections import List
from std.io import Writable, Writer

from .extrapolation import (
    ExtrapolationPolicy,
    _domain_error_message,
    _integrate_exterior_tail,
    _validate_integration_bounds,
)
from .search import _is_finite, _locate_interval_in_domain, _validate_knots


def _segment_parameter(x0: Float64, x1: Float64, x: Float64) -> Float64:
    """Return `(x - x0) / (x1 - x0)` without avoidable overflow."""
    var numerator = x - x0
    var denominator = x1 - x0
    if _is_finite(numerator) and _is_finite(denominator):
        return numerator / denominator

    # A difference of opposite-sign finite values can overflow even when its
    # ratio is small. Scaling all coordinates by the same positive value keeps
    # the ratio unchanged while bounding both differences by two.
    var scale = max(abs(x), max(abs(x0), abs(x1)))
    var scaled_denominator = x1 / scale - x0 / scale
    if scaled_denominator != 0.0:
        return (x / scale - x0 / scale) / scaled_denominator

    # This case requires an extrapolated parameter beyond Float64 resolution:
    # the finite knots became indistinguishable only after scaling by a much
    # larger query. Dividing before subtracting produces the required signed
    # infinity without an infinity-minus-infinity operation.
    return x / denominator - x0 / denominator


def _stable_linear_value(y0: Float64, y1: Float64, parameter: Float64) -> Float64:
    """Evaluate an affine combination through the safest representable form."""
    if y0 == y1:
        return y0

    # A directly representable ordinate delta is the most accurate path and is
    # essential when a huge extrapolation parameter multiplies tiny ordinates.
    # Accept it only if every intermediate and the result remain finite; a
    # scaled fallback can recover cancellation cases that overflow this form.
    var delta = y1 - y0
    if _is_finite(delta):
        var offset = parameter * delta
        if _is_finite(offset):
            var direct = y0 + offset
            if _is_finite(direct):
                return direct

    var scale = max(abs(y0), abs(y1))
    if scale == 0.0:
        return 0.0

    # Forming y1 - y0 can overflow for opposite-sign finite endpoints. Scaling
    # first keeps the delta in [-2, 2]. The final multiplication overflows to a
    # signed infinity exactly when the represented affine result exceeds the
    # finite Float64 range.
    var normalized_y0 = y0 / scale
    var normalized_y1 = y1 / scale
    var normalized = normalized_y0 + parameter * (normalized_y1 - normalized_y0)
    return scale * normalized


def _stable_secant(
    x0: Float64,
    x1: Float64,
    y0: Float64,
    y1: Float64,
) -> Float64:
    """Return the segment slope without an avoidable infinity/infinity."""
    var delta_x = x1 - x0
    var delta_y = y1 - y0
    if _is_finite(delta_x) and _is_finite(delta_y):
        return delta_y / delta_x

    var y_scale = max(abs(y0), abs(y1))
    if y_scale == 0.0:
        return 0.0
    var x_scale = max(abs(x0), abs(x1))
    var normalized_x = x1 / x_scale - x0 / x_scale
    var normalized_y = y1 / y_scale - y0 / y_scale
    return (normalized_y / normalized_x) * (y_scale / x_scale)


def _validate_table(knots: List[Float64], values: List[Float64]) raises:
    """Validate every invariant required by a linear interpolation table."""
    _validate_knots(knots)
    if len(values) != len(knots):
        raise Error(
            String(
                "knot and value sequences must have equal length: len(knots) = ",
                len(knots),
                ", len(values) = ",
                len(values),
            )
        )
    for index in range(len(values)):
        if not _is_finite(values[index]):
            raise Error(
                String(
                    "interpolation values must be finite: values[",
                    index,
                    "] is ",
                    values[index],
                )
            )


def _validate_sorted_queries(queries: Span[Float64, _]) raises:
    """Validate finite, nondecreasing queries before a sorted batch write."""
    for index in range(len(queries)):
        var query = queries[index]
        if not _is_finite(query):
            raise Error(String("query must be finite: received ", query))
        if index > 0 and query < queries[index - 1]:
            raise Error(
                String(
                    "queries must be sorted in nondecreasing order: queries[",
                    index,
                    "] = ",
                    query,
                    " < queries[",
                    index - 1,
                    "] = ",
                    queries[index - 1],
                )
            )


struct LinearInterpolator(Copyable, Equatable, Writable):
    """An owning piecewise-linear interpolant over finite `Float64` data.

    Construction validates the table, which read-only methods trust thereafter.
    Underscore-prefixed fields are private by convention; mutating them directly
    is outside the contract. Call `validate()` for an explicit checkpoint after
    unusual direct access.
    """

    var _knots: List[Float64]
    var _values: List[Float64]
    var _extrapolation: ExtrapolationPolicy

    def __init__(
        out self,
        var knots: List[Float64],
        var values: List[Float64],
        extrapolation: ExtrapolationPolicy = ExtrapolationPolicy.ERROR,
    ) raises:
        """Validate and take ownership of a knot/value table."""
        _validate_table(knots, values)

        self._knots = knots^
        self._values = values^
        self._extrapolation = extrapolation

    def validate(self) raises:
        """Explicitly revalidate the stored knot and value table."""
        _validate_table(self._knots, self._values)

    def knots(self) -> Span[Float64, origin_of(self._knots)]:
        """Read-only view of the validated knot sequence."""
        return Span(self._knots)

    def values(self) -> Span[Float64, origin_of(self._values)]:
        """Read-only view of the validated value sequence."""
        return Span(self._values)

    def knot_count(self) -> Int:
        """Return the knot count."""
        return len(self._knots)

    def domain_start(self) -> Float64:
        """Return the lower domain bound."""
        return self._knots[0]

    def domain_end(self) -> Float64:
        """Return the upper domain bound."""
        return self._knots[len(self._knots) - 1]

    def _evaluate_segment(self, index: Int, x: Float64) -> Float64:
        var x0 = self._knots[index]
        var x1 = self._knots[index + 1]
        var y0 = self._values[index]
        var y1 = self._values[index + 1]

        # Preserve every tabulated value exactly rather than allowing a
        # multiply/add sequence to perturb it by one rounding step.
        if x == x0:
            return y0
        if x == x1:
            return y1
        return _stable_linear_value(y0, y1, _segment_parameter(x0, x1, x))

    def _segment_slope(self, index: Int) -> Float64:
        return _stable_secant(
            self._knots[index],
            self._knots[index + 1],
            self._values[index],
            self._values[index + 1],
        )

    def _integrate_segment(self, index: Int, start: Float64, end: Float64) -> Float64:
        var left = self._evaluate_segment(index, start)
        var right = self._evaluate_segment(index, end)
        return (0.5 * left + 0.5 * right) * (end - start)

    def evaluate(self, x: Float64) raises -> Float64:
        """Evaluate one finite query under the configured extrapolation policy.

        In-domain results avoid intermediate overflow for finite endpoint data.
        `LINEAR` extrapolation returns signed infinity when its represented
        result exceeds the finite `Float64` range. `FILL` returns its payload
        outside the knot domain.
        """
        if not _is_finite(x):
            raise Error(String("query must be finite: received ", x))

        if x < self._knots[0]:
            if self._extrapolation == ExtrapolationPolicy.ERROR:
                raise Error(
                    _domain_error_message(
                        x, self._knots[0], self._knots[len(self._knots) - 1]
                    )
                )
            if self._extrapolation == ExtrapolationPolicy.CLAMP:
                return self._values[0]
            if self._extrapolation.is_fill():
                return self._extrapolation.fill_value()
            return self._evaluate_segment(0, x)

        var final_index = len(self._knots) - 1
        if x > self._knots[final_index]:
            if self._extrapolation == ExtrapolationPolicy.ERROR:
                raise Error(
                    _domain_error_message(x, self._knots[0], self._knots[final_index])
                )
            if self._extrapolation == ExtrapolationPolicy.CLAMP:
                return self._values[final_index]
            if self._extrapolation.is_fill():
                return self._extrapolation.fill_value()
            return self._evaluate_segment(final_index - 1, x)

        return self._evaluate_segment(_locate_interval_in_domain(self._knots, x), x)

    def __call__(self, x: Float64) raises -> Float64:
        """Call `evaluate`; `evaluate` is the primary documented name."""
        return self.evaluate(x)

    def derivative(self, x: Float64) raises -> Float64:
        """Return the piecewise-constant first derivative.

        Interior knots select the segment to their right, matching evaluation's
        interval search. The final knot selects the final segment. Outside the
        domain, `CLAMP` returns zero, `LINEAR` returns the endpoint segment
        slope, `FILL` returns its payload, and `ERROR` rejects the query.
        """
        if not _is_finite(x):
            raise Error(String("query must be finite: received ", x))

        if x < self._knots[0]:
            if self._extrapolation == ExtrapolationPolicy.ERROR:
                raise Error(
                    _domain_error_message(
                        x, self._knots[0], self._knots[len(self._knots) - 1]
                    )
                )
            if self._extrapolation == ExtrapolationPolicy.CLAMP:
                return 0.0
            if self._extrapolation.is_fill():
                return self._extrapolation.fill_value()
            return self._segment_slope(0)

        var final_knot = len(self._knots) - 1
        if x > self._knots[final_knot]:
            if self._extrapolation == ExtrapolationPolicy.ERROR:
                raise Error(
                    _domain_error_message(x, self._knots[0], self._knots[final_knot])
                )
            if self._extrapolation == ExtrapolationPolicy.CLAMP:
                return 0.0
            if self._extrapolation.is_fill():
                return self._extrapolation.fill_value()
            return self._segment_slope(final_knot - 1)

        return self._segment_slope(_locate_interval_in_domain(self._knots, x))

    def derivative(self, queries: Span[Float64, _]) raises -> List[Float64]:
        """Allocate and return first derivatives for finite queries in order.

        Raises on the first offending query (a non-finite query under every
        policy, or an out-of-domain query under `ERROR`) and returns no partial
        results. Empty input returns an empty list.
        """
        var results = List[Float64](length=len(queries), fill=0.0)
        self.derivative_into(queries, results)
        return results^

    def derivative_into(
        self,
        queries: Span[Float64, _],
        results: Span[mut=True, Float64, _],
    ) raises:
        """Evaluate each first derivative into `results` without allocating.

        Raises if the buffer lengths differ or on the first offending query;
        results contents are unspecified after a raise.
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
            results[index] = self.derivative(queries[index])

    def integrate(self, a: Float64, b: Float64) raises -> Float64:
        """Definite integral over `[a, b]` (sign-flipped when `a > b`).

        Segments are integrated in closed form. Portions outside the domain
        follow the configured extrapolation policy; `ERROR` rejects them, and
        `FILL` contributes `fill_value * width` (NaN-propagating by default).
        """
        _validate_integration_bounds(a, b)
        if a == b:
            return 0.0
        if a > b:
            return -self.integrate(b, a)

        var domain_start = self._knots[0]
        var final_knot = len(self._knots) - 1
        var domain_end = self._knots[final_knot]
        var result = 0.0
        if a < domain_start:
            var tail_end = min(b, domain_start)
            if a < tail_end:
                result += _integrate_exterior_tail(
                    self._extrapolation,
                    a,
                    tail_end,
                    domain_start,
                    self._values[0],
                    self._segment_slope(0),
                    domain_start,
                    domain_end,
                )
        if b > domain_end:
            var tail_start = max(a, domain_end)
            if tail_start < b:
                result += _integrate_exterior_tail(
                    self._extrapolation,
                    tail_start,
                    b,
                    domain_end,
                    self._values[final_knot],
                    self._segment_slope(final_knot - 1),
                    domain_start,
                    domain_end,
                )

        var interior_start = max(a, domain_start)
        var interior_end = min(b, domain_end)
        if interior_start >= interior_end:
            return result

        var index = _locate_interval_in_domain(self._knots, interior_start)
        var position = interior_start
        while position < interior_end:
            var segment_end = min(interior_end, self._knots[index + 1])
            result += self._integrate_segment(index, position, segment_end)
            position = segment_end
            index += 1
        return result

    def evaluate(self, queries: Span[Float64, _]) raises -> List[Float64]:
        """Allocate and return results for finite queries in order.

        Raises on the first offending query (a non-finite query under every
        policy, or an out-of-domain query under `ERROR`) and returns no partial
        results. Empty input returns an empty list.
        """
        var results = List[Float64](length=len(queries), fill=0.0)
        self.evaluate_into(queries, results)
        return results^

    def __call__(self, queries: Span[Float64, _]) raises -> List[Float64]:
        """Call `evaluate`; `evaluate` is the primary documented name."""
        return self.evaluate(queries)

    def evaluate_into(
        self,
        queries: Span[Float64, _],
        results: Span[mut=True, Float64, _],
    ) raises:
        """Evaluate `queries[i]` into `results[i]` without allocating.

        Raises if the buffer lengths differ or on the first offending query;
        results contents are unspecified after a raise.
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
            results[index] = self.evaluate(queries[index])

    def evaluate_sorted(self, queries: Span[Float64, _]) raises -> List[Float64]:
        """Allocate and evaluate finite, nondecreasing queries monotonically.

        The result is numerically identical to scalar `evaluate` for every
        query. Unlike the order-agnostic batch API, this method exploits the
        monotone query contract to binary-search the first interior interval,
        then traverse only the remaining spanned intervals. Duplicate queries
        and exact knots are supported.
        """
        var results = List[Float64](length=len(queries), fill=0.0)
        self.evaluate_sorted_into(queries, results)
        return results^

    def evaluate_sorted_into(
        self,
        queries: Span[Float64, _],
        results: Span[mut=True, Float64, _],
    ) raises:
        """Evaluate sorted queries into caller-owned storage without allocation.

        Queries must be finite and sorted in nondecreasing order. The complete
        query sequence and ERROR-policy domain bounds are validated before any
        output is written. Extrapolated prefixes and suffixes follow the
        interpolator's configured policy.
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
        _validate_sorted_queries(queries)
        if len(queries) == 0:
            return

        var domain_start = self._knots[0]
        var final_knot = len(self._knots) - 1
        var domain_end = self._knots[final_knot]
        var query_index = 0

        # ERROR is transactional for the caller-owned buffer. Sortedness makes
        # the first query sufficient for the lower bound. Find the first upper
        # offender before writing so its diagnostic still matches scalar order.
        if self._extrapolation == ExtrapolationPolicy.ERROR:
            if queries[0] < domain_start:
                raise Error(_domain_error_message(queries[0], domain_start, domain_end))
            if queries[len(queries) - 1] > domain_end:
                var lower = 0
                var upper = len(queries)
                while lower < upper:
                    var middle = lower + (upper - lower) // 2
                    if queries[middle] <= domain_end:
                        lower = middle + 1
                    else:
                        upper = middle
                raise Error(
                    _domain_error_message(queries[lower], domain_start, domain_end)
                )

        # A sorted exterior prefix uses one endpoint policy and no interval
        # searches. ERROR still rejects the first offending query.
        while query_index < len(queries) and queries[query_index] < domain_start:
            var query = queries[query_index]
            if self._extrapolation == ExtrapolationPolicy.ERROR:
                raise Error(_domain_error_message(query, domain_start, domain_end))
            if self._extrapolation == ExtrapolationPolicy.CLAMP:
                results[query_index] = self._values[0]
            elif self._extrapolation.is_fill():
                results[query_index] = self._extrapolation.fill_value()
            else:
                results[query_index] = self._evaluate_segment(0, query)
            query_index += 1

        # Exact interior knots remain right-biased, except the final knot,
        # matching `_locate_interval_in_domain` and scalar evaluation exactly.
        var final_segment = final_knot - 1
        var segment = 0
        if query_index < len(queries) and queries[query_index] <= domain_end:
            segment = _locate_interval_in_domain(self._knots, queries[query_index])
        while query_index < len(queries) and queries[query_index] <= domain_end:
            var query = queries[query_index]
            while segment < final_segment and query >= self._knots[segment + 1]:
                segment += 1
            results[query_index] = self._evaluate_segment(segment, query)
            query_index += 1

        # Any remaining sorted suffix lies above the domain.
        while query_index < len(queries):
            var query = queries[query_index]
            if self._extrapolation == ExtrapolationPolicy.ERROR:
                raise Error(_domain_error_message(query, domain_start, domain_end))
            if self._extrapolation == ExtrapolationPolicy.CLAMP:
                results[query_index] = self._values[final_knot]
            elif self._extrapolation.is_fill():
                results[query_index] = self._extrapolation.fill_value()
            else:
                results[query_index] = self._evaluate_segment(final_segment, query)
            query_index += 1

    def __eq__(self, other: Self) -> Bool:
        """Compare validated tables and policy exactly.

        Validated knot and value tables never contain NaN, so elementwise
        `Float64` equality is sound.
        """
        if (
            len(self._knots) != len(other._knots)
            or len(self._values) != len(other._values)
            or self._extrapolation != other._extrapolation
        ):
            return False
        for index in range(len(self._knots)):
            if (
                self._knots[index] != other._knots[index]
                or self._values[index] != other._values[index]
            ):
                return False
        return True

    def __str__(self) -> String:
        var result = String()
        self.write_to(result)
        return result^

    def write_to[W: Writer](self, mut writer: W):
        writer.write(
            "LinearInterpolator(",
            len(self._knots),
            " knots on [",
            self._knots[0],
            ", ",
            self._knots[len(self._knots) - 1],
            "], extrapolation=",
            self._extrapolation,
            ")",
        )
