"""Shared piecewise-cubic Hermite interpolation machinery."""

from std.collections import List
from std.io import Writable, Writer

from .extrapolation import (
    ExtrapolationPolicy,
    _domain_error_message,
    _integrate_exterior_tail,
    _validate_integration_bounds,
)
from .linear import _segment_parameter, _validate_table
from .search import _is_finite, _locate_interval_in_domain


def _validate_hermite_table(
    knots: List[Float64],
    values: List[Float64],
    slopes: List[Float64],
) raises:
    """Validate the knot, value, and knot-slope representation."""
    _validate_table(knots, values)
    if len(slopes) != len(knots):
        raise Error(
            String(
                "knot, value, and slope sequences must have equal length: ",
                "len(knots) = ",
                len(knots),
                ", len(values) = ",
                len(values),
                ", len(slopes) = ",
                len(slopes),
            )
        )
    for index in range(len(slopes)):
        if not _is_finite(slopes[index]):
            raise Error(
                String(
                    "interpolation slopes must be finite: slopes[",
                    index,
                    "] is ",
                    slopes[index],
                )
            )


def _evaluate_segment(
    knots: List[Float64],
    values: List[Float64],
    slopes: List[Float64],
    index: Int,
    x: Float64,
) -> Float64:
    """Evaluate one interval in cubic Hermite basis form."""
    var x0 = knots[index]
    var x1 = knots[index + 1]
    if x == x0:
        return values[index]
    if x == x1:
        return values[index + 1]

    var width = x1 - x0
    var parameter = _segment_parameter(x0, x1, x)
    var parameter_squared = parameter * parameter
    var parameter_cubed = parameter_squared * parameter
    var left_value_basis = 2.0 * parameter_cubed - 3.0 * parameter_squared + 1.0
    var left_slope_basis = parameter_cubed - 2.0 * parameter_squared + parameter
    var right_value_basis = -2.0 * parameter_cubed + 3.0 * parameter_squared
    var right_slope_basis = parameter_cubed - parameter_squared
    return (
        left_value_basis * values[index]
        + left_slope_basis * width * slopes[index]
        + right_value_basis * values[index + 1]
        + right_slope_basis * width * slopes[index + 1]
    )


def _derivative_segment(
    knots: List[Float64],
    values: List[Float64],
    slopes: List[Float64],
    index: Int,
    x: Float64,
) -> Float64:
    """Evaluate the first derivative of one cubic Hermite interval."""
    var x0 = knots[index]
    var x1 = knots[index + 1]
    if x == x0:
        return slopes[index]
    if x == x1:
        return slopes[index + 1]

    var width = x1 - x0
    var parameter = _segment_parameter(x0, x1, x)
    var parameter_squared = parameter * parameter
    return (
        (6.0 * parameter_squared - 6.0 * parameter)
        * (values[index] - values[index + 1])
        / width
        + (3.0 * parameter_squared - 4.0 * parameter + 1.0) * slopes[index]
        + (3.0 * parameter_squared - 2.0 * parameter) * slopes[index + 1]
    )


def _second_derivative_segment(
    knots: List[Float64],
    values: List[Float64],
    slopes: List[Float64],
    index: Int,
    x: Float64,
) -> Float64:
    """Evaluate the second derivative of one cubic Hermite interval."""
    var x0 = knots[index]
    var x1 = knots[index + 1]
    var width = x1 - x0
    var parameter = _segment_parameter(x0, x1, x)
    return (
        (12.0 * parameter - 6.0) * (values[index] - values[index + 1]) / width
        + (6.0 * parameter - 4.0) * slopes[index]
        + (6.0 * parameter - 2.0) * slopes[index + 1]
    ) / width


def _hermite_antiderivative(
    parameter: Float64,
    width: Float64,
    left_value: Float64,
    right_value: Float64,
    left_slope: Float64,
    right_slope: Float64,
) -> Float64:
    """Evaluate an interval's Hermite-basis antiderivative from its left end."""
    var parameter_squared = parameter * parameter
    var parameter_cubed = parameter_squared * parameter
    var parameter_fourth = parameter_cubed * parameter
    var left_value_basis = 0.5 * parameter_fourth - parameter_cubed + parameter
    var left_slope_basis = (
        0.25 * parameter_fourth
        - (2.0 / 3.0) * parameter_cubed
        + 0.5 * parameter_squared
    )
    var right_value_basis = -0.5 * parameter_fourth + parameter_cubed
    var right_slope_basis = 0.25 * parameter_fourth - (1.0 / 3.0) * parameter_cubed
    return width * (
        left_value_basis * left_value
        + left_slope_basis * width * left_slope
        + right_value_basis * right_value
        + right_slope_basis * width * right_slope
    )


def _integrate_segment(
    knots: List[Float64],
    values: List[Float64],
    slopes: List[Float64],
    index: Int,
    start: Float64,
    end: Float64,
) -> Float64:
    """Integrate one partial cubic Hermite interval exactly."""
    var x0 = knots[index]
    var x1 = knots[index + 1]
    var width = x1 - x0
    var start_parameter = _segment_parameter(x0, x1, start)
    var end_parameter = _segment_parameter(x0, x1, end)
    return _hermite_antiderivative(
        end_parameter,
        width,
        values[index],
        values[index + 1],
        slopes[index],
        slopes[index + 1],
    ) - _hermite_antiderivative(
        start_parameter,
        width,
        values[index],
        values[index + 1],
        slopes[index],
        slopes[index + 1],
    )


def _left_slope(slopes: List[Float64]) -> Float64:
    """Return the left endpoint tangent slope."""
    return slopes[0]


def _right_slope(slopes: List[Float64]) -> Float64:
    """Return the right endpoint tangent slope."""
    return slopes[len(slopes) - 1]


def _evaluate_linear_ray(
    endpoint_x: Float64,
    endpoint_y: Float64,
    endpoint_slope: Float64,
    x: Float64,
) -> Float64:
    """Evaluate the tangent ray rooted at one endpoint."""
    if endpoint_slope == 0.0:
        return endpoint_y
    return endpoint_y + endpoint_slope * (x - endpoint_x)


struct CubicHermiteInterpolator(Copyable, Equatable, Writable):
    """An owning C1 piecewise-cubic Hermite interpolant.

    Each supplied slope is `dy/dx` at the corresponding knot. Construction
    validates the knot, value, and slope table in O(n); each finite query takes
    O(log n) for interval search and O(1) segment evaluation. Exact knot queries
    reproduce the stored values, and knot derivative queries reproduce the
    supplied slopes.

    Construction establishes the representation invariant, which read-only
    methods trust thereafter. Underscore-prefixed fields are private by
    convention; mutating them directly is outside the contract. Call
    `validate()` for an explicit checkpoint after unusual direct access.

    Outside the domain, `LINEAR` follows the endpoint tangent ray, preserving
    C1 continuity at the join. `CLAMP` differentiates to zero, and `FILL`
    returns its payload for both value and derivative queries.
    """

    var _knots: List[Float64]
    var _values: List[Float64]
    var _slopes: List[Float64]
    var _extrapolation: ExtrapolationPolicy

    def __init__(
        out self,
        var knots: List[Float64],
        var values: List[Float64],
        var slopes: List[Float64],
        extrapolation: ExtrapolationPolicy = ExtrapolationPolicy.ERROR,
    ) raises:
        """Validate and take ownership of values and `dy/dx` knot slopes."""
        _validate_hermite_table(knots, values, slopes)

        self._knots = knots^
        self._values = values^
        self._slopes = slopes^
        self._extrapolation = extrapolation

    def validate(self) raises:
        """Explicitly revalidate the stored knot, value, and slope table."""
        _validate_hermite_table(self._knots, self._values, self._slopes)

    def knots(self) -> Span[Float64, origin_of(self._knots)]:
        """Read-only view of the validated knot sequence."""
        return Span(self._knots)

    def values(self) -> Span[Float64, origin_of(self._values)]:
        """Read-only view of the validated value sequence."""
        return Span(self._values)

    def slopes(self) -> Span[Float64, origin_of(self._slopes)]:
        """Read-only view of the validated slope sequence."""
        return Span(self._slopes)

    def knot_count(self) -> Int:
        """Return the knot count."""
        return len(self._knots)

    def domain_start(self) -> Float64:
        """Return the lower domain bound."""
        return self._knots[0]

    def domain_end(self) -> Float64:
        """Return the upper domain bound."""
        return self._knots[len(self._knots) - 1]

    def evaluate(self, x: Float64) raises -> Float64:
        """Evaluate one finite query under the configured policy.

        Exact knot queries return the stored ordinate. `LINEAR` extrapolation
        uses the appropriate endpoint tangent ray, while `FILL` returns its
        payload outside the domain.
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
            return _evaluate_linear_ray(
                self._knots[0], self._values[0], _left_slope(self._slopes), x
            )

        var final_knot = len(self._knots) - 1
        if x > self._knots[final_knot]:
            if self._extrapolation == ExtrapolationPolicy.ERROR:
                raise Error(
                    _domain_error_message(x, self._knots[0], self._knots[final_knot])
                )
            if self._extrapolation == ExtrapolationPolicy.CLAMP:
                return self._values[final_knot]
            if self._extrapolation.is_fill():
                return self._extrapolation.fill_value()
            return _evaluate_linear_ray(
                self._knots[final_knot],
                self._values[final_knot],
                _right_slope(self._slopes),
                x,
            )

        return _evaluate_segment(
            self._knots,
            self._values,
            self._slopes,
            _locate_interval_in_domain(self._knots, x),
            x,
        )

    def __call__(self, x: Float64) raises -> Float64:
        """Call `evaluate`; `evaluate` is the primary documented name."""
        return self.evaluate(x)

    def derivative(self, x: Float64) raises -> Float64:
        """Evaluate the first derivative under the configured policy.

        Outside the domain, `CLAMP` returns zero, `LINEAR` returns the endpoint
        tangent slope, `FILL` returns its payload, and `ERROR` rejects the
        query. Non-finite queries always raise.
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
            return _left_slope(self._slopes)

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
            return _right_slope(self._slopes)

        return _derivative_segment(
            self._knots,
            self._values,
            self._slopes,
            _locate_interval_in_domain(self._knots, x),
            x,
        )

    def derivative(self, queries: Span[Float64, _]) raises -> List[Float64]:
        """Allocate and return first derivatives for finite queries in order.

        Raises on the first offending query and returns no partial results.
        Empty input returns an empty list.
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

    def second_derivative(self, x: Float64) raises -> Float64:
        """Evaluate the second derivative under the configured policy.

        Interior knots select the segment to their right, matching evaluation's
        interval search; the final knot selects the final segment. Outside the
        domain, `CLAMP` returns zero, the affine `LINEAR` tangent-ray extension
        also returns zero, `FILL` returns its payload, and `ERROR` rejects the
        query. Non-finite queries always raise.
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
            if self._extrapolation.is_fill():
                return self._extrapolation.fill_value()
            return 0.0

        var final_knot = len(self._knots) - 1
        if x > self._knots[final_knot]:
            if self._extrapolation == ExtrapolationPolicy.ERROR:
                raise Error(
                    _domain_error_message(x, self._knots[0], self._knots[final_knot])
                )
            if self._extrapolation.is_fill():
                return self._extrapolation.fill_value()
            return 0.0

        return _second_derivative_segment(
            self._knots,
            self._values,
            self._slopes,
            _locate_interval_in_domain(self._knots, x),
            x,
        )

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
                    _left_slope(self._slopes),
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
                    _right_slope(self._slopes),
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
            result += _integrate_segment(
                self._knots,
                self._values,
                self._slopes,
                index,
                position,
                segment_end,
            )
            position = segment_end
            index += 1
        return result

    def evaluate(self, queries: Span[Float64, _]) raises -> List[Float64]:
        """Allocate and return results for finite queries in order.

        Raises on the first offending query and returns no partial results.
        Empty input returns an empty list.
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

    def __eq__(self, other: Self) -> Bool:
        """Compare validated tables, stored slopes, and policy exactly.

        Validation excludes NaN from every table, so elementwise `Float64`
        equality is sound.
        """
        if (
            len(self._knots) != len(other._knots)
            or len(self._values) != len(other._values)
            or len(self._slopes) != len(other._slopes)
            or self._extrapolation != other._extrapolation
        ):
            return False
        for index in range(len(self._knots)):
            if (
                self._knots[index] != other._knots[index]
                or self._values[index] != other._values[index]
                or self._slopes[index] != other._slopes[index]
            ):
                return False
        return True

    def __str__(self) -> String:
        var result = String()
        self.write_to(result)
        return result^

    def write_to[W: Writer](self, mut writer: W):
        writer.write(
            "CubicHermiteInterpolator(",
            len(self._knots),
            " knots on [",
            self._knots[0],
            ", ",
            self._knots[len(self._knots) - 1],
            "], extrapolation=",
            self._extrapolation,
            ")",
        )
