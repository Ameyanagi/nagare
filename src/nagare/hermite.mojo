"""Shared piecewise-cubic Hermite interpolation machinery."""

from std.collections import List

from .extrapolation import ExtrapolationPolicy, _domain_error_message
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


struct CubicHermiteInterpolator(Copyable):
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
            raise Error("query must be finite")

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

    def derivative(self, x: Float64) raises -> Float64:
        """Evaluate the first derivative under the configured policy.

        Outside the domain, `CLAMP` returns zero, `LINEAR` returns the endpoint
        tangent slope, `FILL` returns its payload, and `ERROR` rejects the
        query. Non-finite queries always raise.
        """
        if not _is_finite(x):
            raise Error("query must be finite")

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

    def evaluate(self, queries: Span[Float64, _]) raises -> List[Float64]:
        """Allocate and return results for finite queries in order.

        Raises on the first offending query and returns no partial results.
        Empty input returns an empty list.
        """
        var results = List[Float64](length=len(queries), fill=0.0)
        self.evaluate_into(queries, results)
        return results^

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
