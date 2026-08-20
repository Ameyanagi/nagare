"""Validated one-dimensional natural cubic spline interpolation."""

from std.collections import List

from .extrapolation import ExtrapolationPolicy, _domain_error_message
from .linear import _validate_table
from .search import _is_finite, _locate_interval_in_domain


def _solve_natural_moments(
    knots: List[Float64], values: List[Float64]
) -> List[Float64]:
    """Return knot second derivatives from the natural tridiagonal system."""
    var knot_count = len(knots)
    var moments = List[Float64](capacity=knot_count)
    for _ in range(knot_count):
        moments.append(0.0)

    var interior_count = knot_count - 2
    if interior_count == 0:
        return moments^

    # Thomas elimination over only the interior unknowns. The endpoint moments
    # are already exactly zero. Strictly increasing knots make this natural
    # spline system strictly diagonally dominant, so no singularity path exists.
    var modified_upper = List[Float64](capacity=interior_count)
    var modified_rhs = List[Float64](capacity=interior_count)
    for knot_index in range(1, knot_count - 1):
        var left_width = knots[knot_index] - knots[knot_index - 1]
        var right_width = knots[knot_index + 1] - knots[knot_index]
        var diagonal = 2.0 * (left_width + right_width)
        var rhs = 6.0 * (
            (values[knot_index + 1] - values[knot_index]) / right_width
            - (values[knot_index] - values[knot_index - 1]) / left_width
        )

        var denominator = diagonal
        var reduced_rhs = rhs
        if knot_index > 1:
            var previous = knot_index - 2
            denominator -= left_width * modified_upper[previous]
            reduced_rhs -= left_width * modified_rhs[previous]

        modified_upper.append(right_width / denominator)
        modified_rhs.append(reduced_rhs / denominator)

    var unknown_index = interior_count - 1
    moments[unknown_index + 1] = modified_rhs[unknown_index]
    while unknown_index > 0:
        unknown_index -= 1
        moments[unknown_index + 1] = (
            modified_rhs[unknown_index]
            - modified_upper[unknown_index] * moments[unknown_index + 2]
        )
    return moments^


def _validate_coefficient_buffers(
    a: List[Float64],
    b: List[Float64],
    c: List[Float64],
    d: List[Float64],
    interval_count: Int,
) raises:
    if (
        len(a) != interval_count
        or len(b) != interval_count
        or len(c) != interval_count
        or len(d) != interval_count
    ):
        raise Error("spline coefficient buffers must match interval count")

    for index in range(interval_count):
        if (
            not _is_finite(a[index])
            or not _is_finite(b[index])
            or not _is_finite(c[index])
            or not _is_finite(d[index])
        ):
            raise Error("spline coefficients must be finite")


struct CubicSplineInterpolator(Copyable):
    """An owning natural cubic spline over finite `Float64` data.

    Each interval stores contiguous `a`, `b`, `c`, and `d` buffers for
    `S_i(dx) = a_i + b_i dx + c_i dx^2 + d_i dx^3`, where
    `dx = x - knots[i]`. Construction solves the natural tridiagonal system,
    so the endpoint second derivatives are exactly zero.

    Construction validates the table and coefficient representation, which
    read-only methods trust thereafter. Underscore-prefixed fields are private
    by convention; mutating them directly is outside the contract. Call
    `validate()` for an explicit checkpoint after unusual direct access.

    Outside the domain, `LINEAR` follows the endpoint tangent ray rather than
    extending a boundary cubic. Because a natural spline has zero endpoint
    second derivative, that ray is C2-continuous with the spline at the join.
    `FILL` returns its payload for values and both derivative orders outside the
    domain.
    """

    var _knots: List[Float64]
    var _values: List[Float64]
    var _a: List[Float64]
    var _b: List[Float64]
    var _c: List[Float64]
    var _d: List[Float64]
    var _extrapolation: ExtrapolationPolicy

    def __init__(
        out self,
        var knots: List[Float64],
        var values: List[Float64],
        extrapolation: ExtrapolationPolicy = ExtrapolationPolicy.ERROR,
    ) raises:
        """Validate, solve, and take ownership of a knot/value table."""
        _validate_table(knots, values)

        var moments = _solve_natural_moments(knots, values)
        var interval_count = len(knots) - 1
        var a = List[Float64](capacity=interval_count)
        var b = List[Float64](capacity=interval_count)
        var c = List[Float64](capacity=interval_count)
        var d = List[Float64](capacity=interval_count)
        for index in range(interval_count):
            var width = knots[index + 1] - knots[index]
            var left_moment = moments[index]
            var right_moment = moments[index + 1]
            a.append(values[index])
            b.append(
                (values[index + 1] - values[index]) / width
                - width * (2.0 * left_moment + right_moment) / 6.0
            )
            c.append(left_moment / 2.0)
            d.append((right_moment - left_moment) / (6.0 * width))

        # Non-finite coefficients mean finite inputs exceeded the usable range
        # of this Float64 representation. This is a genuine construction
        # failure, separate from the mathematically nonsingular Thomas solve.
        _validate_coefficient_buffers(a, b, c, d, interval_count)

        self._knots = knots^
        self._values = values^
        self._a = a^
        self._b = b^
        self._c = c^
        self._d = d^
        self._extrapolation = extrapolation

    def validate(self) raises:
        """Explicitly revalidate the stored table and coefficient buffers."""
        _validate_table(self._knots, self._values)
        _validate_coefficient_buffers(
            self._a,
            self._b,
            self._c,
            self._d,
            len(self._knots) - 1,
        )

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
        if x == x0:
            return self._values[index]
        if x == x1:
            return self._values[index + 1]

        var dx = x - x0
        return (
            (self._d[index] * dx + self._c[index]) * dx + self._b[index]
        ) * dx + self._a[index]

    def _derivative_segment(self, index: Int, x: Float64) -> Float64:
        var dx = x - self._knots[index]
        return (3.0 * self._d[index] * dx + 2.0 * self._c[index]) * dx + self._b[index]

    def _second_derivative_segment(self, index: Int, x: Float64) -> Float64:
        var dx = x - self._knots[index]
        return 6.0 * self._d[index] * dx + 2.0 * self._c[index]

    def _left_slope(self) -> Float64:
        return self._b[0]

    def _right_slope(self) -> Float64:
        var final_interval = len(self._a) - 1
        return self._derivative_segment(
            final_interval, self._knots[len(self._knots) - 1]
        )

    def _evaluate_linear_ray(
        self,
        endpoint_x: Float64,
        endpoint_y: Float64,
        endpoint_slope: Float64,
        x: Float64,
    ) -> Float64:
        if endpoint_slope == 0.0:
            return endpoint_y
        return endpoint_y + endpoint_slope * (x - endpoint_x)

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
            return self._evaluate_linear_ray(
                self._knots[0], self._values[0], self._left_slope(), x
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
            return self._evaluate_linear_ray(
                self._knots[final_knot],
                self._values[final_knot],
                self._right_slope(),
                x,
            )

        return self._evaluate_segment(_locate_interval_in_domain(self._knots, x), x)

    def derivative(self, x: Float64) raises -> Float64:
        """Evaluate the first derivative under the configured policy.

        Outside the domain, `CLAMP` differentiates the clamped constant and
        returns zero. `LINEAR` returns the constant endpoint-ray slope, while
        `FILL` returns its payload and `ERROR` rejects the query. Non-finite
        queries always raise.
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
            return self._left_slope()

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
            return self._right_slope()

        return self._derivative_segment(_locate_interval_in_domain(self._knots, x), x)

    def second_derivative(self, x: Float64) raises -> Float64:
        """Evaluate the second derivative under the configured policy.

        Both `CLAMP` and `LINEAR` are constant/linear outside the domain and
        therefore return zero there. `FILL` returns its payload, `ERROR` rejects
        exterior queries, and non-finite queries always raise.
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

        if x == self._knots[0] or x == self._knots[final_knot]:
            return 0.0
        return self._second_derivative_segment(
            _locate_interval_in_domain(self._knots, x), x
        )

    def evaluate(self, queries: Span[Float64, _]) raises -> List[Float64]:
        """Allocate and return results for finite queries in order.

        Raises on the first offending query (a non-finite query under every
        policy, or an out-of-domain query under `ERROR`) and returns no partial
        results. Empty input returns an empty list.
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
