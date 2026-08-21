"""Validated one-dimensional cubic spline interpolation."""

from std.collections import List
from std.io import Writable, Writer
from std.memory import bitcast

from .extrapolation import (
    ExtrapolationPolicy,
    _NAN,
    _domain_error_message,
    _integrate_exterior_tail,
    _validate_integration_bounds,
)
from .linear import _validate_table
from .search import _is_finite, _locate_interval_in_domain


struct BoundaryCondition(Copyable, Equatable, ImplicitlyCopyable, Writable):
    """Nominal endpoint constraints for a cubic spline.

    `NOT_A_KNOT` is the scipy-compatible default. `NATURAL` fixes both endpoint
    second derivatives to zero. `PERIODIC` matches endpoint first and second
    derivatives and requires equal endpoint values. `clamped()` pins the two
    endpoint first derivatives to finite payload slopes.
    """

    var _value: Int
    var _start_slope: Float64
    var _end_slope: Float64

    comptime NOT_A_KNOT = BoundaryCondition(
        _value=0, _start_slope=_NAN, _end_slope=_NAN
    )
    comptime NATURAL = BoundaryCondition(_value=1, _start_slope=_NAN, _end_slope=_NAN)
    comptime PERIODIC = BoundaryCondition(_value=2, _start_slope=_NAN, _end_slope=_NAN)

    def __init__(
        out self,
        *,
        _value: Int,
        _start_slope: Float64,
        _end_slope: Float64,
    ):
        self._value = _value
        self._start_slope = _start_slope
        self._end_slope = _end_slope

    @staticmethod
    def clamped(start_slope: Float64, end_slope: Float64) raises -> BoundaryCondition:
        """Return a clamped condition with finite endpoint derivatives."""
        if not _is_finite(start_slope) or not _is_finite(end_slope):
            raise Error(
                String(
                    "clamped boundary slopes must be finite: start_slope = ",
                    start_slope,
                    ", end_slope = ",
                    end_slope,
                )
            )
        return BoundaryCondition(
            _value=3,
            _start_slope=start_slope,
            _end_slope=end_slope,
        )

    def is_clamped(self) -> Bool:
        """Return whether this condition carries endpoint slopes."""
        return self._value == 3

    def start_slope(self) -> Float64:
        """Return the stored start slope, ignored by non-clamped variants."""
        return self._start_slope

    def end_slope(self) -> Float64:
        """Return the stored end slope, ignored by non-clamped variants."""
        return self._end_slope

    def __eq__(self, other: Self) -> Bool:
        return (
            self._value == other._value
            and bitcast[DType.uint64](self._start_slope)
            == bitcast[DType.uint64](other._start_slope)
            and bitcast[DType.uint64](self._end_slope)
            == bitcast[DType.uint64](other._end_slope)
        )

    def __str__(self) -> String:
        var result = String()
        self.write_to(result)
        return result^

    def write_to[W: Writer](self, mut writer: W):
        """Write the stable public spelling, including clamped payloads."""
        if self == Self.NOT_A_KNOT:
            writer.write("NOT_A_KNOT")
        elif self == Self.NATURAL:
            writer.write("NATURAL")
        elif self == Self.PERIODIC:
            writer.write("PERIODIC")
        else:
            writer.write(
                "CLAMPED(start=",
                self._start_slope,
                ", end=",
                self._end_slope,
                ")",
            )


def _validate_boundary_table(values: List[Float64], boundary: BoundaryCondition) raises:
    """Validate boundary-specific construction invariants."""
    if boundary.is_clamped():
        if not _is_finite(boundary.start_slope()) or not _is_finite(
            boundary.end_slope()
        ):
            raise Error(
                String(
                    "clamped boundary slopes must be finite: start_slope = ",
                    boundary.start_slope(),
                    ", end_slope = ",
                    boundary.end_slope(),
                )
            )
        return
    if boundary == BoundaryCondition.PERIODIC:
        if len(values) < 3:
            raise Error(
                String(
                    "periodic cubic spline requires at least three knots: received ",
                    len(values),
                )
            )
        var final_index = len(values) - 1
        if values[0] != values[final_index]:
            raise Error(
                String(
                    "periodic boundary requires equal endpoint values: values[0] = ",
                    values[0],
                    ", values[",
                    final_index,
                    "] = ",
                    values[final_index],
                )
            )
        return
    if (
        boundary != BoundaryCondition.NOT_A_KNOT
        and boundary != BoundaryCondition.NATURAL
    ):
        raise Error("boundary condition is invalid")


def _solve_tridiagonal(
    var lower: List[Float64],
    var diagonal: List[Float64],
    var upper: List[Float64],
    var rhs: List[Float64],
) -> List[Float64]:
    """Solve a trusted nonsingular tridiagonal system with Thomas elimination."""
    var count = len(diagonal)
    for index in range(1, count):
        var multiplier = lower[index] / diagonal[index - 1]
        diagonal[index] -= multiplier * upper[index - 1]
        rhs[index] -= multiplier * rhs[index - 1]

    var solution = List[Float64](length=count, fill=0.0)
    solution[count - 1] = rhs[count - 1] / diagonal[count - 1]
    var index = count - 1
    while index > 0:
        index -= 1
        solution[index] = (rhs[index] - upper[index] * solution[index + 1]) / diagonal[
            index
        ]
    return solution^


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


def _solve_clamped_moments(
    knots: List[Float64],
    values: List[Float64],
    start_slope: Float64,
    end_slope: Float64,
) -> List[Float64]:
    """Return knot second derivatives for pinned endpoint first derivatives."""
    var knot_count = len(knots)
    var lower = List[Float64](length=knot_count, fill=0.0)
    var diagonal = List[Float64](length=knot_count, fill=0.0)
    var upper = List[Float64](length=knot_count, fill=0.0)
    var rhs = List[Float64](length=knot_count, fill=0.0)

    var first_width = knots[1] - knots[0]
    var first_secant = (values[1] - values[0]) / first_width
    diagonal[0] = 2.0 * first_width
    upper[0] = first_width
    rhs[0] = 6.0 * (first_secant - start_slope)

    for knot_index in range(1, knot_count - 1):
        var left_width = knots[knot_index] - knots[knot_index - 1]
        var right_width = knots[knot_index + 1] - knots[knot_index]
        lower[knot_index] = left_width
        diagonal[knot_index] = 2.0 * (left_width + right_width)
        upper[knot_index] = right_width
        rhs[knot_index] = 6.0 * (
            (values[knot_index + 1] - values[knot_index]) / right_width
            - (values[knot_index] - values[knot_index - 1]) / left_width
        )

    var final_index = knot_count - 1
    var final_width = knots[final_index] - knots[final_index - 1]
    var final_secant = (values[final_index] - values[final_index - 1]) / final_width
    lower[final_index] = final_width
    diagonal[final_index] = 2.0 * final_width
    rhs[final_index] = 6.0 * (end_slope - final_secant)
    return _solve_tridiagonal(lower^, diagonal^, upper^, rhs^)


def _solve_not_a_knot_moments(
    knots: List[Float64], values: List[Float64]
) -> List[Float64]:
    """Return scipy-compatible not-a-knot moments.

    Two knots degrade to a line. Three knots use the one parabola through all
    samples. For larger tables, the endpoint third-derivative constraints are
    substituted into the first and last interior rows, restoring a standard
    tridiagonal system before Thomas elimination.
    """
    var knot_count = len(knots)
    if knot_count == 2:
        return [0.0, 0.0]
    if knot_count == 3:
        var left_width = knots[1] - knots[0]
        var right_width = knots[2] - knots[1]
        var left_secant = (values[1] - values[0]) / left_width
        var right_secant = (values[2] - values[1]) / right_width
        var moment = 2.0 * (right_secant - left_secant) / (left_width + right_width)
        return [moment, moment, moment]

    var interior_count = knot_count - 2
    var lower = List[Float64](length=interior_count, fill=0.0)
    var diagonal = List[Float64](length=interior_count, fill=0.0)
    var upper = List[Float64](length=interior_count, fill=0.0)
    var rhs = List[Float64](length=interior_count, fill=0.0)
    for unknown_index in range(interior_count):
        var knot_index = unknown_index + 1
        var left_width = knots[knot_index] - knots[knot_index - 1]
        var right_width = knots[knot_index + 1] - knots[knot_index]
        if unknown_index > 0:
            lower[unknown_index] = left_width
        diagonal[unknown_index] = 2.0 * (left_width + right_width)
        if unknown_index < interior_count - 1:
            upper[unknown_index] = right_width
        rhs[unknown_index] = 6.0 * (
            (values[knot_index + 1] - values[knot_index]) / right_width
            - (values[knot_index] - values[knot_index - 1]) / left_width
        )

    var first_width = knots[1] - knots[0]
    var second_width = knots[2] - knots[1]
    diagonal[0] = (first_width + second_width) * (first_width / second_width + 2.0)
    upper[0] = second_width - first_width * first_width / second_width

    var final_unknown = interior_count - 1
    var penultimate_width = knots[knot_count - 2] - knots[knot_count - 3]
    var final_width = knots[knot_count - 1] - knots[knot_count - 2]
    lower[final_unknown] = (
        penultimate_width - final_width * final_width / penultimate_width
    )
    diagonal[final_unknown] = (penultimate_width + final_width) * (
        2.0 + final_width / penultimate_width
    )

    var interior = _solve_tridiagonal(lower^, diagonal^, upper^, rhs^)
    var moments = List[Float64](length=knot_count, fill=0.0)
    for index in range(interior_count):
        moments[index + 1] = interior[index]
    moments[0] = (
        (first_width + second_width) * moments[1] - first_width * moments[2]
    ) / second_width
    moments[knot_count - 1] = (
        (final_width + penultimate_width) * moments[knot_count - 2]
        - final_width * moments[knot_count - 3]
    ) / penultimate_width
    return moments^


def _solve_periodic_moments(
    knots: List[Float64], values: List[Float64]
) -> List[Float64]:
    """Return periodic moments from a cyclic tridiagonal system."""
    var knot_count = len(knots)
    var first_width = knots[1] - knots[0]
    var final_width = knots[knot_count - 1] - knots[knot_count - 2]
    var first_secant = (values[1] - values[0]) / first_width
    var final_secant = (values[knot_count - 1] - values[knot_count - 2]) / final_width

    # With three knots the cyclic corners overlap the ordinary off-diagonals,
    # so solve the resulting symmetric 2-by-2 system explicitly.
    if knot_count == 3:
        var total_width = first_width + final_width
        var rhs = 6.0 * (first_secant - final_secant)
        var first_moment = rhs / total_width
        var middle_moment = -first_moment
        return [first_moment, middle_moment, first_moment]

    # M[n - 1] equals M[0], leaving n - 1 unique unknowns. The derivative
    # matching equation is row zero; the remaining rows are interior C2
    # equations. Sherman-Morrison removes the two cyclic corner entries.
    var unique_count = knot_count - 1
    var lower = List[Float64](length=unique_count, fill=0.0)
    var diagonal = List[Float64](length=unique_count, fill=0.0)
    var upper = List[Float64](length=unique_count, fill=0.0)
    var rhs = List[Float64](length=unique_count, fill=0.0)
    diagonal[0] = 2.0 * (first_width + final_width)
    upper[0] = first_width
    rhs[0] = 6.0 * (first_secant - final_secant)

    for knot_index in range(1, unique_count):
        var left_width = knots[knot_index] - knots[knot_index - 1]
        var right_width = knots[knot_index + 1] - knots[knot_index]
        lower[knot_index] = left_width
        diagonal[knot_index] = 2.0 * (left_width + right_width)
        if knot_index < unique_count - 1:
            upper[knot_index] = right_width
        rhs[knot_index] = 6.0 * (
            (values[knot_index + 1] - values[knot_index]) / right_width
            - (values[knot_index] - values[knot_index - 1]) / left_width
        )

    var alpha = final_width
    var beta = final_width
    var gamma = -diagonal[0]
    diagonal[0] -= gamma
    diagonal[unique_count - 1] -= alpha * beta / gamma

    var lower_x = lower.copy()
    var diagonal_x = diagonal.copy()
    var upper_x = upper.copy()
    var rhs_x = rhs.copy()
    var solution = _solve_tridiagonal(lower_x^, diagonal_x^, upper_x^, rhs_x^)

    var correction_rhs = List[Float64](length=unique_count, fill=0.0)
    correction_rhs[0] = gamma
    correction_rhs[unique_count - 1] = alpha
    var correction = _solve_tridiagonal(lower^, diagonal^, upper^, correction_rhs^)
    var factor = (solution[0] + beta * solution[unique_count - 1] / gamma) / (
        1.0 + correction[0] + beta * correction[unique_count - 1] / gamma
    )
    for index in range(unique_count):
        solution[index] -= factor * correction[index]

    var moments = List[Float64](length=knot_count, fill=0.0)
    for index in range(unique_count):
        moments[index] = solution[index]
    moments[knot_count - 1] = moments[0]
    return moments^


def _solve_moments(
    knots: List[Float64],
    values: List[Float64],
    boundary: BoundaryCondition,
) -> List[Float64]:
    """Dispatch to the validated boundary-condition solver."""
    if boundary == BoundaryCondition.NATURAL:
        return _solve_natural_moments(knots, values)
    if boundary == BoundaryCondition.NOT_A_KNOT:
        return _solve_not_a_knot_moments(knots, values)
    if boundary == BoundaryCondition.PERIODIC:
        return _solve_periodic_moments(knots, values)
    return _solve_clamped_moments(
        knots,
        values,
        boundary.start_slope(),
        boundary.end_slope(),
    )


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


struct CubicSplineInterpolator(Copyable, Equatable, Writable):
    """An owning C2 cubic spline over finite `Float64` data.

    Each interval stores contiguous `a`, `b`, `c`, and `d` buffers for
    `S_i(dx) = a_i + b_i dx + c_i dx^2 + d_i dx^3`, where
    `dx = x - knots[i]`. Construction is O(n), and each interval query is
    O(log n). Boundary choices are scipy-compatible `NOT_A_KNOT` (the default),
    `NATURAL`, `PERIODIC`, and slope-carrying `clamped(start, end)`.

    Not-a-knot makes the third derivative continuous across the first and last
    interior knots. As in scipy, two knots degrade to a line and three knots
    fit the single parabola through all samples. Natural fixes endpoint second
    derivatives to zero, clamped pins endpoint first derivatives, and periodic
    matches endpoint first and second derivatives while requiring equal values.

    Construction validates the table and coefficient representation, which
    read-only methods trust thereafter. Underscore-prefixed fields are private
    by convention; mutating them directly is outside the contract. Call
    `validate()` for an explicit checkpoint after unusual direct access.

    Outside the domain, `LINEAR` follows the endpoint tangent ray rather than
    extending a boundary cubic. The ray is C1-continuous at each join and is
    C2-continuous there only under `NATURAL`, whose endpoint second derivatives
    are zero. `FILL` returns its payload for values and both derivative orders
    outside the domain.
    """

    var _knots: List[Float64]
    var _values: List[Float64]
    var _a: List[Float64]
    var _b: List[Float64]
    var _c: List[Float64]
    var _d: List[Float64]
    var _extrapolation: ExtrapolationPolicy
    var _boundary: BoundaryCondition

    def __init__(
        out self,
        var knots: List[Float64],
        var values: List[Float64],
        extrapolation: ExtrapolationPolicy = ExtrapolationPolicy.ERROR,
        boundary: BoundaryCondition = BoundaryCondition.NOT_A_KNOT,
    ) raises:
        """Validate, solve, and take ownership of a knot/value table."""
        _validate_table(knots, values)
        _validate_boundary_table(values, boundary)

        var moments = _solve_moments(knots, values, boundary)
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
        self._boundary = boundary

    def validate(self) raises:
        """Explicitly revalidate the stored table and coefficient buffers."""
        _validate_table(self._knots, self._values)
        _validate_boundary_table(self._values, self._boundary)
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

    def boundary(self) -> BoundaryCondition:
        """Return the configured boundary condition."""
        return self._boundary

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

    def _segment_antiderivative(self, index: Int, x: Float64) -> Float64:
        """Evaluate one shifted cubic's antiderivative from its left knot."""
        var dx = x - self._knots[index]
        return dx * (
            self._a[index]
            + dx
            * (
                0.5 * self._b[index]
                + dx * (self._c[index] / 3.0 + dx * 0.25 * self._d[index])
            )
        )

    def _integrate_segment(self, index: Int, start: Float64, end: Float64) -> Float64:
        return self._segment_antiderivative(index, end) - self._segment_antiderivative(
            index, start
        )

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

        return self._second_derivative_segment(
            _locate_interval_in_domain(self._knots, x), x
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
                    self._left_slope(),
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
                    self._right_slope(),
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
        """Compare validated tables, boundary condition, and policy exactly.

        Validation excludes NaN from knot and value tables, so elementwise
        `Float64` equality is sound. Coefficients are deterministic derived
        storage and are therefore not a separate value identity component.
        """
        if (
            len(self._knots) != len(other._knots)
            or len(self._values) != len(other._values)
            or self._boundary != other._boundary
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
            "CubicSplineInterpolator(",
            len(self._knots),
            " knots on [",
            self._knots[0],
            ", ",
            self._knots[len(self._knots) - 1],
            "], boundary=",
            self._boundary,
            ", extrapolation=",
            self._extrapolation,
            ")",
        )
