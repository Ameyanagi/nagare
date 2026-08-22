"""Modified Akima (makima) piecewise-cubic Hermite interpolation."""

from std.collections import List
from std.io import Writable, Writer

from .extrapolation import ExtrapolationPolicy
from .hermite import CubicHermiteInterpolator
from .linear import _validate_table


def _validate_makima_table(knots: List[Float64], values: List[Float64]) raises:
    """Validate the table and makima's four-knot minimum."""
    if len(knots) < 4:
        raise Error(
            String(
                "makima interpolation requires at least four knots: ",
                "received ",
                len(knots),
            )
        )
    _validate_table(knots, values)


def _makima_slopes(knots: List[Float64], values: List[Float64]) -> List[Float64]:
    """Derive scipy-compatible modified-Akima knot slopes."""
    var knot_count = len(knots)
    var extended = List[Float64](length=knot_count + 3, fill=0.0)
    for index in range(knot_count - 1):
        extended[index + 2] = (values[index + 1] - values[index]) / (
            knots[index + 1] - knots[index]
        )

    extended[1] = 2.0 * extended[2] - extended[3]
    extended[0] = 2.0 * extended[1] - extended[2]
    extended[knot_count + 1] = 2.0 * extended[knot_count] - extended[knot_count - 1]
    extended[knot_count + 2] = 2.0 * extended[knot_count + 1] - extended[knot_count]

    var difference_weights = List[Float64](capacity=knot_count + 2)
    var sum_weights = List[Float64](capacity=knot_count + 2)
    for index in range(knot_count + 2):
        difference_weights.append(abs(extended[index + 1] - extended[index]))
        sum_weights.append(abs(extended[index + 1] + extended[index]))

    var right_weights = List[Float64](capacity=knot_count)
    var left_weights = List[Float64](capacity=knot_count)
    var combined_weights = List[Float64](capacity=knot_count)
    var maximum_combined_weight = 0.0
    for knot_index in range(knot_count):
        var right_weight = (
            difference_weights[knot_index + 2] + 0.5 * sum_weights[knot_index + 2]
        )
        var left_weight = difference_weights[knot_index] + 0.5 * sum_weights[knot_index]
        var combined_weight = right_weight + left_weight
        right_weights.append(right_weight)
        left_weights.append(left_weight)
        combined_weights.append(combined_weight)
        maximum_combined_weight = max(maximum_combined_weight, combined_weight)

    var slopes = List[Float64](capacity=knot_count)
    var relative_threshold = 1e-9 * maximum_combined_weight
    for knot_index in range(knot_count):
        if combined_weights[knot_index] <= relative_threshold:
            slopes.append(0.0)
            continue
        slopes.append(
            (
                right_weights[knot_index] * extended[knot_index + 1]
                + left_weights[knot_index] * extended[knot_index + 2]
            )
            / combined_weights[knot_index]
        )
    return slopes^


struct MakimaInterpolator(Copyable, Equatable, Writable):
    """Modified Akima (makima) C1 piecewise-cubic interpolant.

    Matches scipy `Akima1DInterpolator(method="makima")`, not the default
    `method="akima"`. Classic Akima is deliberately not offered.

    Construction requires at least four knots, validates the table, and derives
    local modified-Akima slopes in O(n). Each query takes O(log n). Changing one
    data point perturbs only nearby segments, while the modified weights reduce
    overshoot and define equal-difference cases.

    The shared `CubicHermiteInterpolator` owns the table and is the sole value,
    derivative, batch, and extrapolation engine. Its validated representation
    is trusted by read-only methods; call `validate()` for an explicit
    checkpoint after unusual direct access.
    """

    var _engine: CubicHermiteInterpolator

    def __init__(
        out self,
        var knots: List[Float64],
        var values: List[Float64],
        extrapolation: ExtrapolationPolicy = ExtrapolationPolicy.ERROR,
    ) raises:
        """Validate, derive modified-Akima slopes, and take ownership."""
        _validate_makima_table(knots, values)
        var slopes = _makima_slopes(knots, values)
        self._engine = CubicHermiteInterpolator(
            knots^,
            values^,
            slopes^,
            extrapolation=extrapolation,
        )

    def validate(self) raises:
        """Explicitly revalidate the stored Hermite table."""
        _validate_makima_table(self._engine._knots, self._engine._values)
        self._engine.validate()

    def knots(self) -> Span[Float64, origin_of(self._engine._knots)]:
        """Read-only view of the validated knot sequence."""
        return Span(self._engine._knots)

    def values(self) -> Span[Float64, origin_of(self._engine._values)]:
        """Read-only view of the validated value sequence."""
        return Span(self._engine._values)

    def knot_count(self) -> Int:
        """Return the knot count."""
        return self._engine.knot_count()

    def domain_start(self) -> Float64:
        """Return the lower domain bound."""
        return self._engine.domain_start()

    def domain_end(self) -> Float64:
        """Return the upper domain bound."""
        return self._engine.domain_end()

    def evaluate(self, x: Float64) raises -> Float64:
        """Evaluate one finite query under the configured policy.

        Exact knot queries return the stored ordinate. `LINEAR` extrapolation
        uses the derived endpoint tangent ray, while `FILL` returns its payload
        outside the domain.
        """
        return self._engine.evaluate(x)

    def __call__(self, x: Float64) raises -> Float64:
        """Call `evaluate`; `evaluate` is the primary documented name."""
        return self.evaluate(x)

    def derivative(self, x: Float64) raises -> Float64:
        """Evaluate the first derivative under the configured policy.

        Outside the domain, `CLAMP` returns zero, `LINEAR` returns the derived
        endpoint slope, `FILL` returns its payload, and `ERROR` rejects the
        query. Non-finite queries always raise.
        """
        return self._engine.derivative(x)

    def derivative(self, queries: Span[Float64, _]) raises -> List[Float64]:
        """Allocate and return first derivatives for finite queries in order.

        Raises on the first offending query and returns no partial results.
        Empty input returns an empty list.
        """
        return self._engine.derivative(queries)

    def derivative_into(
        self,
        queries: Span[Float64, _],
        results: Span[mut=True, Float64, _],
    ) raises:
        """Evaluate each first derivative into `results` without allocating.

        Raises if the buffer lengths differ or on the first offending query;
        results contents are unspecified after a raise.
        """
        self._engine.derivative_into(queries, results)

    def second_derivative(self, x: Float64) raises -> Float64:
        """Evaluate the policy-aware second derivative.

        At an interior knot, the Hermite engine selects the segment to its
        right; the final knot selects the final segment.
        """
        return self._engine.second_derivative(x)

    def integrate(self, a: Float64, b: Float64) raises -> Float64:
        """Definite integral over `[a, b]` (sign-flipped when `a > b`).

        Segments are integrated in closed form. Portions outside the domain
        follow the configured extrapolation policy; `ERROR` rejects them, and
        `FILL` contributes `fill_value * width` (NaN-propagating by default).
        """
        return self._engine.integrate(a, b)

    def evaluate(self, queries: Span[Float64, _]) raises -> List[Float64]:
        """Allocate and return results for finite queries in order.

        Raises on the first offending query and returns no partial results.
        Empty input returns an empty list.
        """
        return self._engine.evaluate(queries)

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
        self._engine.evaluate_into(queries, results)

    def __eq__(self, other: Self) -> Bool:
        """Compare the validated table, derived slopes, and policy exactly.

        Engine validation excludes NaN, so elementwise `Float64` equality is
        sound.
        """
        return self._engine == other._engine

    def __str__(self) -> String:
        var result = String()
        self.write_to(result)
        return result^

    def write_to[W: Writer](self, mut writer: W):
        writer.write(
            "MakimaInterpolator(",
            self.knot_count(),
            " knots on [",
            self.domain_start(),
            ", ",
            self.domain_end(),
            "], extrapolation=",
            self._engine._extrapolation,
            ")",
        )
