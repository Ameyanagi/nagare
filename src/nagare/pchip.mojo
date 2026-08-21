"""Shape-preserving piecewise-cubic Hermite interpolation."""

from std.collections import List
from std.io import Writable, Writer

from .extrapolation import ExtrapolationPolicy
from .hermite import CubicHermiteInterpolator
from .linear import _validate_table


def _same_sign(left: Float64, right: Float64) -> Bool:
    """Return whether two values have the same three-way sign."""
    if left == 0.0:
        return right == 0.0
    if right == 0.0:
        return False
    return (left < 0.0) == (right < 0.0)


def _edge_slope(
    first_width: Float64,
    second_width: Float64,
    first_secant: Float64,
    second_secant: Float64,
) -> Float64:
    """Return scipy's shape-preserving one-sided endpoint estimate."""
    var slope = (
        (2.0 * first_width + second_width) * first_secant - first_width * second_secant
    ) / (first_width + second_width)
    if not _same_sign(slope, first_secant):
        return 0.0
    if not _same_sign(first_secant, second_secant) and abs(slope) > 3.0 * abs(
        first_secant
    ):
        return 3.0 * first_secant
    return slope


def _pchip_slopes(knots: List[Float64], values: List[Float64]) -> List[Float64]:
    """Derive Fritsch-Carlson knot slopes from a validated table."""
    var interval_count = len(knots) - 1
    var widths = List[Float64](capacity=interval_count)
    var secants = List[Float64](capacity=interval_count)
    for index in range(interval_count):
        var width = knots[index + 1] - knots[index]
        widths.append(width)
        secants.append((values[index + 1] - values[index]) / width)

    if interval_count == 1:
        return [secants[0], secants[0]]

    var slopes = List[Float64](capacity=len(knots))
    slopes.append(_edge_slope(widths[0], widths[1], secants[0], secants[1]))
    for knot_index in range(1, len(knots) - 1):
        var left_secant = secants[knot_index - 1]
        var right_secant = secants[knot_index]
        if (
            left_secant == 0.0
            or right_secant == 0.0
            or not _same_sign(left_secant, right_secant)
        ):
            slopes.append(0.0)
            continue

        var left_width = widths[knot_index - 1]
        var right_width = widths[knot_index]
        var first_weight = 2.0 * right_width + left_width
        var second_weight = right_width + 2.0 * left_width
        slopes.append(
            (first_weight + second_weight)
            / (first_weight / left_secant + second_weight / right_secant)
        )
    slopes.append(
        _edge_slope(
            widths[interval_count - 1],
            widths[interval_count - 2],
            secants[interval_count - 1],
            secants[interval_count - 2],
        )
    )
    return slopes^


struct PchipInterpolator(Copyable, Equatable, Writable):
    """Shape-preserving C1 piecewise-cubic Hermite interpolant (Fritsch–Carlson).

    Construction validates the table and derives monotonicity-preserving knot
    slopes in O(n). Each query takes O(log n). The interpolant never overshoots
    locally monotone data, reproduces every knot exactly, and has a continuous
    first derivative.

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
        """Validate, derive Fritsch-Carlson slopes, and take ownership."""
        _validate_table(knots, values)
        var slopes = _pchip_slopes(knots, values)
        self._engine = CubicHermiteInterpolator(
            knots^,
            values^,
            slopes^,
            extrapolation=extrapolation,
        )

    def validate(self) raises:
        """Explicitly revalidate the stored Hermite table."""
        self._engine.validate()

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

    def derivative(self, x: Float64) raises -> Float64:
        """Evaluate the first derivative under the configured policy.

        Outside the domain, `CLAMP` returns zero, `LINEAR` returns the derived
        endpoint slope, `FILL` returns its payload, and `ERROR` rejects the
        query. Non-finite queries always raise.
        """
        return self._engine.derivative(x)

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
            "PchipInterpolator(",
            self.knot_count(),
            " knots on [",
            self.domain_start(),
            ", ",
            self.domain_end(),
            "], extrapolation=",
            self._engine._extrapolation,
            ")",
        )
