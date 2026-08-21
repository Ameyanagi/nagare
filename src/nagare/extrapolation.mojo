"""Explicit out-of-domain behavior for interpolants."""

from std.io import Writable, Writer
from std.memory import bitcast


comptime _NAN = bitcast[DType.float64](UInt64(0x7FF8_0000_0000_0000))


def _domain_error_message(
    query: Float64, domain_start: Float64, domain_end: Float64
) -> String:
    """Format the shared teaching diagnostic for rejected extrapolation."""
    return String(
        "query ",
        query,
        " is outside the knot domain [",
        domain_start,
        ", ",
        domain_end,
        "]; pass extrapolation= to allow this",
    )


def _validate_integration_bounds(a: Float64, b: Float64) raises:
    """Reject non-finite integration bounds with both values in context."""
    if a != a or a - a != 0.0 or b != b or b - b != 0.0:
        raise Error(String("integration bounds must be finite: a = ", a, ", b = ", b))


struct ExtrapolationPolicy(Copyable, Equatable, ImplicitlyCopyable, Writable):
    """Nominal policy for queries outside an interpolant's knot domain.

    `ERROR` rejects an out-of-domain query, `CLAMP` returns the nearest endpoint
    value, and `LINEAR` follows the interpolant's endpoint tangent as a linear
    ray. `FILL` returns its stored `Float64` payload, which defaults to NaN. For
    a piecewise-linear interpolant, the linear ray extends the endpoint segment.

    The constants and `fill()` are the public construction surface. The integer
    discriminant, payload, and underscore-prefixed initializer arguments are
    private by convention. Equality compares both fields by bit pattern, so the
    default NaN payload compares equal to itself and custom fills remain exact.
    """

    var _value: Int
    var _fill: Float64

    comptime ERROR = ExtrapolationPolicy(_value=0, _fill=_NAN)
    comptime CLAMP = ExtrapolationPolicy(_value=1, _fill=_NAN)
    comptime LINEAR = ExtrapolationPolicy(_value=2, _fill=_NAN)
    comptime FILL = ExtrapolationPolicy(_value=3, _fill=_NAN)

    def __init__(out self, *, _value: Int, _fill: Float64):
        self._value = _value
        self._fill = _fill

    @staticmethod
    def fill(value: Float64) -> ExtrapolationPolicy:
        """Return a fill policy carrying `value`.

        Any NaN payload is canonicalized to the default NaN, so every NaN
        fill compares equal to `FILL` under the bitwise payload equality.
        """
        if value != value:
            return ExtrapolationPolicy(_value=3, _fill=_NAN)
        return ExtrapolationPolicy(_value=3, _fill=value)

    def is_fill(self) -> Bool:
        """Return True for the FILL variant regardless of payload."""
        return self._value == 3

    def fill_value(self) -> Float64:
        """Return the stored payload, which is ignored by non-fill variants."""
        return self._fill

    def __eq__(self, other: Self) -> Bool:
        return self._value == other._value and bitcast[DType.uint64](
            self._fill
        ) == bitcast[DType.uint64](other._fill)

    def __str__(self) -> String:
        var result = String()
        self.write_to(result)
        return result^

    def write_to[W: Writer](self, mut writer: W):
        """Write the stable public spelling, including a fill payload."""
        if self == Self.ERROR:
            writer.write("ERROR")
        elif self == Self.CLAMP:
            writer.write("CLAMP")
        elif self == Self.LINEAR:
            writer.write("LINEAR")
        else:
            writer.write("FILL(value=", self._fill, ")")


def _integrate_exterior_tail(
    policy: ExtrapolationPolicy,
    start: Float64,
    end: Float64,
    endpoint_x: Float64,
    endpoint_y: Float64,
    endpoint_slope: Float64,
    domain_start: Float64,
    domain_end: Float64,
) raises -> Float64:
    """Integrate one known non-empty exterior tail under `policy`."""
    if policy == ExtrapolationPolicy.ERROR:
        var rejected = end
        if start < domain_start:
            rejected = start
        raise Error(_domain_error_message(rejected, domain_start, domain_end))

    var width = end - start
    if policy == ExtrapolationPolicy.CLAMP:
        return endpoint_y * width
    if policy.is_fill():
        return policy.fill_value() * width

    var left_offset = start - endpoint_x
    var right_offset = end - endpoint_x
    var midpoint_offset = 0.5 * left_offset + 0.5 * right_offset
    return width * (endpoint_y + endpoint_slope * midpoint_offset)
