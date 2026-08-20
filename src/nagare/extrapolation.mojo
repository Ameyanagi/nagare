"""Explicit out-of-domain behavior for interpolants."""

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


struct ExtrapolationPolicy(Copyable, Equatable, ImplicitlyCopyable):
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
