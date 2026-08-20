"""Explicit out-of-domain behavior for interpolants."""


struct ExtrapolationPolicy(Copyable, Equatable, ImplicitlyCopyable):
    """Nominal policy for queries outside an interpolant's knot domain.

    `ERROR` rejects an out-of-domain query, `CLAMP` returns the nearest endpoint
    value, and `LINEAR` follows the interpolant's endpoint tangent as a linear
    ray. For a piecewise-linear interpolant, that ray extends the endpoint
    segment.

    The constants are the public construction surface. The integer
    discriminant and its underscore-prefixed initializer argument are private
    by convention.
    """

    comptime ERROR = ExtrapolationPolicy(_value=0)
    comptime CLAMP = ExtrapolationPolicy(_value=1)
    comptime LINEAR = ExtrapolationPolicy(_value=2)

    var _value: Int

    def __init__(out self, *, _value: Int):
        self._value = _value

    def __eq__(self, other: Self) -> Bool:
        return self._value == other._value
