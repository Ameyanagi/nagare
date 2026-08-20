"""Explicit out-of-domain behavior for interpolants."""

from std.collections import Optional


struct ExtrapolationPolicy(Copyable, Equatable, ImplicitlyCopyable):
    """Validated policy for queries outside an interpolant's knot domain.

    `ERROR` rejects an out-of-domain query, `CLAMP` returns the nearest endpoint
    value, and `LINEAR` extends the nearest endpoint segment.
    """

    comptime ERROR = ExtrapolationPolicy()
    comptime CLAMP = ExtrapolationPolicy(False)
    comptime LINEAR = ExtrapolationPolicy(True)

    # Optional[Bool] has exactly three states, so every representable value is a
    # valid policy: None = error, False = clamp, True = linear.
    var _extends_linearly: Optional[Bool]

    def __init__(out self, extends_linearly: Optional[Bool] = None):
        """Construct one of the three representable policies.

        Public callers should normally use `ERROR`, `CLAMP`, or `LINEAR`.
        `None`, `False`, and `True` are their respective typed representations.
        """
        self._extends_linearly = extends_linearly.copy()

    def __eq__(self, other: Self) -> Bool:
        if self._extends_linearly:
            if other._extends_linearly:
                return self._extends_linearly.value() == other._extends_linearly.value()
            return False
        return not other._extends_linearly
