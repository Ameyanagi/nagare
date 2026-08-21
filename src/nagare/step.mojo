"""Validated one-dimensional piecewise-constant interpolation."""

from std.collections import List
from std.io import Writable, Writer

from .extrapolation import ExtrapolationPolicy, _domain_error_message
from .linear import _validate_table
from .search import _is_finite, _locate_interval_in_domain


struct StepMode(Copyable, Equatable, ImplicitlyCopyable, Writable):
    """Nominal selection rule for values between step-interpolator knots."""

    var _value: Int

    comptime PREVIOUS = StepMode(_value=0)
    comptime NEXT = StepMode(_value=1)
    comptime NEAREST = StepMode(_value=2)

    def __init__(out self, *, _value: Int):
        self._value = _value

    def __eq__(self, other: Self) -> Bool:
        return self._value == other._value

    def __str__(self) -> String:
        var result = String()
        self.write_to(result)
        return result^

    def write_to[W: Writer](self, mut writer: W):
        """Write the stable public mode spelling."""
        if self == Self.PREVIOUS:
            writer.write("PREVIOUS")
        elif self == Self.NEXT:
            writer.write("NEXT")
        else:
            writer.write("NEAREST")


def _validate_step_extrapolation(policy: ExtrapolationPolicy) raises:
    if policy == ExtrapolationPolicy.LINEAR:
        raise Error(
            "LINEAR extrapolation is undefined for step interpolation; "
            "use ERROR, CLAMP, or FILL"
        )


struct StepInterpolator(Copyable, Equatable, Writable):
    """An owning piecewise-constant interpolant over finite `Float64` data.

    `PREVIOUS`, the scientific default and zero-order hold, returns `values[i]`
    for queries in `[knots[i], knots[i + 1])`. `NEXT` returns `values[i + 1]`
    for queries in `(knots[i], knots[i + 1]]`. Both return the exact tabulated
    value at every knot. `NEAREST` selects the closer knot, with an exact
    midpoint tie going left. The result is piecewise constant and C^-1 at knots
    where adjacent values differ.

    Construction validates the table in O(n), and each finite query takes
    O(log n). Underscore-prefixed fields are private by convention; mutating
    them directly is outside the contract. Call `validate()` for an explicit
    checkpoint after unusual direct access.

    Outside the domain, `ERROR` rejects the query, `CLAMP` returns the nearest
    endpoint value, and `FILL` returns its payload. `LINEAR` is rejected at
    construction because linear extrapolation is undefined for a step function.
    """

    var _knots: List[Float64]
    var _values: List[Float64]
    var _mode: StepMode
    var _extrapolation: ExtrapolationPolicy

    def __init__(
        out self,
        var knots: List[Float64],
        var values: List[Float64],
        mode: StepMode = StepMode.PREVIOUS,
        extrapolation: ExtrapolationPolicy = ExtrapolationPolicy.ERROR,
    ) raises:
        """Validate and take ownership of a step interpolation table."""
        _validate_table(knots, values)
        _validate_step_extrapolation(extrapolation)

        self._knots = knots^
        self._values = values^
        self._mode = mode
        self._extrapolation = extrapolation

    def validate(self) raises:
        """Explicitly revalidate the stored table and extrapolation policy."""
        _validate_table(self._knots, self._values)
        _validate_step_extrapolation(self._extrapolation)

    def knots(self) -> Span[Float64, origin_of(self._knots)]:
        """Read-only view of the validated knot sequence."""
        return Span(self._knots)

    def values(self) -> Span[Float64, origin_of(self._values)]:
        """Read-only view of the validated value sequence."""
        return Span(self._values)

    def knot_count(self) -> Int:
        """Return the knot count."""
        return len(self._knots)

    def domain_start(self) -> Float64:
        """Return the lower domain bound."""
        return self._knots[0]

    def domain_end(self) -> Float64:
        """Return the upper domain bound."""
        return self._knots[len(self._knots) - 1]

    def _evaluate_interval(self, index: Int, x: Float64) -> Float64:
        if x == self._knots[index]:
            return self._values[index]
        if x == self._knots[index + 1]:
            return self._values[index + 1]
        if self._mode == StepMode.PREVIOUS:
            return self._values[index]
        if self._mode == StepMode.NEXT:
            return self._values[index + 1]
        if x - self._knots[index] <= self._knots[index + 1] - x:
            return self._values[index]
        return self._values[index + 1]

    def evaluate(self, x: Float64) raises -> Float64:
        """Evaluate one finite query under the configured mode and policy."""
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
            return self._extrapolation.fill_value()

        var final_knot = len(self._knots) - 1
        if x > self._knots[final_knot]:
            if self._extrapolation == ExtrapolationPolicy.ERROR:
                raise Error(
                    _domain_error_message(x, self._knots[0], self._knots[final_knot])
                )
            if self._extrapolation == ExtrapolationPolicy.CLAMP:
                return self._values[final_knot]
            return self._extrapolation.fill_value()

        return self._evaluate_interval(_locate_interval_in_domain(self._knots, x), x)

    def __call__(self, x: Float64) raises -> Float64:
        """Call `evaluate`; `evaluate` is the primary documented name."""
        return self.evaluate(x)

    def evaluate(self, queries: Span[Float64, _]) raises -> List[Float64]:
        """Allocate and return results for finite queries in order."""
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
        """Compare validated tables, mode, and policy exactly.

        Validation excludes NaN from tables, so elementwise `Float64` equality
        is sound.
        """
        if (
            len(self._knots) != len(other._knots)
            or len(self._values) != len(other._values)
            or self._mode != other._mode
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
            "StepInterpolator(",
            len(self._knots),
            " knots on [",
            self._knots[0],
            ", ",
            self._knots[len(self._knots) - 1],
            "], mode=",
            self._mode,
            ", extrapolation=",
            self._extrapolation,
            ")",
        )
