"""Deterministic, dependency-free fixtures for linear benchmarks."""

from std.collections import List


struct LinearBenchmarkFixture(Copyable):
    """One owned knot/value table used by the benchmark matrix."""

    var knots: List[Float64]
    var values: List[Float64]

    def __init__(out self, var knots: List[Float64], var values: List[Float64]):
        self.knots = knots^
        self.values = values^


def make_linear_fixture(count: Int, irregular: Bool) raises -> LinearBenchmarkFixture:
    """Build a deterministic moderate-magnitude table with `count` knots."""
    if count < 2:
        raise Error("benchmark fixtures require at least two knots")

    var knots = List[Float64]()
    var values = List[Float64]()
    var current = -2.0
    for index in range(count):
        if index > 0:
            var step = 0.75
            if irregular:
                step += Float64((index * 17) % 5) * 0.125
            current += step
        knots.append(current)
        var perturbation = Float64((index * 7) % 3 - 1) * 0.25
        values.append(1.5 * current - 2.0 + perturbation)
    return LinearBenchmarkFixture(knots^, values^)


def make_tiny_ordinate_fixture() -> LinearBenchmarkFixture:
    """Build the extreme extrapolation fixture whose finite result is near two."""
    return LinearBenchmarkFixture(
        [0.0, 1e-308],
        [-1e-308, 1e-308],
    )


def make_full_range_fixture() -> LinearBenchmarkFixture:
    """Build the documented ill-conditioned full-range in-domain fixture."""
    return LinearBenchmarkFixture(
        [-Float64.MAX_FINITE, Float64.MAX_FINITE],
        [-Float64.MAX_FINITE, Float64.MAX_FINITE],
    )
