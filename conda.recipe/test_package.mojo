from nagare import (
    BoundaryCondition,
    CubicHermiteInterpolator,
    CubicSplineInterpolator,
    ExtrapolationPolicy,
    LinearInterpolator,
    MakimaInterpolator,
    PchipInterpolator,
    StepInterpolator,
    StepMode,
    locate_interval,
)
from std.collections import List
from std.testing import assert_equal


def main() raises:
    assert_equal(locate_interval([0.0, 1.0, 3.0], 2.0), 1)
    var interpolator = LinearInterpolator(
        [0.0, 1.0, 3.0],
        [2.0, 4.0, 8.0],
        extrapolation=ExtrapolationPolicy.CLAMP,
    )
    assert_equal(interpolator.evaluate(2.0), 6.0)
    assert_equal(interpolator.evaluate(4.0), 8.0)
    var queries: List[Float64] = [-1.0, 0.0, 2.0, 4.0]
    var results = List[Float64](length=len(queries), fill=0.0)
    interpolator.evaluate_sorted_into(queries, results)
    assert_equal(results, [2.0, 2.0, 6.0, 8.0])

    var step = StepInterpolator([0.0, 1.0], [10.0, 20.0], mode=StepMode.NEXT)
    assert_equal(step.evaluate(0.5), 20.0)

    var hermite = CubicHermiteInterpolator([0.0, 1.0], [0.0, 1.0], [1.0, 1.0])
    assert_equal(hermite.evaluate(0.5), 0.5)

    var cubic = CubicSplineInterpolator(
        [0.0, 1.0, 2.0],
        [0.0, 1.0, 4.0],
        boundary=BoundaryCondition.NATURAL,
    )
    assert_equal(cubic.evaluate(1.0), 1.0)

    var pchip = PchipInterpolator([0.0, 1.0, 2.0], [0.0, 1.0, 4.0])
    assert_equal(pchip.evaluate(1.0), 1.0)

    var makima = MakimaInterpolator(
        [0.0, 1.0, 2.0, 3.0, 4.0],
        [0.0, 1.0, 4.0, 9.0, 16.0],
    )
    assert_equal(makima.evaluate(2.0), 4.0)
