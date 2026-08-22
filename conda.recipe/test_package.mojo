from nagare import ExtrapolationPolicy, LinearInterpolator, locate_interval
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
