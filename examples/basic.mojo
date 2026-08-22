from nagare import ExtrapolationPolicy, LinearInterpolator
from std.collections import List


def main() raises:
    var temperature = LinearInterpolator(
        [0.0, 10.0, 30.0],
        [12.0, 18.0, 27.0],
        extrapolation=ExtrapolationPolicy.CLAMP,
    )

    print("temperature at t=5:", temperature.evaluate(5.0))
    print("clamped temperature at t=40:", temperature.evaluate(40.0))

    var query_times: List[Float64] = [-5.0, 0.0, 5.0, 20.0, 40.0]
    var temperatures = List[Float64](length=len(query_times), fill=0.0)
    temperature.evaluate_sorted_into(query_times, temperatures)
    print("sorted batch:", temperatures)
