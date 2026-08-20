from nagare import CubicSplineInterpolator, ExtrapolationPolicy, PchipInterpolator
from std.collections import List


def main() raises:
    var sample_times: List[Float64] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    var readings: List[Float64] = [0.0, 0.1, 0.2, 5.0, 9.9, 10.0]
    var pchip = PchipInterpolator(
        sample_times.copy(),
        readings.copy(),
        extrapolation=ExtrapolationPolicy.FILL,
    )
    var cubic = CubicSplineInterpolator(sample_times^, readings^)

    var uniform_times = List[Float64](capacity=23)
    for step in range(23):
        uniform_times.append(-0.25 + 0.25 * Float64(step))
    var resampled = pchip.evaluate(uniform_times)

    print("Uniform 0.25 s sensor resampling; NaN marks unavailable edges.")
    print("Measured range: [0.0, 10.0]. Cubic overshoots; PCHIP stays in range:")
    for index in range(len(uniform_times)):
        var query = uniform_times[index]
        var pchip_value = resampled[index]
        if pchip_value != pchip_value:
            print("t =", query, ": PCHIP = NaN (outside sampled time range)")
            continue

        var cubic_value = cubic.evaluate(query)
        if cubic_value < 0.0 or cubic_value > 10.0:
            print(
                "t =",
                query,
                ": cubic =",
                cubic_value,
                ", PCHIP =",
                pchip_value,
            )
