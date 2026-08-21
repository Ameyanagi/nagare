from nagare import ExtrapolationPolicy, LinearInterpolator


def main() raises:
    var temperature = LinearInterpolator(
        [0.0, 10.0, 30.0],
        [12.0, 18.0, 27.0],
        extrapolation=ExtrapolationPolicy.CLAMP,
    )

    print("temperature at t=5:", temperature.evaluate(5.0))
    print("clamped temperature at t=40:", temperature.evaluate(40.0))
