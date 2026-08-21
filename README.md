# Nagare

> **Experimental — API not yet released.**

Interpolation and spline algorithms for Mojo.

## Which interpolator?

| Choose | When it fits |
| --- | --- |
| `LinearInterpolator` | A robust baseline with no cubic overshoot. |
| `PchipInterpolator` | Measured monotone data where shape preservation and no overshoot matter. |
| `AkimaInterpolator` | A pleasing smooth makima curve with local response to data changes. |
| `CubicSplineInterpolator` | The smoothest C2 curve; not-a-knot defaults match scipy. |
| `StepInterpolator` with `PREVIOUS` | A zero-order hold for sampled controls or state. |
| `CubicHermiteInterpolator` | Slopes are already known, such as from an ODE solver or instrument. |

For shape-preserving sensor resampling, unavailable edges can stay explicit:

```mojo
from nagare import ExtrapolationPolicy, PchipInterpolator
from std.collections import List

var samples = PchipInterpolator(
    [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
    [0.0, 0.1, 0.2, 5.0, 9.9, 10.0],
    extrapolation=ExtrapolationPolicy.FILL,
)
var query_times: List[Float64] = [-0.25, 0.0, 0.25, 5.0, 5.25]
var resampled = samples.evaluate(query_times)  # NaN outside [0.0, 5.0]
```

## Scope

Nagare is a focused numerical interpolation toolkit with explicit extrapolation and continuity contracts.

The implemented scope is intentionally narrow: validated knot preparation,
binary-search interval location, step and linear interpolation, and a focused
family of one-dimensional cubic interpolants with independently derived and
scipy-cross-checked reference fixtures.
The project is independently installable and does not require any application
from the wider ecosystem.

## Development

Install [Pixi](https://pixi.sh/), then run:

```sh
pixi install --locked
pixi run check
pixi run example
pixi run bench-linear
```

The exact stable Mojo compiler and all development dependencies are captured in
`pixi.lock`. Runtime and library code is Mojo-first and pure Mojo wherever
practical. Build-time data generation may use another language when justified,
but generated outputs must be deterministic, checksum-pinned, licensed, and
documented.

## Package

The Mojo import is `nagare`. The eventual Conda distribution is
`mojo-nagare`. Source lives under `src/nagare/`, whose
`__init__.mojo` defines the package boundary.

The package exposes validated owning interpolators with scalar, batch,
derivative, and closed-form integration surfaces where mathematically defined.

```mojo
from nagare import ExtrapolationPolicy, LinearInterpolator, locate_interval

var line = LinearInterpolator(
    [0.0, 1.0, 3.0],
    [2.0, 4.0, 8.0],
    extrapolation=ExtrapolationPolicy.CLAMP,
)
print(line.evaluate(2.0))
```

Inputs must contain at least two finite, strictly increasing knots and an equal
number of finite values. The default extrapolation policy raises instead of
silently extending data beyond its domain. In-domain evaluation avoids
intermediate overflow for finite endpoint data. Explicit linear extrapolation
returns signed infinity when its represented result exceeds the finite
`Float64` range. Because Mojo 1.0 permits external mutation of
underscore-prefixed fields, those fields are private by convention and direct
mutation is outside the contract. Construction validates the stored lengths,
finite values, and knot ordering; read-only methods trust that invariant, and
`validate()` provides an explicit checkpoint after unusual direct access.
Intervals spanning nearly the full finite range produce finite,
endpoint-bounded in-domain results, but sub-ULP central offsets are not
promised exact affine recovery. See the
[v0.1 execution plan](docs/v0.1-plan.md) for the complete numerical contracts.

## Extrapolation and spline boundaries

Every interpolator defaults to `ExtrapolationPolicy.ERROR`. `CLAMP` holds the
nearest endpoint value, `LINEAR` follows the endpoint tangent ray, and `FILL`
returns a payload outside the domain (`NaN` by default, or a custom value from
`ExtrapolationPolicy.fill(value)`). Step interpolation rejects `LINEAR` because
a step function has no endpoint tangent. Definite integrals apply the same
policy to exterior tails.

`CubicSplineInterpolator` defaults to `BoundaryCondition.NOT_A_KNOT`, matching
scipy. `NATURAL` fixes both endpoint second derivatives to zero, `PERIODIC`
matches endpoint first and second derivatives and requires equal endpoint
values, and `BoundaryCondition.clamped(start_slope, end_slope)` pins the two
endpoint first derivatives.

## Repository map

- `src/nagare/`: library or application source
- `tests/`: TestSuite unit, reference-value, and invariant tests
- `examples/`: small compilable usage programs
- `benchmarks/`: versioned linear benchmark, deterministic fixtures, and
  reproducible methodology
- `docs/`: architecture, design, compatibility, roadmap, and release policy
- `conda.recipe/`: local Rattler build recipe

See [the architecture](docs/architecture.md), [design principles](docs/design.md),
and [roadmap](docs/roadmap.md) before proposing a new dependency or feature.

## License

Licensed under either Apache-2.0 or MIT, at your option.
