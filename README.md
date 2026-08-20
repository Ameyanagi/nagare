# Nagare

> **Experimental — API not yet released.**

Interpolation and spline algorithms for Mojo.

## Scope

Nagare is a focused numerical interpolation toolkit with explicit extrapolation and continuity contracts.

The first implementation milestone is intentionally narrow: implement validated knot preparation, binary-search interval location, linear interpolation, and cubic splines with independently generated reference fixtures.
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

The first precise vertical slice is available: validated interval location and
an owning one-dimensional linear interpolator with explicit error, clamp, and
linear extrapolation policies.

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
