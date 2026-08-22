# Nagare

Interpolation and spline algorithms for Mojo.

> **Experimental — API not yet released.**

## Which interpolator?

| Choose | When it fits |
| --- | --- |
| `LinearInterpolator` | A robust baseline with no cubic overshoot. |
| `PchipInterpolator` | Measured monotone data where shape preservation and no overshoot matter. |
| `MakimaInterpolator` | A pleasing smooth makima curve with local response to data changes. |
| `CubicSplineInterpolator` | The smoothest C2 curve; not-a-knot defaults match scipy. |
| `StepInterpolator` with `PREVIOUS` | A zero-order hold for sampled controls or state. |
| `CubicHermiteInterpolator` | Slopes are already known, such as from an ODE solver or instrument. |

## Install

In your own Pixi project, add the hosted Nagare channel alongside the Mojo and
conda-forge channels in `pixi.toml`:

```toml
[workspace]
channels = [
  "https://ameyanagi.github.io/mojo-channel",
  "https://conda.modular.com/max",
  "conda-forge",
]
```

Then add Nagare:

```sh
pixi add mojo-nagare
```

Alternatively, work from a source checkout:

```sh
git clone https://github.com/Ameyanagi/nagare.git
cd nagare
pixi install --locked
```

From the repository root, run your own file with the source package on the
import path:

```sh
pixi run mojo run -I src my_file.mojo
```

## Quickstart

```mojo
from nagare import ExtrapolationPolicy, PchipInterpolator
from std.collections import List


def main() raises:
    var samples = PchipInterpolator(
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        [0.0, 0.1, 0.2, 5.0, 9.9, 10.0],
        extrapolation=ExtrapolationPolicy.FILL,
    )
    var query_times: List[Float64] = [-0.25, 0.0, 0.25, 5.0, 5.25]
    var resampled = samples.evaluate(query_times)  # NaN outside [0.0, 5.0]
    print(resampled)
```

Run it with `pixi run mojo run -I src my_file.mojo` from the repository root.

## Data reuse

Interpolator constructors take ownership of their `List[Float64]` arguments.
To reuse the same data for two interpolators, pass `.copy()` on earlier uses and
move the lists with `^` on the last use:

```mojo
from nagare import CubicSplineInterpolator, PchipInterpolator
from std.collections import List


def main() raises:
    var sample_times: List[Float64] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    var readings: List[Float64] = [0.0, 0.1, 0.2, 5.0, 9.9, 10.0]
    var pchip = PchipInterpolator(sample_times.copy(), readings.copy())
    var cubic = CubicSplineInterpolator(sample_times^, readings^)
    print("at 2.5:", pchip.evaluate(2.5), cubic.evaluate(2.5))
```

## Sorted resampling without repeated searches

When linear-interpolation query coordinates are already finite and
nondecreasing, use `evaluate_sorted` to binary-search the first interior query
once, then traverse only the remaining spanned intervals instead of searching
for every query. Duplicate coordinates and exact knots are supported, and
extrapolated prefixes and suffixes use the interpolator's configured policy.
For repeated workloads, `evaluate_sorted_into` reuses caller-owned storage:

```mojo
from nagare import ExtrapolationPolicy, LinearInterpolator
from std.collections import List


def main() raises:
    var samples = LinearInterpolator(
        [0.0, 1.0, 3.0],
        [10.0, 12.0, 18.0],
        extrapolation=ExtrapolationPolicy.CLAMP,
    )
    var queries: List[Float64] = [-1.0, 0.0, 0.5, 1.0, 2.0, 4.0]
    var output = List[Float64](length=len(queries), fill=0.0)
    samples.evaluate_sorted_into(queries, output)
    print(output)
```

The order-agnostic `evaluate` and `evaluate_into` APIs remain the right choice
for unsorted inputs.

## Sensor resampling example

[`examples/resample_sensor.mojo`](examples/resample_sensor.mojo) resamples
sensor readings on a uniform grid and shows how PCHIP preserves shape where a
cubic spline overshoots. Run it from the repository root:

```sh
pixi run mojo run -I src examples/resample_sensor.mojo
```

## Scope

Nagare is a focused numerical interpolation toolkit with explicit extrapolation
and continuity contracts. It provides validated knot preparation, binary-search
interval location, step and linear interpolation, and a focused family of
one-dimensional cubic interpolants with independently derived and
scipy-cross-checked reference fixtures. It is independently installable and
does not require an application from the wider ecosystem.

## Package and numerical contract

The Mojo import is `nagare`, the Conda distribution is `mojo-nagare`, and
source lives under `src/nagare/`. The package exposes owning interpolators with
scalar and batch evaluation, callable sugar, read-only knot/value views, scalar
and batch first derivatives, scalar second derivatives on cubic interpolants,
and closed-form integration surfaces where mathematically defined.

Inputs must contain at least two finite, strictly increasing knots and an equal
number of finite values. The default extrapolation policy raises instead of
silently extending data beyond its domain. In-domain evaluation avoids
intermediate overflow for finite endpoint data. Explicit linear extrapolation
returns signed infinity when its represented result exceeds the finite
`Float64` range.

Because Mojo 1.0 permits external mutation of underscore-prefixed fields, those
fields are private by convention and direct mutation is outside the contract.
Construction validates stored lengths, finite values, and knot ordering;
read-only methods trust that invariant, and `validate()` provides an explicit
checkpoint after unusual direct access. Intervals spanning nearly the full
finite range produce finite, endpoint-bounded in-domain results, but sub-ULP
central offsets are not promised exact affine recovery. See the
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

## Development

Install [Pixi](https://pixi.sh/), then run:

```sh
pixi install --locked
pixi run check
pixi run example
pixi run bench-linear
pixi run profile-linear  # macOS sample or Linux perf
```

The pinned Mojo compiler and all development dependencies are captured in
`pixi.lock`. Runtime and library code is Mojo-first and pure Mojo wherever
practical. Build-time data generation may use another language when justified,
but generated outputs must be deterministic, checksum-pinned, licensed, and
documented.

## Repository map

- `src/nagare/`: library source
- `tests/`: TestSuite unit, reference-value, and invariant tests
- `examples/`: small compilable usage programs
- `benchmarks/`: versioned benchmark, fixtures, and methodology
- `docs/`: architecture, design, compatibility, roadmap, and release policy
- `conda.recipe/`: local Rattler build recipe

See [the architecture](docs/architecture.md), [design principles](docs/design.md),
and [roadmap](docs/roadmap.md) before proposing a new dependency or feature.

## License

Licensed under either Apache-2.0 or MIT, at your option.
