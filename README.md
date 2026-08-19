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

The current scaffold includes only an internal smoke marker. Nothing is
re-exported as a stable public API yet.

## Repository map

- `src/nagare/`: library or application source
- `tests/`: TestSuite unit, reference-value, and invariant tests
- `examples/`: small compilable usage programs
- `benchmarks/`: reproducible methodology and later benchmark programs
- `docs/`: architecture, design, compatibility, roadmap, and release policy
- `conda.recipe/`: local Rattler build recipe

See [the architecture](docs/architecture.md), [design principles](docs/design.md),
and [roadmap](docs/roadmap.md) before proposing a new dependency or feature.

## License

Licensed under either Apache-2.0 or MIT, at your option.
