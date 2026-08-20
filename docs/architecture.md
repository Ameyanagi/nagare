# Architecture

Nagare owns Interval search, interpolant construction and evaluation, spline coefficients, continuity, and extrapolation policy.

## Dependency boundary

Allowed ecosystem dependencies: Mojo standard library only for the initial release.
Expected downstream consumers: Scientific applications, data-processing libraries, and plotting clients that need interpolation.

Dependencies point from applications and higher-level packages toward smaller
foundations. This repository must never import a downstream consumer. New
dependencies require a documented need and must not force unrelated users to
install an application, renderer, language layer, or scientific stack.

## Layers

Planned implementation areas: search, linear and polynomial interpolation, PCHIP, Akima, cubic and B-splines, smoothing splines, and extrapolation policy.

The implemented foundation is layered as follows:

```text
public root
├── ExtrapolationPolicy       nominal out-of-domain behavior
├── locate_interval           validated standalone lookup
└── LinearInterpolator        validated owning interpolant
       ├── search             internal lookup over construction-validated data
       └── extrapolation      shared semantic policy
```

The unchecked interval locator is internal and operates on the finite,
strictly increasing knot invariant established at interpolant construction.
The standalone public locator still accepts raw user input and validates it on
each call.

The package root exports only the small documented public surface. Algorithms,
generated tables, platform details, and backend implementations remain in
their owning modules. Generic Mojo-native buffers, spans, strings, and
collections are preferred over an ecosystem-specific universal container.

## Data flow

Input validation occurs at the relevant public boundary: interpolant
construction validates stored tables, standalone search validates its raw knot
input, and evaluation validates each query. Internal layers operate on trusted
typed values and produce deterministic outputs for deterministic inputs. I/O,
clocks, randomness, terminal queries, filesystem access, and accelerator
selection stay at explicit effect or backend boundaries.

`LinearInterpolator` owns its `List[Float64]` inputs, but Mojo 1.0 does not make
underscore-prefixed fields private. Nagare treats those fields as private by
convention, and direct mutation is outside the contract. Read-only metadata is
`O(1)` and evaluation performs `O(log n)` interval search without rescanning
the table. Callers who perform unusual direct access can request an explicit
`O(n)` checkpoint through `validate()`.

Benchmark fixture generation and timing stay under `benchmarks/`; they are not
root exports and are not installed as library modules. The benchmark calls only
the documented root API, so it measures the same validation and arithmetic
contract available to downstream consumers.
