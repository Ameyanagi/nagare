# Roadmap

## v0.1 — Foundation

- Implement validated knot preparation, binary-search interval location, linear interpolation, and cubic splines with independently generated reference fixtures.
- Define the smallest useful public API and its invariants.
- Add unit, reference-value, and property/invariant coverage.
- Build and test the precompiled package on supported targets.

The issue-sized sequence, numerical contracts, validation gates, and explicit
non-goals are maintained in [the v0.1 execution plan](v0.1-plan.md).

## v0.2 — Usability

- Add ergonomic APIs only after v0.1 usage demonstrates repeated friction.
- Expand examples and integration fixtures.
- Publish the first modular-community recipe when the package is useful alone.

## v0.3 — Performance

- Extend the checked-in linear baseline only when a representative new workload
  exposes a decision the existing matrix cannot answer.
- Optimize measured bottlenecks without weakening correctness or API clarity.
- Add SIMD or specialized backends only behind the same semantic contract.

## v1.0 — Stability

- Document every public symbol and error contract.
- Provide a compatibility and deprecation policy.
- Support the declared OS and architecture matrix in CI.
- Require downstream proof from at least one independent consumer.

## Not planned

Plotting, dataframes, file I/O, optimization, general linear algebra, and signal-processing policy are outside this package.
