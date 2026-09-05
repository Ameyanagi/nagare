# Design

## Principles

- Mojo is the runtime implementation language.
- Prefer pure Mojo and safe standard-library APIs.
- Keep the root API small, typed, documented, and testable.
- Separate semantic contracts from optimized CPU, SIMD, GPU, terminal, or
  rendering backends.
- Establish correctness and reference fixtures before optimization.
- Make invalid public configuration unrepresentable when practical; otherwise
  reject it explicitly.
- Preserve source mappings, numerical tolerances, ownership, and provenance as
  first-class data when the domain requires them.
- Do not add a framework-wide array, executor, renderer, or application model.

## Tradeoffs

The project accepts a narrower initial feature set in exchange for reviewable
contracts and sparse dependencies. Generated tables are acceptable when their
sources, Unicode or data version, licenses, checksums, and deterministic update
procedure are committed. Consumers must not need the generator toolchain.

Construction establishes table invariants, and read-only methods trust them
thereafter. Mojo 1.0 exposes underscore-prefixed fields to callers, but Nagare
treats them as private by convention; direct mutation is outside the contract.
Each owning interpolator provides one explicit `validate()` checkpoint for
callers performing unusual direct access without imposing an `O(n)` scan on
metadata access or repeated evaluation.

`ExtrapolationPolicy` follows the standard-library nominal-enum pattern: a
private-by-convention integer discriminant with `ERROR`, `CLAMP`, `LINEAR`, and
payload-carrying `FILL` as its public construction surface. This keeps policy
selection typed while preserving exact payload identity.

Scaled coordinate arithmetic keeps interpolation finite across most of the
`Float64` range. On an interval spanning nearly `[-MAX_FINITE, +MAX_FINITE]`, a
small central offset can be below one ULP after normalization. The contract in
that ill-conditioned case is a finite result inside the endpoint hull, not
exact affine recovery.

## Custom interval lookup

`KnotIndex` is the small owned counterpart to `locate_interval`. It validates
finite, strictly increasing knots at construction and trusts its private-by-
convention storage thereafter. `validate()` provides an explicit checkpoint.
Metadata and the read-only `knots()` view are non-raising and O(1); `locate(x)`
validates the new finite, in-domain query and uses the same right-biased binary
search as every interpolant. Construction costs O(n); Q queries cost
O(Q log n), avoiding the one-shot API's O(Q n) knot validation.

`locate_into` accepts spans and writes into caller-owned integer storage without
allocating. It supports unsorted and repeated queries, checks matching buffer
lengths, and has the same error contract as interpolation's `evaluate_into`:
results are unspecified after an invalid query raises.

## Cubic spline boundaries

`CubicSplineInterpolator` defaults to scipy-compatible not-a-knot boundaries.
It also supports natural, clamped, and periodic conditions through the nominal
`BoundaryCondition` value. Natural and clamped systems use Thomas elimination;
not-a-knot removes its two corner terms before Thomas elimination, and periodic
uses a Sherman-Morrison cyclic solve. The scipy degradations are explicit:
not-a-knot over two knots is a line, while three knots form the single parabola
through the samples.

Each interval is stored in four contiguous coefficient buffers for a cubic in
the shifted coordinate `dx = x - knots[i]`; queries and closed-form integration
do not rebuild or revalidate that representation. All boundary choices produce
a C2 spline, and construction is O(n).

`ExtrapolationPolicy.LINEAR` always means a linear ray from an interpolant's
endpoint value with its endpoint first derivative. It never means extending a
boundary cubic. The ray matches the value and first derivative at the join for
every boundary condition. Its zero second derivative makes that join
C2-continuous only for a natural spline. Under `CLAMP`, exterior derivatives
are those of the clamped constant: both first and second derivative are zero.

## Out of scope

Plotting, dataframes, file I/O, optimization, general linear algebra, and signal-processing policy are outside this package.
