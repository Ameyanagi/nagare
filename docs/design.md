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
`LinearInterpolator.validate()` provides one explicit validation checkpoint
for callers performing unusual direct access without imposing an `O(n)` scan
on metadata access or repeated evaluation.

`ExtrapolationPolicy` follows the standard-library nominal-enum pattern: a
private-by-convention integer discriminant with `ERROR`, `CLAMP`, and `LINEAR`
constants as its public construction surface. This keeps policy selection
typed and makes equality one integer comparison.

Scaled coordinate arithmetic keeps interpolation finite across most of the
`Float64` range. On an interval spanning nearly `[-MAX_FINITE, +MAX_FINITE]`, a
small central offset can be below one ULP after normalization. The contract in
that ill-conditioned case is a finite result inside the endpoint hull, not
exact affine recovery.

## Natural cubic spline

`CubicSplineInterpolator` solves the natural tridiagonal system once at
construction with the Thomas algorithm. The endpoint second derivatives are
fixed to zero, and each interval is stored in four contiguous coefficient
buffers for a cubic in the shifted coordinate `dx = x - knots[i]`; queries do
not rebuild or revalidate that representation.

`ExtrapolationPolicy.LINEAR` always means a linear ray from an interpolant's
endpoint value with its endpoint first derivative. It never means extending a
boundary cubic. For a natural spline, the ray matches the value and first
derivative by construction and its zero second derivative matches the natural
boundary, so the join is C2-continuous. Under `CLAMP`, exterior derivatives are
those of the clamped constant: both the first and second derivative are zero.

## Out of scope

Plotting, dataframes, file I/O, optimization, general linear algebra, and signal-processing policy are outside this package.
