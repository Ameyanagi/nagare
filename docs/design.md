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

Mojo 1.0 exposes underscore-prefixed fields to callers, so ownership alone
does not preserve a table invariant. Public observations revalidate stored
knots and values before indexing. This favors defined failure after external
mutation over an inaccurately advertised repeated-lookup cost.

`ExtrapolationPolicy` uses the three total states of `Optional[Bool]`. Even if a
caller mutates its exposed representation, every reachable state retains
defined error, clamp, or linear semantics.

Scaled coordinate arithmetic keeps interpolation finite across most of the
`Float64` range. On an interval spanning nearly `[-MAX_FINITE, +MAX_FINITE]`, a
small central offset can be below one ULP after normalization. The contract in
that ill-conditioned case is a finite result inside the endpoint hull, not
exact affine recovery.

## Out of scope

Plotting, dataframes, file I/O, optimization, general linear algebra, and signal-processing policy are outside this package.
