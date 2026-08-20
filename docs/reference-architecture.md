# Reference architecture for one-dimensional interpolation

This document records architecture research for Nagare's one-dimensional
interpolation stack. It is a design input, not an implementation port. No
upstream source, test vector, generated output, or dependency is copied into
Nagare. Mathematical identities will be implemented independently in pure
Mojo and tested from independently derived fixtures.

The v0.1 boundary remains the one in [the execution plan](v0.1-plan.md):
validated interval search, linear interpolation, and a natural cubic spline.
PCHIP and Akima are included here to ensure that the cubic foundation can
support them later without a rewrite; they remain explicit v0.1 non-goals.

## Research record

The repositories below were shallow-cloned outside the Nagare repository and
reviewed on 2026-08-20. The clone directories are research inputs only and are
not package, build, test, or benchmark dependencies.

| Project | Exact revision and version | License | API and algorithm surface reviewed |
| --- | --- | --- | --- |
| [enterpolation](https://github.com/NicolasKlenert/enterpolation) | `4d0cb2fdeb4ed7413a785968658c518a6cefc047`, crate `0.3.0` | `MIT OR Apache-2.0` in `Cargo.toml`, `LICENSE-MIT`, and `LICENSE-APACHE` | `src/linear/{mod,builder,error}.rs`, `src/base/{list,signal,adaptors,space}.rs`, and `src/bspline/{mod,builder,error}.rs` |
| [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) | `1564d003955a2c5aa84076a095e3cdf729fa9603`, package `0.16.3` | MIT in `LICENSE.md` | `src/Interpolations.jl`, `src/convenience-constructors.jl`, `src/gridded/`, `src/extrapolation/`, `src/monotonic/`, and `src/b-splines/` |
| [SciPy](https://github.com/scipy/scipy) | `a2c4d68b3cab98a0e5773ea0016e7c4109da6511`, source version `2.0.0.dev0` | BSD-3-Clause in `LICENSE.txt` | `scipy/interpolate/__init__.py`, `_cubic.py`, `_interpolate.py`, and the corresponding `tests/test_polyint.py` contracts |

These exact revisions can be recovered without relying on a moving branch:

```text
git clone --filter=blob:none --no-checkout <repository>
git -C <clone> fetch --depth=1 origin <revision>
git -C <clone> checkout --detach <revision>
```

### Primary-source observations

**enterpolation.** The crate root divides common signal/list/space machinery
from feature-gated linear, Bezier, and B-spline modules. Its
[`Linear::builder`](https://github.com/NicolasKlenert/enterpolation/blob/4d0cb2fdeb4ed7413a785968658c518a6cefc047/src/linear/mod.rs)
accepts elements plus explicit or equidistant knots and builds a scalar
`Signal::eval`/`Curve::domain` object. Sorted chains own upper-bound lookup;
`sample` and `extract` adapt iterators rather than add a separate array API.
Linear evaluation calculates one interval and merge factor per query.
[`BSpline`](https://github.com/NicolasKlenert/enterpolation/blob/4d0cb2fdeb4ed7413a785968658c518a6cefc047/src/bspline/mod.rs)
owns knots, elements, degree, and a static or dynamic workspace strategy, then
performs a de Boor-style local recurrence per query. Open, clamped, and legacy
builders report specific degree, count, ordering, and workspace failures. This
is strong evidence for separating construction from evaluation, but its
generic curve algebra and nonstandard B-spline knot modes are not a Nagare API
template.

**Interpolations.jl.** The
[`Interpolations` module](https://github.com/JuliaMath/Interpolations.jl/blob/1564d003955a2c5aa84076a095e3cdf729fa9603/src/Interpolations.jl)
layers `interpolate`, coordinate `scale`, and `extrapolate` wrappers over an
abstract-array hierarchy. Gridded interpolation owns knot vectors and data
coefficients, while cubic B-splines
[`prefilter`](https://github.com/JuliaMath/Interpolations.jl/blob/1564d003955a2c5aa84076a095e3cdf729fa9603/src/b-splines/prefiltering.jl)
the input at construction using tridiagonal or boundary-corrected systems.
The
[`MonotonicInterpolation`](https://github.com/JuliaMath/Interpolations.jl/blob/1564d003955a2c5aa84076a095e3cdf729fa9603/src/monotonic/monotonic.jl)
path is closer to Nagare's future PCHIP/Akima needs: method-specific tangent
construction emits owned linear, quadratic, and cubic segment coefficients,
and value/gradient/Hessian calls share local-polynomial evaluation. Its
[`Extrapolation`](https://github.com/JuliaMath/Interpolations.jl/blob/1564d003955a2c5aa84076a095e3cdf729fa9603/src/extrapolation/extrapolation.jl)
wrapper evaluates flat or linear continuation from a boundary value and
gradient. Nagare adopts those separations without the multidimensional type
and dependency hierarchy.

**SciPy.** The public
[`scipy.interpolate`](https://github.com/scipy/scipy/blob/a2c4d68b3cab98a0e5773ea0016e7c4109da6511/scipy/interpolate/__init__.py)
surface distinguishes high-level `CubicSpline`, `PchipInterpolator`, and
`Akima1DInterpolator` constructors from low-level `PPoly`, `BPoly`, and
`BSpline` representations. In
[`_cubic.py`](https://github.com/scipy/scipy/blob/a2c4d68b3cab98a0e5773ea0016e7c4109da6511/scipy/interpolate/_cubic.py),
one input-preparation path checks dimension, count, finiteness, and strict
ordering. PCHIP computes shape-preserving Fritsch-Butland tangents, Akima
computes local weighted tangents from extended secants, and cubic spline solves
for derivatives under several boundary modes; all three emit a piecewise
power-basis object. The
[`test_polyint.py`](https://github.com/scipy/scipy/blob/a2c4d68b3cab98a0e5773ea0016e7c4109da6511/scipy/interpolate/tests/test_polyint.py)
organization validates knot passage, derivative continuity, boundary
conditions, two-point reductions, shapes, and invalid inputs. Nagare adopts
the separation and categories of evidence, not SciPy's fixtures, array
semantics, broad boundary matrix, dependencies, or extrapolation defaults.

### What the references validate

| Concern | enterpolation | Interpolations.jl | SciPy | Nagare decision |
| --- | --- | --- | --- | --- |
| Construction | Consistent typed builders; named length, ordering, degree, knot, and workspace errors | Core interpolation plus scale and extrapolation wrappers; array-oriented constructors | High-level constructors validate through a shared `prepare_input` path | One direct constructor per algorithm over one internal validated sample contract |
| Sample validation | Builders wrap sorted knot chains, but non-decreasing knots and unchecked constructors are available | Checks axes and sorting; gridded construction can mutate duplicate knots with `nextfloat` | Cubic constructors require one-dimensional, finite, strictly increasing coordinates and finite values | Reject duplicates, descending knots, and all non-finite samples; never sort, deduplicate, or repair input |
| Query lookup | Sorted-chain upper-bound search, with specialized equidistant lookup | `searchsortedfirst`; specialized ordered-vector batch lookup exists | Shared piecewise-polynomial evaluation accepts scalar and array queries | Preserve the implemented right-biased binary search and exact-final-knot rule |
| Extrapolation | Linear extrapolation is the ordinary behavior; clamping is an adaptor | `extrapolate` is a wrapper with throw, flat, line, periodic, reflect, and directional variants | Policy varies by class; cubic piecewise-polynomial extrapolation continues an endpoint polynomial | Keep the closed `ERROR`, `CLAMP`, and `LINEAR` policy, with `ERROR` as the default |
| Scalar and batch | Scalar `eval`; iterator-based `sample` and `extract` compose batches | Scalar calls and vector arguments share the callable object | One call accepts scalar or array input | Scalar is the semantic kernel; a later explicit `evaluate_many` is a thin ordered loop |
| Cubic ownership | B-spline objects own knots, elements, degree, and workspace strategy; de Boor work is performed per query | Monotonic cubic objects own precomputed linear, quadratic, and cubic coefficients; B-splines prefilter data at construction | PCHIP, Akima, and cubic spline constructors produce a shared `PPoly` coefficient representation | Constructors precompute and own per-segment coefficients; evaluation never solves for them |
| Differentiation | Not a common linear/B-spline surface | Monotonic interpolation directly evaluates first and second derivatives | `PPoly` owns derivative and antiderivative operations | v0.1 cubic exposes first and second derivatives needed to verify its C2 contract; broader polynomial algebra waits |

The comparison is architectural, not an assertion that the projects have
identical semantics. In particular, enterpolation's B-spline is a generic
control-point curve, Interpolations.jl's cubic B-spline prefilters sampled
values, and SciPy's `CubicSpline` constructs a data-interpolating piecewise
polynomial. Nagare must not hide those different mathematical objects behind
one ambiguous `spline` mode.

## Validated sample table

All one-dimensional data interpolators share an internal owning sample-table
contract:

```text
_ValidatedSamples1D
├── knots:  List[Float64]
└── values: List[Float64]
```

It is not a root export and is not a new Nagare array type. Its constructor
takes ownership only after validating:

- at least two knots;
- equal knot and value counts;
- finite knots and values; and
- strictly increasing knots.

The contract deliberately rejects the two repair behaviors seen in broader
libraries: automatic sorting and duplicate displacement. Those operations can
change which observation belongs to which coordinate and cannot be a silent
property of interpolation construction.

Mojo 1.0 does not make underscore-prefixed fields private. A caller can reach
and mutate an interpolator's stored lists after construction, even if the type
is intended to be owning. Every public observation therefore revalidates the
sample table before indexing it. A prepared cubic also validates coefficient
counts, finiteness, knot passage, continuity, boundary conditions, and the
algorithm-specific coefficient relation. This cross-state validation rejects
a sample mutation that would otherwise leave stale coefficients, not just a
shape mutation that could make indexing unsafe. Coordinated mutation that
still satisfies every invariant describes another valid interpolator, so its
semantics remain defined. Until Mojo provides an enforceable encapsulation
mechanism, scalar evaluation is consequently `O(n)` validation plus `O(log n)`
search rather than advertised as `O(log n)` alone.

The first internal refactor should make linear and cubic construction call one
sample validation function. It must not expose `_ValidatedSamples1D` from
`src/nagare/__init__.mojo` or make public callers unwrap an internal table just
to construct an interpolator.

## Query and extrapolation contract

### Interval selection

All algorithms use the existing interval definition and locator:

- interval `i` means `[knots[i], knots[i + 1]]`;
- an exact interior knot selects the interval on its right;
- the final knot selects the final interval; and
- finite out-of-domain queries are classified before the internal locator is
  called.

An exact-knot fast path returns the original tabulated ordinate. It is both a
user-visible exactness guarantee and protection against a final rounding step
in a coefficient polynomial.

Specialized uniform-grid or monotonically ordered batch search may be added
only behind the same interval contract. It must not be encoded in a generic
curve trait or exposed as a second definition of knot boundaries.

### Extrapolation

`ExtrapolationPolicy` remains the single owner of out-of-domain behavior:

| Policy | Below the domain | Above the domain |
| --- | --- | --- |
| `ERROR` | raise | raise |
| `CLAMP` | return the first ordinate | return the final ordinate |
| `LINEAR` | extend the value at the first knot using its right derivative | extend the value at the final knot using its left derivative |

For a cubic interpolator, `LINEAR` means endpoint-tangent continuation. It does
not continue the first or last cubic polynomial outside its segment. This
preserves the meaning of the policy across linear, natural cubic, PCHIP, and
Akima interpolators and avoids cubic growth that a caller did not request.

Non-finite queries always raise before policy dispatch. Finite linear
extrapolation may return IEEE signed infinity when the represented endpoint
line is outside the finite `Float64` range, but it must never create NaN from
finite state.

## Scalar and batch evaluation

The scalar operation is the canonical semantic and testing surface:

```mojo
var value = interpolator.evaluate(query)
```

After the scalar implementation of an algorithm is stable, an explicit batch
operation may be added:

```mojo
var values = interpolator.evaluate_many(queries)
```

`evaluate_many` will:

- preserve query order and return an owned `List[Float64]`;
- be exactly equivalent to repeated scalar calls;
- fail at the first invalid query without returning a partial result; and
- avoid scalar/collection overloading or a Nagare-specific array abstraction.

The first implementation is a thin loop. Reusing interval state for sorted
queries, SIMD, parallel execution, output buffers, and multi-dimensional
broadcasting require separate evidence and are not part of the initial batch
contract.

## Algorithm layering

The stable dependency direction is:

```text
root API
├── LinearInterpolator
│     ├── validated samples
│     ├── interval search
│     └── extrapolation policy
│
├── NaturalCubicSpline                         v0.1
│     ├── natural tridiagonal construction
│     └── _PiecewiseCubic1D
│           ├── validated samples
│           ├── interval search
│           └── extrapolation policy
│
├── PchipInterpolator                         after v0.1
│     ├── shape-preserving tangent construction
│     ├── Hermite coefficient construction
│     └── _PiecewiseCubic1D
│
└── AkimaInterpolator                         after PCHIP
      ├── local Akima tangent construction
      ├── Hermite coefficient construction
      └── _PiecewiseCubic1D
```

`_PiecewiseCubic1D` is an internal evaluated representation, analogous in role
but not implementation to SciPy's `PPoly` and Interpolations.jl's monotonic
coefficient object. For segment `i`, it stores four coefficients for Horner
evaluation in a dimensionless local coordinate `t`:

```text
t = (x - knots[i]) / (knots[i + 1] - knots[i])
p_i(t) = ((d_i * t + c_i) * t + b_i) * t + a_i
```

The stable segment-parameter routine already used by linear interpolation is
the sole owner of overflow-resistant `t` calculation. First- and
second-derivative evaluation applies the interval scale explicitly; value
evaluation does not form powers of a potentially extreme knot difference.

PCHIP and Akima differ only in how knot tangents are chosen. Both feed the
same independently implemented Hermite-to-power coefficient builder. Natural
cubic construction instead performs a global tridiagonal solve, then emits the
same normalized per-segment representation. This split gives coefficient
construction and coefficient evaluation separate tests and benchmarks.

B-splines do not belong underneath `_PiecewiseCubic1D`. A future control-point
B-spline needs knot multiplicity, degree, and de Boor workspace semantics; a
future interpolating B-spline needs a coefficient prefilter or solve. Those are
separate designs after the data-interpolating v0.1 proof.

## Precomputation and ownership

Each algorithm constructor owns all state required by later evaluation:

| Algorithm | Construction | Stored state | Scalar observation |
| --- | --- | --- | --- |
| Linear | validate and take ownership, `O(n)` | samples and policy | revalidate `O(n)`, locate `O(log n)`, evaluate one segment |
| Natural cubic | validate, scale, solve one tridiagonal system, and form coefficients, `O(n)` | samples, four coefficient arrays, endpoint derivatives, scale state, and policy | revalidate samples and cross-state cubic invariants `O(n)`, locate `O(log n)`, Horner evaluation |
| PCHIP | validate, compute secants and shape-preserving tangents, form Hermite coefficients, `O(n)` | same evaluated state as natural cubic | same cubic evaluator |
| Akima | validate, extend local secants, compute weighted tangents, form Hermite coefficients, `O(n)` | same evaluated state as natural cubic | same cubic evaluator |

No solve, coefficient allocation, or scratch-workspace allocation occurs in
scalar evaluation. Construction and evaluation remain separate benchmark
families. A future enforceably sealed representation may allow construction
validation to be amortized; current Mojo behavior does not.

Coefficient arrays have exactly `knot_count - 1` elements each. Endpoint
derivatives are stored or deterministically recoverable from the endpoint
coefficients so `LINEAR` extrapolation does not call algorithm-specific
construction logic.

## Errors and numerical conditioning

### Error ownership

Construction raises before producing an interpolator for:

- an invalid sample table;
- an algorithm-specific minimum-size failure;
- a zero or non-finite normalized interval width;
- an unusable pivot or non-finite recurrence during the natural tridiagonal
  solve; or
- any non-finite stored scale or coefficient.

Observation raises for a non-finite query, `ERROR` extrapolation, or externally
mutated sample/coefficient state that violates structural, numerical, or
algorithm-specific cross-state invariants. Natural cubic validation checks
knot passage, C1/C2 joins, and natural boundaries. PCHIP and Akima validation
recomputes their local tangent relation before trusting exposed coefficient
storage. No constructor returns a partially prepared object, and no invalid
input is repaired.

Natural cubic accepts two samples as the unique straight line satisfying the
natural endpoint conditions. PCHIP also accepts two samples and reduces to
that line. Classical Akima can synthesize its two exterior secants, but the
first Nagare Akima issue must choose and test one explicit minimum; it must not
obtain a small-sample behavior accidentally through out-of-bounds indexing.

### Conditioning strategy

Finite inputs do not imply that every intermediate in a naive formula is
representable. Cubic work therefore follows these gates:

1. Compute each local parameter with the existing scaled-coordinate path;
   never require `x1 - x0` to be finite merely to evaluate a segment.
2. Construct with dimensionless knot widths derived from one table-wide
   positive coordinate scale. Use direct subtraction when representable and a
   consistently scaled subtraction when it is not. Reject a width that becomes
   zero at working precision rather than divide by it.
3. Normalize ordinates by one finite positive table scale before building or
   solving for cubic coefficients. Store the output scale with dimensionless
   coefficients instead of forming avoidably overflowing ordinate deltas or
   powers of interval widths.
4. Scale each tridiagonal row before elimination, verify every pivot before
   division, and verify the residual of the solved natural-spline system.
5. Form PCHIP weighted harmonic tangents without multiplying two large
   secants. Sign changes and zero secants select a zero tangent before any
   reciprocal calculation.
6. Form Akima weights with scaled absolute differences. When both weights are
   zero, use the documented symmetric average rather than a `0 / 0` result.
7. Validate every prepared coefficient and scale. A numerically unusable table
   is an explicit construction error, not a NaN-bearing interpolator.

The exact residual tolerance and scaling equations are deliverables of the
coefficient-construction issue, recorded before implementation fixtures are
accepted. They may not be widened solely because an upstream implementation
accepts a case.

For an in-domain query with a finite mathematical result, evaluation should
avoid intermediate overflow. As already documented for linear interpolation,
tables spanning nearly the entire `Float64` range can lose a small offset after
normalization. Tests in that region require finite, bounded, deterministic
behavior; they do not promise exact affine recovery below available precision.

## Adopted and rejected ideas

| Reference idea | Decision | Reason |
| --- | --- | --- |
| enterpolation's consistent construction/evaluation split and named build failures | Adopt | It makes invalid configuration fail before evaluation and keeps errors attributable |
| enterpolation's specialized sorted/equidistant lookup | Defer behind the common locator | It is useful only after a workload proves the additional representation is needed |
| enterpolation's generic `Signal`/`Curve`/chain/adaptor algebra | Reject for the initial API | It expands the abstraction surface beyond precise scalar numeric interpolation |
| enterpolation's ordinary extrapolation and unchecked constructors | Reject | Nagare defaults to an error and preserves validated public state |
| Interpolations.jl's separate interpolation and extrapolation responsibilities | Adopt semantically | Nagare uses one closed policy value rather than a large wrapper hierarchy |
| Interpolations.jl's precomputed monotonic cubic coefficients and shared local polynomial evaluator | Adopt | PCHIP and Akima should vary tangent construction, not duplicate evaluation |
| Interpolations.jl's duplicate-knot displacement | Reject | Moving a sample coordinate silently changes user data and may create extreme conditioning |
| Interpolations.jl's multidimensional abstract-array and weighted-index hierarchy | Reject for this architecture | It introduces data-container, broadcasting, and dimensional policy outside the one-dimensional mission |
| SciPy's shared strict finite/increasing validation for cubic constructors | Adopt | It matches Nagare's existing table contract |
| SciPy's high-level PCHIP, Akima, and cubic constructors over one piecewise-polynomial representation | Adopt in narrower form | It provides the desired coefficient/evaluation separation without exposing general polynomial algebra |
| SciPy's class-dependent default extrapolation and endpoint-polynomial continuation | Reject | One policy must mean the same operation for every Nagare interpolator |
| SciPy's axis, N-dimensional ordinate, complex, and array-namespace behavior | Reject for the initial API | Nagare v0.1 is one-dimensional `Float64` interpolation with no framework array dependency |
| Any upstream test output as a golden fixture | Reject | Reference values must be independently derived; upstream libraries may only cross-check them |

## Minimal Mojo API

The package root remains class-oriented and small:

```mojo
from nagare import (
    ExtrapolationPolicy,
    LinearInterpolator,
    NaturalCubicSpline,
    locate_interval,
)

var spline = NaturalCubicSpline(
    knots^,
    values^,
    extrapolation=ExtrapolationPolicy.ERROR,
)
var y = spline.evaluate(1.25)
var dy = spline.derivative(1.25)
var d2y = spline.second_derivative(1.25)
```

After v0.1, PCHIP and Akima add only their constructor names and the same
observation methods:

```mojo
from nagare import AkimaInterpolator, PchipInterpolator
```

There is no public `_PiecewiseCubic1D`, coefficient table, method enum,
builder/director typestate graph, general `Spline`, generic output vector
space, `interp1d` factory, or universal array. `NaturalCubicSpline` names its
boundary condition explicitly; additional cubic boundary conditions require a
new reviewed public contract rather than a stringly typed option.

`evaluate_many` is not a prerequisite for any scalar constructor. When added,
it appears consistently on each interpolator rather than as a second root
factory.

## Verification corpus

### Reference-value tests

Every checked-in expected table is small, hand-derived or generated from an
independent mathematical script whose equations and provenance are recorded.
SciPy or Interpolations.jl may be used as a secondary cross-check, never as the
sole oracle.

- **Shared:** exact knots; constant and affine data; irregular spacing; two
  samples; each extrapolation policy; non-finite rejection; and externally
  mutated sample/coefficient state.
- **Natural cubic:** a hand-solved three- and four-knot system, tridiagonal
  equation residuals, natural endpoint second derivatives, and independently
  evaluated interior points.
- **PCHIP:** two-point linear reduction, zero tangents at a plateau or secant
  sign change, one-sided endpoint selection, and weighted-harmonic interior
  tangents.
- **Akima:** extended endpoint secants, a nonzero-weight local example, an
  equal-weight fallback, and a case that distinguishes classical Akima from
  modified Akima. The initial type implements only classical Akima.

### Properties and invariants

- all interpolators reproduce every stored ordinate exactly;
- natural cubic is C0, C1, and C2 at every interior knot and has zero endpoint
  second derivatives;
- PCHIP is C1 and does not overshoot the endpoint hull on each monotone
  segment;
- Akima is C1, preserves constant and affine data, and changing a distant
  sample outside the local stencil does not change an unaffected tangent;
- batch output equals repeated scalar evaluation exactly under the same query
  order;
- translating or positively scaling well-conditioned coordinates and values
  produces the correspondingly transformed interpolant within a declared
  tolerance; and
- deterministic extreme tables cover subnormal spacing, adjacent large
  values, mixed-sign maximum-scale coordinates, uneven interval ratios,
  overflow-directed linear extrapolation, and mutation after construction.

Continuity tests evaluate the left and right coefficient formulas at the knot;
they do not compare two calls that both select the same right-biased interval.

### Benchmark matrix

Benchmarks extend the current checked-manifest methodology without committing
raw timing output or comparative claims:

| Dimension | Cases |
| --- | --- |
| Algorithm | linear, natural cubic, then PCHIP and Akima when released |
| Phase | construction only; scalar evaluation of an already constructed object; batch wrapper later |
| Table size | small, medium, large |
| Spacing | uniform, irregular, highly uneven, mixed-sign extreme |
| Query pattern | exact knots, interior fixed, deterministic varying, left/right extrapolation by policy |
| Semantic retention | nonzero coefficient/output checksum checked against the versioned manifest |

Construction must not occur inside an evaluation timing loop. Evaluation must
not include fixture generation. Fixed queries use the repository's black-box
or deterministic-varying-query discipline so loop-invariant computation and
dead-code elimination cannot erase the measured work. Run metadata continues
to record the Mojo toolchain, flags, host, Git revision and dirty state, and
`pixi.lock` SHA-256.

## Dependency-ordered issues

The first four issues below are the completed/active v0.1 linear foundation.
Each later item is small enough for its own review and inherits all repository
validation gates from `docs/v0.1-plan.md`.

1. **NAG-001 through NAG-004 — contracts, interval search, linear
   interpolation, and baseline benchmark.** Keep their existing behavior as
   the compatibility foundation.
2. **NAG-RA-001 — reference architecture.** Record exact upstream revisions,
   adopted/rejected decisions, the evaluated cubic representation, numerical
   gates, corpus, and issue ordering. Gate: documentation review plus locked
   repository and package checks.
3. **NAG-005A — shared validated samples and cubic representation.** Extract
   shared table validation without changing root API; add internal
   `_PiecewiseCubic1D`, structural/mutation validation, exact-knot evaluation,
   and first/second derivative formula tests. Gate: hand-authored coefficients
   exercise every segment and policy without a spline solver; stale samples or
   coefficients fail before indexing or evaluation.
4. **NAG-005B — natural cubic system construction.** Implement dimensionless
   width/value scaling and a pure-Mojo tridiagonal solver with pivot, finite,
   and residual gates. Gate: hand-solved systems, affine reduction, invalid and
   unusable-system tests.
5. **NAG-006 — natural cubic public evaluation.** Export
   `NaturalCubicSpline`; connect construction to the internal cubic evaluator
   and all three policies. Gate: exact knots, C0/C1/C2, natural boundaries,
   reference values, extreme finite cases, compile-fail/API tests, and package
   smoke.
6. **NAG-007 — v0.1 downstream proof.** Complete the existing supported CI,
   installed-package, and independent downstream gates. This closes v0.1.
7. **NAG-008 — PCHIP tangent construction.** Independently implement the
   Fritsch-Butland rule and endpoint limiter with scaled arithmetic; emit the
   shared Hermite coefficients. Gate: tangent fixtures, monotone no-overshoot,
   C1, two-point, plateau, sign-change, mutation, and package smoke tests.
8. **NAG-009 — classical Akima tangent construction.** Define the supported
   sample minimum and endpoint extension, implement scaled local weights, and
   reuse the Hermite builder. Gate: local-stencil, equal-weight, affine, C1,
   extreme, mutation, and package smoke tests.
9. **NAG-010 — explicit batch evaluation.** Add `evaluate_many` consistently
   after all scalar contracts are stable. Gate: exact scalar equivalence,
   order preservation, first-error behavior, empty input, and construction vs
   batch benchmark separation.
10. **NAG-011 — cubic benchmark extension.** Add versioned natural cubic,
    PCHIP, and Akima construction/evaluation cases only where the matrix makes
    an architectural decision. Gate: checked manifest shape and checksums,
    optimizer-retention audit, complete run metadata, and no raw timing or
    marketing claim.

PCHIP and Akima issues do not move into v0.1 merely because their architecture
is now known. B-splines, smoothing splines, polynomial interpolation,
multidimensional interpolation, plotting, dataframes, file I/O, optimization,
general linear algebra, and runtime Python/C/Rust dependencies remain outside
this issue sequence.
