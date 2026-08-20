# Changelog

This project follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and uses semantic versioning after the first public release.

## [Unreleased]

### Added

- Initial experimental repository scaffold.
- A validated, right-biased binary-search interval locator.
- An owning linear interpolator with explicit error, clamp, and linear
  extrapolation policies.
- Direct-delta and scale-normalized affine evaluation paths that avoid
  premature overflow for finite in-domain data, preserve tiny-ordinate
  extrapolation when representable, and define signed-infinity extrapolation
  overflow.
- Observation-time validation of externally mutable knot and value storage,
  plus a bounded-result contract for ill-conditioned full-range intervals.
- Exact, independently calculated reference, and property/invariant tests for
  the first numerical vertical slice.
- An issue-sized v0.1 execution plan with numerical contracts and validation
  gates.
- A versioned 26-case linear construction/evaluation benchmark baseline with
  deterministic fixtures, host/compiler metadata, semantic checksums, and no
  comparative performance claim.
