# Data provenance

No generated lookup data is currently committed.

The small linear-interpolation fixtures in `tests/test_linear.mojo` are
independently calculated from the documented two-point line equation. They do
not copy outputs from SciPy or another implementation. The irregular fixture
uses knots `[0, 0.5, 2, 5]`, values `[1, 2, -1, 8]`, and exact hand-computed
queries chosen to exercise every segment.

The benchmark tables are generated deterministically in memory by committed
pure Mojo code. Their formulas, expected sizes, extreme cases, and absence of
external licensing requirements are recorded in
`benchmarks/datasets/README.md`. Result checksums guard semantics during timing;
they are not presented as cryptographic source checksums.

Every future generated artifact must record:

- upstream project and canonical URL;
- upstream version and retrieval date;
- exact file checksums and licenses;
- generator source and command;
- deterministic output checks;
- review notes for semantic or licensing changes.

Generation tools are development dependencies. Consumers install the generated
Mojo data and do not require Python, Rust, C, or another runtime.
