#!/usr/bin/env bash
set -euo pipefail

for test_file in tests/test_*.mojo; do
  mojo run -I src -I benchmarks "$test_file"
done

mkdir -p .pixi/test-bin

assert_compile_failure() {
  local source_file=$1
  local expected_diagnostic=$2
  local fixture_name
  fixture_name=$(basename "$source_file" .mojo)
  local diagnostic_file=".pixi/test-bin/${fixture_name}.log"

  if mojo build -I src "$source_file" -o ".pixi/test-bin/${fixture_name}" \
    >"$diagnostic_file" 2>&1; then
    printf 'Expected compilation to fail: %s\n' "$source_file" >&2
    exit 1
  fi
  local diagnostic
  diagnostic=$(<"$diagnostic_file")
  if [[ "$diagnostic" != *"$expected_diagnostic"* ]]; then
    printf 'Compile-fail fixture produced an unexpected diagnostic: %s\n' \
      "$source_file" >&2
    sed -n '1,120p' "$diagnostic_file" >&2
    exit 1
  fi
}

assert_compile_failure \
  tests/compile_fail/integer_extrapolation_policy.mojo \
  "missing required argument: '_value'"
assert_compile_failure \
  tests/compile_fail/validated_extrapolation_bypass.mojo \
  "missing required argument: '_value'"

mojo build -I src examples/basic.mojo -o .pixi/test-bin/basic
mojo build -I src examples/resample_sensor.mojo -o .pixi/test-bin/resample_sensor

benchmark_manifest=$(
  mojo run --optimization-level 3 -I src -I benchmarks \
    benchmarks/bench_linear.mojo --manifest
)
if ! diff -u benchmarks/linear_manifest.txt \
  <(printf '%s\n' "$benchmark_manifest"); then
  printf '%s\n' 'Benchmark manifest changed; review its schema and methodology.' >&2
  exit 1
fi
