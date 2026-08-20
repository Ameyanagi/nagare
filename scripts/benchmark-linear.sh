#!/usr/bin/env bash
set -euo pipefail

case "$(uname -s)" in
  Darwin)
    benchmark_cpu=$(sysctl -n machdep.cpu.brand_string 2>/dev/null || uname -m)
    ;;
  Linux)
    benchmark_cpu=""
    if command -v lscpu >/dev/null 2>&1; then
      benchmark_cpu=$(lscpu | awk -F: '/Model name/ {sub(/^[ \t]+/, "", $2); print $2; exit}')
    fi
    if [[ -z "$benchmark_cpu" ]]; then
      benchmark_cpu=$(uname -m)
    fi
    ;;
  *)
    benchmark_cpu=$(uname -m)
    ;;
esac

printf '%s\n' 'metadata_schema=nagare-linear-benchmark-metadata-v1'
printf 'run_utc=%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
printf 'cpu=%s\n' "$benchmark_cpu"
printf 'os=%s\n' "$(uname -srv)"
printf 'architecture=%s\n' "$(uname -m)"
printf 'mojo=%s\n' "$(mojo --version)"
printf '%s\n' 'compiler_options=--optimization-level 3 -I src -I benchmarks'
printf '%s\n' 'command=pixi run bench-linear'

mojo run --optimization-level 3 -I src -I benchmarks benchmarks/bench_linear.mojo
