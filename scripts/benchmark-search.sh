#!/usr/bin/env bash
set -euo pipefail

benchmark_binary=.pixi/bench_search
mojo build --optimization-level 3 -I src benchmarks/bench_search.mojo -o "$benchmark_binary"
"$benchmark_binary"
