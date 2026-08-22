#!/usr/bin/env bash
set -euo pipefail

profile_dir=.pixi/profiles/nagare-linear
profile_binary="$profile_dir/profile_linear"
mkdir -p "$profile_dir"

mojo build --optimization-level 3 -I src -I benchmarks \
  benchmarks/profile_linear.mojo -o "$profile_binary"

profile_darwin() {
  local mode=$1
  "$profile_binary" "$mode" >"$profile_dir/$mode.out" &
  local profile_pid=$!
  /usr/bin/sample "$profile_pid" 3 \
    -file "$profile_dir/$mode.sample.txt" >/dev/null
  wait "$profile_pid"
}

profile_linux() {
  local mode=$1
  perf record --call-graph dwarf -o "$profile_dir/$mode.perf.data" -- \
    "$profile_binary" "$mode" >"$profile_dir/$mode.out"
  perf report --stdio -i "$profile_dir/$mode.perf.data" \
    >"$profile_dir/$mode.perf.txt"
}

case "$(uname -s)" in
  Darwin)
    if [[ ! -x /usr/bin/sample ]]; then
      printf '%s\n' 'macOS sample profiler is required' >&2
      exit 1
    fi
    profile_darwin scalar
    profile_darwin sorted
    profile_darwin sparse-scalar
    profile_darwin sparse-sorted
    ;;
  Linux)
    if ! command -v perf >/dev/null 2>&1; then
      printf '%s\n' 'Linux perf profiler is required' >&2
      exit 1
    fi
    profile_linux scalar
    profile_linux sorted
    profile_linux sparse-scalar
    profile_linux sparse-sorted
    ;;
  *)
    printf 'unsupported profiling platform: %s\n' "$(uname -s)" >&2
    exit 1
    ;;
esac

printf 'profiles=%s\n' "$profile_dir"
