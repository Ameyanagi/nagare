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

if ! command -v git >/dev/null 2>&1; then
  printf '%s\n' 'git is required to identify benchmark source state' >&2
  exit 1
fi
benchmark_git_head=$(git rev-parse --verify HEAD)
if [[ -n $(git status --porcelain --untracked-files=normal) ]]; then
  benchmark_git_state=dirty
else
  benchmark_git_state=clean
fi

if command -v sha256sum >/dev/null 2>&1; then
  benchmark_lock_sha256=$(sha256sum pixi.lock | awk '{print $1}')
elif command -v shasum >/dev/null 2>&1; then
  benchmark_lock_sha256=$(shasum -a 256 pixi.lock | awk '{print $1}')
else
  printf '%s\n' 'sha256sum or shasum is required for benchmark metadata' >&2
  exit 1
fi

printf '%s\n' 'metadata_schema=nagare-linear-benchmark-metadata-v3'
printf 'run_utc=%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
printf 'git_head=%s\n' "$benchmark_git_head"
printf 'git_state=%s\n' "$benchmark_git_state"
printf 'pixi_lock_sha256=%s\n' "$benchmark_lock_sha256"
printf 'cpu=%s\n' "$benchmark_cpu"
printf 'os=%s\n' "$(uname -srv)"
printf 'architecture=%s\n' "$(uname -m)"
printf 'mojo=%s\n' "$(mojo --version)"
printf '%s\n' 'compiler_options=--optimization-level 3 -I src -I benchmarks'
printf '%s\n' 'command=pixi run bench-linear'

mojo run --optimization-level 3 -I src -I benchmarks benchmarks/bench_linear.mojo
