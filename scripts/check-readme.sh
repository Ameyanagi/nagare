#!/usr/bin/env bash
set -euo pipefail

blocks_dir=.pixi/readme-blocks
mkdir -p "$blocks_dir"
rm -f "$blocks_dir"/block_*

awk -v blocks_dir="$blocks_dir" '
  $0 ~ /^```mojo[[:space:]]*$/ {
    block_count++
    in_block = 1
    output_file = sprintf("%s/block_%d.mojo", blocks_dir, block_count)
    printf "%s", "" > output_file
    close(output_file)
    next
  }

  in_block && $0 ~ /^```[[:space:]]*$/ {
    close(output_file)
    in_block = 0
    next
  }

  in_block {
    print > output_file
  }

  END {
    if (in_block) {
      print "Unclosed fenced Mojo block in README.md" > "/dev/stderr"
      exit 1
    }
  }
' README.md

shopt -s nullglob
readme_blocks=("$blocks_dir"/block_*.mojo)
if ((${#readme_blocks[@]} == 0)); then
  printf '%s\n' 'README.md contains no fenced Mojo blocks.' >&2
  exit 1
fi

for source_file in "${readme_blocks[@]}"; do
  block_name=$(basename "$source_file" .mojo)
  diagnostic_file="$blocks_dir/${block_name}.log"
  output_file="$blocks_dir/$block_name"

  if ! mojo build -I src "$source_file" -o "$output_file" \
    >"$diagnostic_file" 2>&1; then
    printf 'README Mojo block failed to compile: %s\n' "$source_file" >&2
    cat "$diagnostic_file" >&2
    exit 1
  fi
done
