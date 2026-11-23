#!/usr/bin/env bash
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 OUTPUT_DIR" >&2
  exit 2
fi
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$1"
output="$(cd "$1" && pwd)"
make -C "$root/crates/spatial-stata" PLUGIN_OUT="$output/spatial_spur.plugin"
cp "$root/stata/"*.ado "$root/stata/"*.pkg "$root/stata/stata.toc" "$output/"
for help in "$root"/stata/*.sthlp; do
  if [[ -f "$help" ]]; then cp "$help" "$output/"; fi
done
cp "$root/LICENSE" "$output/LICENSE"
cp "$root/crates/spatial-py/THIRD_PARTY_NOTICES.txt" "$output/THIRD_PARTY_NOTICES.txt"
cp "$root/crates/spatial-stata/c/SPI-NOTICE" "$output/SPI-NOTICE"
