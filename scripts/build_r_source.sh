#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 OUTPUT_DIR" >&2
  exit 2
fi

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$1"
output="$(cd "$1" && pwd)"
stage="$(mktemp -d)"
trap 'rm -rf "$stage"' EXIT
package="$stage/spatialSpur"
mkdir -p "$package"
rsync -a --exclude='/src/rust-target' --exclude='/src/*.o' --exclude='/src/*.so' \
  --exclude='/src/*.dll' --exclude='*.DS_Store' "$root/r/spatialSpur/" "$package/"
rust="$package/src/rust"
mkdir -p "$rust/crates"
cp -R "$root/crates/spatial-core" "$root/crates/spatial-r" "$rust/crates/"
cp "$root/Cargo.lock" "$rust/Cargo.lock"
cp "$root/LICENSE" "$rust/LICENSE"
cat > "$rust/Cargo.toml" <<'EOF'
[workspace]
resolver = "2"
members = ["crates/spatial-core", "crates/spatial-r"]

[profile.release]
opt-level = 3
lto = "thin"
codegen-units = 1
EOF
(
  cd "$rust"
  # Cargo removes unrelated workspace members from the copied lockfile.
  cargo metadata --format-version 1 >/dev/null
  cargo vendor --locked vendor > vendor-config.toml
  COPYFILE_DISABLE=1 tar -cf - vendor | xz -9e > vendor.tar.xz
  rm -rf vendor
)
(
  cd "$output"
  R CMD build --no-build-vignettes "$package"
)
version="$(sed -n 's/^Version: //p' "$package/DESCRIPTION")"
archive="$output/spatialSpur_$version.tar.gz"
bytes="$(wc -c < "$archive" | tr -d '[:space:]')"
echo "Built $archive ($bytes bytes)"
if [[ "$bytes" -ge 10000000 ]]; then
  echo "The package exceeds the project's 10 MB source archive budget." >&2
  exit 1
fi
