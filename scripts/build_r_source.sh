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
rsync -a --exclude='/src/rust-target' --exclude='/src/rust/target' --exclude='/src/*.o' --exclude='/src/*.so' \
  --exclude='/src/*.dll' --exclude='*.DS_Store' "$root/r/spatialSpur/" "$package/"
rust="$package/src/rust"
mkdir -p "$rust/crates"
rsync -a --exclude='/tests' "$root/crates/spatial-core/" "$rust/crates/spatial-core/"
rsync -a --exclude='/tests' "$root/crates/spatial-r/" "$rust/crates/spatial-r/"
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

[workspace.metadata.vendor-filter]
platforms = [
    "x86_64-unknown-linux-gnu",
    "aarch64-unknown-linux-gnu",
    "x86_64-apple-darwin",
    "aarch64-apple-darwin",
    "x86_64-pc-windows-gnu",
    "aarch64-pc-windows-gnullvm",
    "x86_64-pc-windows-msvc",
]
exclude-crate-paths = [
    { name = "*", exclude = "tests" },
    { name = "*", exclude = "benches" },
    { name = "*", exclude = "examples" },
]
EOF
(
  cd "$rust"
  # Cargo removes unrelated workspace members from the copied lockfile.
  cargo metadata --format-version 1 >/dev/null
  if ! command -v cargo-vendor-filterer >/dev/null 2>&1; then
    echo "cargo-vendor-filterer is required to build < 5 MB CRAN source archives." >&2
    exit 1
  fi
  if cargo vendor-filterer --offline vendor > vendor-config.toml 2>/dev/null; then
    :
  else
    cargo vendor-filterer vendor > vendor-config.toml
  fi
  if [[ "$(tar --version)" == bsdtar* ]]; then
    COPYFILE_DISABLE=1 tar --no-xattrs --no-mac-metadata -cf - vendor | xz -9e > vendor.tar.xz
  else
    tar -cf - vendor | xz -9e > vendor.tar.xz
  fi
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
if [[ "$bytes" -ge 5000000 ]]; then
  echo "The package exceeds the 5 MB CRAN source archive threshold ($bytes bytes)." >&2
  exit 1
fi
