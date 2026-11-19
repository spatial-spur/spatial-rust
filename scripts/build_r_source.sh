#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 R_PACKAGE_DIR OUTPUT_DIR" >&2
  exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
package_dir="$(cd "$1" && pwd)"
package="$(sed -n 's/^Package:[[:space:]]*//p' "$package_dir/DESCRIPTION" | head -n 1 | tr -d '\r')"
if [[ "$package" != scpcR && "$package" != spuR ]]; then
  echo "Expected an scpcR or spuR package (got Package: '${package}')." >&2
  exit 2
fi
mkdir -p "$2"
output_dir="$(cd "$2" && pwd)"
staging="$(mktemp -d)"
trap 'rm -rf "$staging"' EXIT

# Stage working-tree sources without any artifacts from a previous local build.
mkdir "$staging/$package"
rsync -a --exclude=/.git --exclude=/src/rust/target --exclude=/src/rust/vendor \
  --exclude=/src/rust/.cargo --exclude=/src/rust/vendor.tar.xz \
  --exclude=/src/rust/vendor-config.toml --exclude=/src/rust/spatial-core/tests \
  --exclude='/*.tar.gz' --exclude='/*.Rcheck' --exclude='/src/*.o' --exclude='/src/*.so' \
  --exclude='/src/*.dylib' --exclude='/src/*.dll' --exclude='/src/.spatial_rust_mode' \
  --exclude='*.DS_Store' \
  "$package_dir/" "$staging/$package/"
bash "$script_dir/sync_r_crates.sh" "$staging/$package"

(
  cd "$staging/$package/src/rust"
  cargo vendor --locked vendor > vendor-config.toml
  COPYFILE_DISABLE=1 tar -cJf vendor.tar.xz vendor
  rm -rf vendor
)
cd "$output_dir"
R CMD build --no-build-vignettes "$staging/$package"
