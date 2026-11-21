#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 2 ]]; then
  echo "Usage: $0 R_PACKAGE_DIR OUTPUT_DIR" >&2
  exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
orig_pwd="$(pwd)"
package_dir="$(cd "$1" && pwd)"
package="$(sed -n 's/^Package:[[:space:]]*//p' "$package_dir/DESCRIPTION" | head -n 1 | tr -d '\r')"
if [[ "$package" != scpcR && "$package" != spuR ]]; then
  echo "Expected an scpcR or spuR package (got Package: '${package}')." >&2
  exit 2
fi
bash "$script_dir/sync_r_crates.sh" --check "$package_dir"
staging="$(mktemp -d)"
trap 'rm -rf "$staging"' EXIT

# Stage working-tree sources without any artifacts from a previous local build.
mkdir "$staging/$package"
rsync -a --exclude=/.git --exclude='/src/rust/target*' --exclude=/src/rust/vendor \
  --exclude=/src/rust/.cargo --exclude=/src/rust/vendor.tar.xz \
  --exclude=/src/rust/vendor-config.toml --exclude=/src/rust/spatial-core/tests \
  --exclude='/*.tar.gz' --exclude='/*.Rcheck' --exclude='/src/*.o' --exclude='/src/*.so' \
  --exclude='/src/*.dylib' --exclude='/src/*.dll' --exclude='/src/.spatial_rust_mode' \
  --exclude='*.DS_Store' \
  "$package_dir/" "$staging/$package/"
mkdir -p "$2"
output_dir="$(cd "$orig_pwd" && cd "$2" && pwd)"

(
  cd "$staging/$package/src/rust"
  if cargo vendor --offline --locked vendor > vendor-config.toml 2>/dev/null; then
    :
  else
    cargo vendor --locked vendor > vendor-config.toml
  fi

  # Package the complete locked registry sources, including their original checksums.
  if [[ "$(tar --version)" == bsdtar* ]]; then
    if command -v xz >/dev/null 2>&1; then
      COPYFILE_DISABLE=1 tar --no-xattrs --no-mac-metadata -cf - vendor | xz -9e > vendor.tar.xz
    else
      COPYFILE_DISABLE=1 tar --no-xattrs --no-mac-metadata -cJf vendor.tar.xz vendor
    fi
  else
    if command -v xz >/dev/null 2>&1; then
      tar -cf - vendor | xz -9e > vendor.tar.xz
    else
      XZ_OPT=-9e tar -cJf vendor.tar.xz vendor
    fi
  fi
  rm -rf vendor
)
cd "$output_dir"
R CMD build --no-build-vignettes "$staging/$package"
