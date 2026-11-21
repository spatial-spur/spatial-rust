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
  if cargo vendor --offline --locked vendor > vendor-config.toml 2>/dev/null; then
    :
  else
    cargo vendor --locked vendor > vendor-config.toml
  fi

  # Prune unreachable dead-target crates (pulled only via the empty intersection
  # cfg(apple) -> sysctl -> cfg(linux) -> walkdir -> cfg(windows) -> winapi-util -> windows-sys)
  # and strip integration tests/benches/examples so CRAN source tarballs stay well under 5 MB.
  python3 - <<'PY'
import hashlib, json, os, re, shutil

vendor_dir = "vendor"
dead_crates = {
    "windows-sys",
    "windows-link",
    "winapi-util",
    "walkdir",
    "same-file",
    "hermit-abi",
}

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

for entry in sorted(os.listdir(vendor_dir)):
    cdir = os.path.join(vendor_dir, entry)
    if not os.path.isdir(cdir):
        continue
    chk_path = os.path.join(cdir, ".cargo-checksum.json")
    if entry in dead_crates and os.path.exists(chk_path):
        with open(chk_path, "r") as f:
            pkg_hash = json.load(f).get("package")
        keep = {}
        for fn in os.listdir(cdir):
            fp = os.path.join(cdir, fn)
            if (fn == "Cargo.toml" or re.match(r"^(LICENSE|COPYING|UNLICENSE)", fn, re.I)) and os.path.isfile(fp):
                with open(fp, "rb") as rf:
                    keep[fn] = rf.read()
        for fn in os.listdir(cdir):
            fp = os.path.join(cdir, fn)
            if os.path.isdir(fp):
                shutil.rmtree(fp)
            else:
                os.remove(fp)
        for fn, data in keep.items():
            with open(os.path.join(cdir, fn), "wb") as wf:
                wf.write(data)
        os.makedirs(os.path.join(cdir, "src"), exist_ok=True)
        with open(os.path.join(cdir, "src", "lib.rs"), "w") as lf:
            lf.write("//! Unused dead-target stub for CRAN source vendoring.\n")
        files = {}
        for root, _, fns in os.walk(cdir):
            for fn in sorted(fns):
                fp = os.path.join(root, fn)
                rel = os.path.relpath(fp, cdir).replace(os.sep, "/")
                files[rel] = sha256_file(fp)
        with open(chk_path, "w") as cf:
            json.dump({"files": files, "package": pkg_hash}, cf, separators=(",", ":"))
        continue

    changed = False
    for sub in ("tests", "benches", "examples"):
        sp = os.path.join(cdir, sub)
        if os.path.isdir(sp):
            shutil.rmtree(sp)
            changed = True
    if changed and os.path.exists(chk_path):
        with open(chk_path, "r") as f:
            data = json.load(f)
        data["files"] = {
            k: v for k, v in data.get("files", {}).items() if os.path.exists(os.path.join(cdir, k))
        }
        with open(chk_path, "w") as f:
            json.dump(data, f, separators=(",", ":"))
PY

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
