#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SPATIAL_RUST_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORKSPACE_DIR="$(cd "${SPATIAL_RUST_DIR}/.." && pwd)"

CHECK_ONLY=0
PKG_DIR=""
for arg in "$@"; do
  case "${arg}" in
    --check)
      CHECK_ONLY=1
      ;;
    -h|--help)
      echo "Usage: $0 [--check] [R_PACKAGE_DIR]"
      exit 0
      ;;
    *)
      if [[ -n "${PKG_DIR}" ]]; then
        echo "Unexpected extra argument: ${arg}" >&2
        exit 2
      fi
      PKG_DIR="${arg}"
      ;;
  esac
done

if [[ -z "${PKG_DIR}" ]]; then
  PKG_DIR="${WORKSPACE_DIR}/scpcR"
fi

DEST_RUST_DIR="${PKG_DIR}/src/rust"
mkdir -p "${DEST_RUST_DIR}"

sync_crate() {
  local crate_name="$1"
  local src="${SPATIAL_RUST_DIR}/crates/${crate_name}"
  local dst="${DEST_RUST_DIR}/${crate_name}"

  if [[ ! -d "${src}" ]]; then
    echo "Source crate directory missing: ${src}" >&2
    exit 1
  fi

  if [[ "${CHECK_ONLY}" -eq 0 ]]; then
    rm -rf "${dst}"
    mkdir -p "${dst}"
    (cd "${src}" && tar --exclude='target*' --exclude='.DS_Store' -cf - .) | (cd "${dst}" && tar -xf -)
  fi

  if [[ ! -d "${dst}" ]]; then
    echo "Destination crate directory missing: ${dst}" >&2
    exit 1
  fi

  diff -rq -x 'target*' -x '.DS_Store' "${src}" "${dst}"
}

sync_crate "spatial-core"
sync_crate "spatial-r"

echo "Verified byte-for-byte sync for spatial-core and spatial-r into ${DEST_RUST_DIR}"
