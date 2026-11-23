"""Bind built artifacts to source and methods versions."""

import argparse
import hashlib
import json
import subprocess
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", nargs="+", type=Path)
    args = parser.parse_args()
    if git("status", "--porcelain"):
        parser.error("commit the source tree before recording release artifacts")
    version = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    spec = (ROOT / "spec" / "VERSION").read_text().strip()
    manifest = {
        "source_revision": git("rev-parse", "HEAD"),
        "distribution_version": version,
        "methods_version": spec,
        "artifacts": [],
    }
    for path in sorted(args.artifacts):
        if not path.is_file():
            parser.error(f"artifact {path} is not a file")
        with path.open("rb") as source:
            digest = hashlib.file_digest(source, "sha256").hexdigest()
        manifest["artifacts"].append(
            {"filename": path.name, "bytes": path.stat().st_size, "sha256": digest}
        )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
